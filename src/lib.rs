const NU: f64 = 0.3; // ポアソン比

#[derive(Debug, Clone)]
pub struct PhysicalParams {
    pub ly: f64,       // 板の高さ (m)
    pub aspect: f64, // 板のアスペクト比
    pub h: f64,        // 板の厚さ (m)
    pub t0: f64,       // 板の張力 (N/m)
    pub rho: f64,      // 板の密度 (kg/m^3)
    pub e: f64,        // ヤング率 (Pa)

    pub t60_dc: f64,   // 低域の減衰時間 (s)
    pub t60_f1: f64,   // 高域の減衰時間 (s)

    pub op_x: f64,     // 出力位置 X (0.0 - 1.0)
    pub op_y: f64,     // 出力位置 Y (0.0 - 1.0)
}

pub struct PlateReverb {
    // パラメータ
    width: usize,
    height: usize,

    // 状態バッファ
    u_current: Vec<f64>,
    u_prev: Vec<f64>,
    u_next: Vec<f64>,

    // 物理パラメータ
    mu: f64,    // 剛性項の係数
    lambda: f64, // 張力項の係数

    // 減衰
    sigma0: f64, // 距離に依存しない減衰 (空気抵抗)
    sigma1: f64, // 距離(曲率)に依存する減衰 (粘性抵抗)
    k: f64,      // 時間ステップ

    // ピックアップ位置
    pickup_index: usize,
}

impl PlateReverb {
    pub fn new(params: PhysicalParams, sample_rate: f64) -> Self {
        let nu = 0.3;              // ポアソン比
        let k = 1.0 / sample_rate; // 時間ステップ(サンプル周期)

        // 曲げ剛性を計算
        let d = (params.e * params.h.powi(3)) / (12.0 * (1.0 - nu * nu));

        // 安定性条件からグリッド間隔を計算する
        let kappa = (d / (params.rho * params.h)).sqrt();

        // 安定性のための最小グリッド間隔を計算する
        // 安全のためにちょっと数値を盛る
        let max_mu :f64 = 0.05;
        let ds_limit_sq = (k * kappa) / max_mu.sqrt();
        let ds_limit = ds_limit_sq.sqrt();

        // グリッドの数を計算する
        let h_grid_count = (params.ly / ds_limit).floor() as usize;
        let w_grid_count = ((params.ly * params.aspect) / ds_limit).floor() as usize;
        let h_grid_count = h_grid_count.max(5);
        let w_grid_count = w_grid_count.max(5);

        // 実際のグリッド間隔を計算する
        let ds = params.ly / h_grid_count as f64;
        let ds_sq = ds * ds;

        // 係数を計算する
        let mu = (k * kappa / ds_sq).powi(2);

        // 張力係数lambdaを計算する
        // 波の速度は c = sqrt(T0 / (rho * h))
        let c_sq = params.t0 / (params.rho * params.h);
        let lambda = (c_sq * k * k) / ds_sq;

        // 安定性の確認
        println!("Grid: {}x{}, mu: {:.5} (stable if <= 0.25)", w_grid_count, h_grid_count, mu);

        // 減衰係数を計算する
        // 低域の減衰
        let sigma0 = 6.907 / params.t60_dc;

        // 高域の減衰
        let omega_high = std::f64::consts::PI * sample_rate;
        let desired_decay_high = 6.907 / params.t60_f1;
        let sigma1 = if desired_decay_high > sigma0 {
            (desired_decay_high - sigma0) / (omega_high * omega_high)
        } else {
            0.0
        };

        // ピックアップ位置を計算する
        let pickup_x = (params.op_x * w_grid_count as f64) as usize;
        let pickup_y = (params.op_y * h_grid_count as f64) as usize;
        let pickup_index = pickup_y * w_grid_count + pickup_x;

        let size = w_grid_count * h_grid_count;

        Self {
            width: w_grid_count,
            height: h_grid_count,
            u_current: vec![0.0; size],
            u_prev: vec![0.0; size],
            u_next: vec![0.0; size],
            mu,
            lambda,
            sigma0,
            sigma1,
            k,
            pickup_index: if pickup_index < size { pickup_index } else { size / 2 },
        }
    }

    /// インパルスを入力する
    pub fn add_impulse(&mut self, amount: f64) {
        // DEBUG: 板の少し中心からずれたところをたたく
        let x = (self.width as f64 * 0.3) as usize;
        let y = (self.height as f64 * 0.3) as usize;
        let index = y * self.width + x;
        if index < self.u_current.len() {
            self.u_current[index] += amount;
            self.u_prev[index] += amount;
        }
    }

    pub fn process_step(&mut self) -> f64 {
        let width = self.width;
        let height = self.height;

        let k = self.k;
        let sigma0 = self.sigma0;
        let sigma1 = self.sigma1;

        let denominator = 1.0 + sigma0 * k;
        let inv_denominator = 1.0 / denominator;

        // 境界は計算せず、固定端とする
        for y in 2..(height - 2) {
            for x in 2..(width - 2) {
                let index = y * width + x;
                let u_cnt = self.u_current[index];
                let u_prv = self.u_prev[index];

                // 5点ステンシル法によるラプラシアンの離散化
                // 現在の時間ステップ
                let s1_curr = self.u_current[index - 1]
                    + self.u_current[index + 1]
                    + self.u_current[index - width]
                    + self.u_current[index + width];
                let lap_curr = s1_curr - 4.0 * u_cnt;

                // 過去の時間ステップ
                let s1_prev = self.u_prev[index - 1]
                    + self.u_prev[index + 1]
                    + self.u_prev[index - width]
                    + self.u_prev[index + width];
                let lap_prev = s1_prev - 4.0 * u_prv;

                // 13点ステンシル法によるラプラシアンの離散化
                let s2 = self.u_current[index - width - 1]
                    + self.u_current[index - width + 1]
                    + self.u_current[index + width - 1]
                    + self.u_current[index + width + 1];

                let s3 = self.u_current[index - 2]
                    + self.u_current[index + 2]
                    + self.u_current[index - 2 * width]
                    + self.u_current[index + 2 * width];

                // ビハーモニック演算子の離散化
                let biharmonic = 20.0 * u_cnt - 8.0 * s1_curr + 2.0 * s2 + 1.0 * s3;

                // 時間発展の方程式
                // u[n + 1] = 2u[n] - u[n - 1] + mu * ∇^4 u[n]
                // +で減衰をかける
                // 2.0 * center - self.u_prev[index]は慣性の法則
                // self.mu * biharmonicは剛性による復元力(曲がった板が元に戻ろうとする力)
                let inertia = 2.0 * u_cnt - u_prv * (1.0 - sigma0 * k);
                let forces = -self.mu * biharmonic + self.lambda * lap_curr;
                let viscosity = 2.0 * sigma1 * k * (lap_curr - lap_prev);

                self.u_next[index] = (inertia + forces + viscosity) * inv_denominator;
            }
        }

        // バッファ更新
        self.u_prev.copy_from_slice(&self.u_current);
        self.u_current.copy_from_slice(&self.u_next);

        self.u_current[self.pickup_index]
    }
}