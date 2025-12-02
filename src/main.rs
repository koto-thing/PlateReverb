use hound;

use plate_reverb::{PlateReverb, PhysicalParams};

fn main() {
    // 設定
    let sample_rate = 44100.0;
    let params = PhysicalParams {
        ly: 0.5,
        aspect: 1.0,
        h: 0.005,
        t0: 1500.0,
        rho: 7850.0,
        e: 2.0e11,
        t60_dc: 4.0,
        t60_f1: 0.05,
        op_x: 0.3,
        op_y: 0.3,
    };

    // インスタンスを作成
    let mut reverb = PlateReverb::new(params, sample_rate);

    // インパルスを与える
    reverb.add_impulse(1.0);

    // WAVファイルの設定
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate: sample_rate as u32,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };

    let mut writer = hound::WavWriter::create("output.wav", spec).unwrap();

    // シミュレーションループ
    println!("Simulating...");
    for _ in 0..(sample_rate * 5.0) as usize {
        let s = reverb.process_step();
        let output = (s * 100.0).clamp(-1.0, 1.0); // ちょっと音量を上げる
        writer.write_sample((output * i16::MAX as f64) as i16).unwrap();
    }

    println!("Done!");
}
