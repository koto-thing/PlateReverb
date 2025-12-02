import numpy as np
import soundfile as sf
from scipy import signal

def apply_reverb(dry_file, ir_file, output_file, mix = 0.5):
    print(f"Loading {dry_file} and {ir_file}...")

    # オーディオファイルを読み込む
    dry_data, sr_dry = sf.read(dry_file)
    ir_data, sr_ir = sf.read(ir_file)

    # サンプリングレートが異なる場合はエラーを出す
    if sr_dry != sr_ir:
        raise ValueError(f"Sample rates do not match! Dry: {sr_dry}, IR: {sr_ir}")

    # ステレオに対応にしとく
    if dry_data.ndim == 2:
        print("Processing stereo dry signal...")
        wet_data_left = signal.fftconvolve(dry_data[:, 0], ir_data, mode='full')
        wet_data_right = signal.fftconvolve(dry_data[:, 1], ir_data, mode='full')

        # 長さをそろえる
        len_result = len(wet_data_left)
        wet_data = np.vstack((wet_data_left, wet_data_right)).T
    else:
        print("Processing mono dry signal...")
        wet_data = signal.fftconvolve(dry_data, ir_data, mode='full')

    # ドライとウェットのミックス
    pad_length = len(wet_data) - len(dry_data)
    if dry_data.ndim == 2:
        dry_padded = np.pad(dry_data, ((0, pad_length), (0, 0)))
    else:
        dry_padded = np.pad(dry_data, (0, pad_length))

    # ノーマライズ
    wet_data = wet_data / np.max(np.abs(wet_data))

    # ミックス
    output_data = (dry_padded * (1.0 - mix)) + (wet_data * mix)

    # クリッピング
    max_val = np.max(np.abs(output_data))
    if max_val > 1.0:
        output_data /= max_val

    print(f"Saving to {output_file}...")
    sf.write(output_file, output_data, sr_dry)
    print("Done!")

if __name__ == "__main__":
    apply_reverb("input.wav", "ir.wav", "output_reverb.wav", mix=0.5)