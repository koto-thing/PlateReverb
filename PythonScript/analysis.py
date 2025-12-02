import numpy as np
import matplotlib.pyplot as plt
import librosa.display

# Wavを読み込む
y, sr = librosa.load('output.wav', sr=None)

# フーリエ変換を行う
n_fft = 2048
hop_length = 128
D = librosa.stft(y, n_fft=n_fft, hop_length=hop_length)
S_db = librosa.amplitude_to_db(np.abs(D), ref=np.max)

# スペクトログラムを表示する
plt.figure(figsize=(10, 6))
librosa.display.specshow(S_db, sr=sr, hop_length=hop_length, x_axis='time', y_axis='linear')
plt.colorbar(format='%+2.0f dB')
plt.title('Spectrogram of Plate Reverb')
plt.show()