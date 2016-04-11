import numpy as np
import griffinLimInverse as gl
import librosa
song = librosa.core.load('song.mp3',duration=5)
spectrum = gl.STFT(song[0],100,5)
reconstructed = gl.griffinLimInverse(spectrum,100,5,50)
librosa.output.write_wav("./reconstructed.wav",reconstructed,22050)