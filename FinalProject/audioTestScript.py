import numpy as np
import griffinLimInverse as gl
import librosa
song = librosa.core.load('song.mp3',duration=5)
bins=4096
hops=1024
iters=20
spectrum = gl.STFT(song[0],bins,hops)
reconstructed = gl.griffinLimInverse(spectrum,bins,hops,iters)
librosa.output.write_wav("./reconstructed2.wav",reconstructed,22050)