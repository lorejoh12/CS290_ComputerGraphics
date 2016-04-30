import numpy as np
import griffinLimInverse as gl
import librosa
# This is a very simple script to test Griffin Lim transformation
song = librosa.core.load('song.mp3',duration=120)
bins=4096
hops=1024
iters=20
spectrum = gl.STFT(song[0],bins,hops)
reconstructed = gl.griffinLimInverse(spectrum,bins,hops,iters)
librosa.output.write_wav("./reconstructedSong20.wav",reconstructed,22050)