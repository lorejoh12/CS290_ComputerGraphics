import sys
import numpy as np 
import matplotlib.pyplot as pyplot
import math

THRESHOLD = 2.2204e-16


# Do the short-time Fourier transform on a signal using a half-sine window
# Inputs:
#   signal: The Nx1 input signal
#   window_size: The window size
#   hop_size: The hop size
# Returns:
#   matrix of the STFT magnitudes
def STFT(signal, window_size, hop_size):
    hops_per_win = window_size/hop_size
    # Make sure we're at an integer multiple of hops per window
    if hops_per_win - math.floor(hops_per_win) > 0:
        print "Warning: Window size is not integer multiple of hop size"
    # Create our window values
    win = []
    for i in range(window_size):
        win.append(math.sin(math.pi*i/(window_size-1)))
    num_windows = math.floor((len(signal)-window_size)/hop_size)+1
    # Create a spectrum of all the values for the window size over the number of windows
    spectrum = np.zeros((num_windows, window_size), dtype=np.complex64)
    for i in range(int(num_windows)):
        sampled_signal = [];
        for k in range(0,window_size):
            sampled_signal.append(signal[k+(i*hop_size)])
        fft_vals = np.fft.fft(np.multiply(sampled_signal, win))
        spectrum[i] = fft_vals
    # Second half of the spectrum is redundant for real signals
    col_index=0
    if W%2 == 0:
        # Even
        col_index = window_size/2+1
    else:
        # Odd
        col_index = (window_size-1)/2+1
    spectrum = spectrum[:,0:col_index]
    return spectrum

# Do the inverse short-time Fourier transform on a signal using a half-sine window
# Inputs:
#   spectrum: The spectrum output of the STFT
#   window_size: The window size
#   hop_size: The hop size
# Returns:
#   Nx1 signal
def iSTFT(spectrum, window_size, hop_size):
    #First, put back redundant STFT
    if W%2==0:
        #Even case
        spectrum2 = spectrum[:,1:np.shape(spectrum)[1]-1]
        spectrum = np.concatenate((spectrum,np.fliplr(np.conjugate(spectrum2))),axis=1)

    else:
        spectrum2 = spectrum[:,1:np.shape(spectrum)[1]]
        spectrum = np.concatenate((spectrum,np.fliplr(np.conjugate(spectrum2))),axis=1)
    # Figure out how long the reconstructed signal really is
    signal_length = window_size+hop_size*(np.shape(spectrum)[0]-1)
    reconstructed_signal = np.zeros(signal_length)
    # Setup the window
    hops_per_win = window_size/hop_size
    # Make sure we're at an integer multiple of hops per window
    if hops_per_win - math.floor(hops_per_win) > 0:
        print "Warning: Window size is not integer multiple of hop size"
    # Create our window values
    win = []
    for i in range(window_size):
        win.append(math.sin(math.pi*i/(window_size-1)))
    win = np.divide(win,(hops_per_win/2))
    # Do overlap/add synthesis
    for i in range(np.shape(spectrum)[0]):
        i1 = 1+(i*hop_size)
        i2 = i1 + window_size-1
        reconstructed_signal[i1-1:i2]= reconstructed_signal[i1-1:i2]+np.real(np.array(win).T*np.fft.ifft(spectrum[i,:]))


    return reconstructed_signal


# Take the inverse griffith lim transformation of a signal transformed using the STFT
# Inputs: 
#   spectrum: the magnitude spectrum from an STFT using a half sine window
#   window_size: the size of the window
#   hop_size: the size of the hop
#   n_iters: the number of times to apply this inversion
# Return: The reconstructed signal

def griffinLimInverse(spectrum, window_size, hop_size, num_iters):
    spec_temp = spectrum
    for i in range(num_iters):
        print "Iteration %d of %d" % (i+1,num_iters)
        spec_temp = STFT(iSTFT(spec_temp,window_size,hop_size), window_size, hop_size)
        norm = np.real(np.sqrt(spec_temp*np.conjugate(spec_temp)))
        norm[norm<THRESHOLD]=1
        spec_temp = np.multiply(spectrum,np.divide(spec_temp,norm))
    return iSTFT(spec_temp, window_size, hop_size)


