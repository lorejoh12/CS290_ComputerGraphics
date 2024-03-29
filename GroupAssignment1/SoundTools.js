//Do fast fft-based convolution
//It is assumed that the first buffer contains only a single channel and the second channel has more than onechannels
function doConvolution(buffer1, buffer2, numChannels) {
    //Allocate space for the output buffer and copy over
    var len1 = buffer1.getChannelData(0).length;
    var len2 = buffer2.getChannelData(0).length;
    //Zeropad to nearest power of 2 above N+M+1
    var NPadded = Math.pow(2, Math.ceil(Math.log(len1+len2+1)/Math.log(2)));
        
    convbuffer = context.createBuffer(numChannels, NPadded, globalFs);

    for(var channel=0; channel<numChannels; channel++){
        var samples1 = buffer1.getChannelData(0);
        M = samples1.length;
        var samples2 = buffer2.getChannelData(channel);
        N = samples2.length;
        
        //Zeropad to nearest power of 2 above N+M+1
        var NPadded = Math.pow(2, Math.ceil(Math.log(N+M+1)/Math.log(2)));
        
        //Zeropad signal 1
        var X = new Float32Array(NPadded);
        for (var i = 0; i < NPadded; i++) {
            if (i < M) {
                X[i] = samples1[i];
            }
            else {
                X[i] = 0;
            }
        }
        //Zeropad signal 2
        var Y = new Float32Array(NPadded);
        for (var i = 0; i < NPadded; i++) {
            if (i < N) {
                Y[i] = samples2[i];
            }
            else {
                Y[i] = 0;
            }
        }
        //For now, assume both sounds are sampled at 44100hz
        var fftX = new FFT(NPadded, NPadded);
        var fftY = new FFT(NPadded, NPadded);
        fftX.forward(X);
        fftY.forward(Y);
        //Multiply both in the frequency domain
        var real = new Float32Array(NPadded);
        var imag = new Float32Array(NPadded);
        for (var i = 0; i < NPadded; i++) {
            real[i] = fftX.real[i]*fftY.real[i] - fftX.imag[i]*fftY.imag[i];
            imag[i] = fftX.real[i]*fftY.imag[i] + fftX.imag[i]*fftY.real[i];
        }
        var fftRes = new FFT(NPadded, NPadded);
        var res = fftRes.inverse(real, imag);
        var maxAbs = 0.0;
        //Normalize output to prevent clipping
        for (var i = 0; i < res.length; i++) {
            if (Math.abs(res[i]) > maxAbs) {
                maxAbs = Math.abs(res[i]);
            }
        }
        //Copy over the output
        var convsamples = convbuffer.getChannelData(channel);
        for (var i = 0; i < res.length; i++) {
            convsamples[i] = res[i];
        }
    }
}
