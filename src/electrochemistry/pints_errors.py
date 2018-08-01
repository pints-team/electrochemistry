import pints

class FourierTransformedError(pints.ProblemErrorMeasure):
    def __init__(self, problem, omega, harmonics):
    """
    Args:

    problem: a pints problem
    omega (Hz): the frequency
    harmonics: list of the harmonics to include (0 = dc, 1 = 1st harmonic, ...)
    """
        super(FourierTransformedError, self).__init__(problem)
        signal_f = omega/(2*pi)
        #i = np.array([0,1,2,3,4,5, 6, 7, 8])
        #w = np.array([1,1,1,1,1,10,10,10,1])
        self._harmonics = harmonics
        hs = harmonics*omega
        bandwidth = omega/5.0
        #bandwidth = 2.0

        fs = 1.0/(t[1]-t[0])
        lowcut = hs - bandwidth/2
        highcut = hs + bandwidth/2
        lowcut[i==0] = 0
        highcut[i==0] = bandwidth/2

        all_taps = [bandpass(l,h,fs) for l,h in zip(lowcut,highcut)]

        #F_I = fft.rfft(np.hanning(len(I))*I)
        F_I = fft.rfft(I)
        freq_responses = [np.abs(signal.freqz(taps,worN=len(F_I))[1]) for taps in all_taps]

        if include_noise:
            weighted_freq_response = np.ones(len(F_I),dtype=complex)/linalg.norm(F_I)
        else:
            weighted_freq_response = np.zeros(len(F_I),dtype=complex)
        if rescale_harmonics:
            for weight,fr in zip(w,freq_responses):
                weighted_freq_response += weight*fr/linalg.norm(fr*F_I)
        else:
            for weight,fr in zip(w,freq_responses):
                weighted_freq_response += weight*fr

        #weighted_freq_response = weighted_freq_response / (np.sum(w)+1)

        I_filtered = F_I*weighted_freq_response

        return weighted_freq_response,I_filtered

    def __call__(self,x):
        return np.sum(((np.sum((self._problem.evaluate(x) - self._values)**2,
                               axis=0) * self._weights) * self._ninv),
                      axis=0)


