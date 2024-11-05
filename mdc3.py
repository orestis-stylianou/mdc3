# Multiscale Detrended Cross-Correlation Coefficient (based on 10.3389/fnins.2024.1422085)
# Author: Orestis Stylianou, orestisstylianou@rocketmail.com
# data: dataset for which mdc3 is to be estimated
# channels_number: number of channels/signals
# fmin: minimum frequency for mdc3 estimation
# fmax: maximum frequency for mdc3 estimation
# fstep: frequency increments from fmin to fmax for mdc3 estimation
# sampling_rate: sampling rate of the recorded signals
# trend: degree of the polynomial detrending
# directed: selection between directed(True) and undirected(False) variant of mdc3 estimation
def mdc3(data, channels_number, fmin, fmax, fstep, sampling_rate, trend, directed):
    import numpy
    import pandas
    import itertools
    import obspy.signal.detrend
    import copy
    import scipy
    import statsmodels.tsa.tsatools as tools
    def nchoosek(startnum, endnum, step=1, n=1):
        c = []
        for i in itertools.combinations(range(startnum,endnum+1,step),n):
            c.append(list(i))
        return c
    def crosscov(signal1, signal2):
       lags = numpy.arange(-signal1.size+1,signal1.size)
       c = numpy.correlate(signal1,signal2,mode = "full") 
       c = c/len(signal1)
       return c, lags
    if not isinstance(directed, bool):
        print("directed must be a logical value")
        return
    if data.shape[0] != channels_number and data.shape[1] != channels_number:
        print("Channel number does not agree with data")
        return
    frequencies = numpy.arange(fmin,fmax+fstep,fstep)
    scales = numpy.round(sampling_rate/frequencies)
    scales = numpy.unique(scales)
    frequencies = sampling_rate/scales
    index_1 = numpy.where(frequencies <= fmax)
    index_1 = index_1[0][0]
    index_2 = numpy.where(frequencies >= fmin)
    index_2 = index_2[0][-1]
    frequencies = frequencies[index_1:index_2+1]
    scales = numpy.round(sampling_rate/frequencies)
    scales_number = len(scales)
    combinations = nchoosek(0, channels_number-1, step=1, n=2)
    if (data.shape[0] == channels_number):
            data = numpy.transpose(data)
    coefficient = numpy.zeros((channels_number,channels_number,scales_number))
    for s in numpy.arange(0,scales_number):
        scale = scales[s]
        finish = 0
        covariance = numpy.zeros((channels_number,channels_number))
        variance = numpy.zeros((1,channels_number))
        for window in numpy.arange(0,numpy.floor(data.shape[0]/scale)):
            start = int(finish)
            finish = int(start + scale)
            tmp = pandas.DataFrame.to_numpy(data.iloc[start:finish,])            
            detrended =  tools.detrend(tmp, order=trend, axis=0)
            if directed == False:
                covariance_values = numpy.cov(detrended,rowvar=False)
            elif directed == True:
                covariance_values = numpy.zeros((channels_number,channels_number))
                for combo in numpy.arange(0,len(combinations)):
                    ch1 = combinations[combo][0]
                    ch2 = combinations[combo][1]
                    signal1 = copy.deepcopy(detrended[:,ch1])
                    signal2 = copy.deepcopy(detrended[:,ch2])
                    [c,lags] = crosscov(signal1, signal2)
                    zero_lag_index = numpy.where(lags == 0)[0][0]
                    leading_covariance = c[0:zero_lag_index]
                    following_covariance = c[zero_lag_index+1:]
                    max_leading_covariance = numpy.max(leading_covariance)
                    max_following_covariance = numpy.max(following_covariance)
                    min_leading_covariance = numpy.min(leading_covariance)
                    min_following_covariance = numpy.min(following_covariance)
                    if abs(max_leading_covariance) > abs(min_leading_covariance):
                        covariance_values[ch2,ch1] = max_leading_covariance
                    elif abs(max_leading_covariance) < abs(min_leading_covariance):
                        covariance_values[ch2,ch1] = min_leading_covariance
                    else:
                        covariance_values[ch2,ch1] = 0
                    if abs(max_following_covariance) > abs(min_following_covariance):
                        covariance_values[ch1,ch2] = max_following_covariance
                    elif abs(max_following_covariance) < abs(min_following_covariance):
                        covariance_values[ch1,ch2] = min_following_covariance
                    else:
                        covariance_values[ch1,ch2] = 0
            for combo in numpy.arange(0,len(combinations)):
                ch1 = combinations[combo][0]
                ch2 = combinations[combo][1]
                covariance[ch1,ch2] = covariance_values[ch1,ch2] + covariance[ch1,ch2]
                covariance[ch2,ch1] = covariance_values[ch2,ch1] + covariance[ch2,ch1]
            window_variance = numpy.cov(detrended,rowvar=False)
            window_variance = numpy.diagonal(window_variance)
            variance = window_variance + variance
        covariance = covariance/numpy.floor(data.shape[0]/scale)
        variance = variance/numpy.floor(data.shape[0]/scale)
        for combo in numpy.arange(0,len(combinations)):
            ch1 = combinations[combo][0]
            ch2 = combinations[combo][1]
            coefficient[ch1,ch2,s] = covariance[ch1,ch2]/numpy.sqrt(variance[0,ch1]*variance[0,ch2])
            coefficient[ch2,ch1,s] = covariance[ch2,ch1]/numpy.sqrt(variance[0,ch1]*variance[0,ch2])
    avg_coefficient = numpy.zeros((channels_number,channels_number))
    for combo in numpy.arange(0,len(combinations)):
        ch1 = combinations[combo][0]
        ch2 = combinations[combo][1]
        signal1 = copy.deepcopy(data.iloc[:,ch1])
        signal2 = copy.deepcopy(data.iloc[:,ch2])
        obspy.signal.detrend.polynomial(signal1, order=trend, plot=False)
        obspy.signal.detrend.polynomial(signal2, order=trend, plot=False)
        p = numpy.ceil(numpy.log2(signal1.size))
        p = 2**p
        nfft = max(256,p)
        [freq,spectrum] = scipy.signal.csd(signal1, signal2, fs=sampling_rate, window='hamming', nperseg=signal1.size/8, noverlap=signal1.size/16, detrend = False,
                                           nfft=nfft, scaling='spectrum',average = 'median')
        spectrum = numpy.absolute(spectrum)
        indexes = numpy.zeros((len(frequencies)))
        for i in range(0,len(frequencies)):
            target_frequency = frequencies[i]
            index = numpy.where(freq == target_frequency)
            if len(index[0]) == 0:
                index_1 = numpy.where(freq > target_frequency)
                index_1 = index_1[0][0]
                diff_1 = abs(freq[index_1] - target_frequency)
                index_2 = numpy.where(freq < target_frequency)
                index_2 = index_2[0][-1]
                diff_2 = abs(freq[index_2] - target_frequency)
                if diff_1 < diff_2:
                    index = index_1
                else:
                    index = index_2
            else:
                index = index[0][0]
            indexes[i] = index 
        indexes = numpy.int_(indexes)
        spectrum = spectrum[indexes]
        weights = spectrum/numpy.sum(spectrum)
        if directed == False:
            connection_coefficients = coefficient[ch1,ch2,:]
            connection_coefficients = numpy.arctanh(connection_coefficients)
            avg_coefficient[ch1,ch2] = numpy.tanh(numpy.sum(connection_coefficients*weights))
            avg_coefficient[ch2,ch1] = numpy.tanh(numpy.sum(connection_coefficients*weights))
        elif directed == True:
            leading_connection_coefficients = coefficient[ch2,ch1,:]
            following_connection_coefficients = coefficient[ch1,ch2,:]
            leading_connection_coefficients = numpy.arctanh(leading_connection_coefficients)
            following_connection_coefficients = numpy.arctanh(following_connection_coefficients)
            avg_coefficient[ch2,ch1] = numpy.tanh(numpy.sum(leading_connection_coefficients*weights)) 
            avg_coefficient[ch1,ch2] = numpy.tanh(numpy.sum(following_connection_coefficients*weights))
        
    return avg_coefficient
    
    
    
    