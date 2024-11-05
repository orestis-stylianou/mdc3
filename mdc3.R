# Multiscale Detrended Cross-Correlation Coefficient (based on 10.3389/fnins.2024.1422085)
# Author: Orestis Stylianou, orestisstylianou@rocketmail.com
# data: dataset for which mdc3 is to be estimated
# channels_number: number of channels/signals
# fmin: minimum frequency for mdc3 estimation
# fmax: maximum frequency for mdc3 estimation
# fstep: frequency increments from fmin to fmax for mdc3 estimation
# sampling_rate: sampling rate of the recorded signals
# trend: degree of the polynomial detrending
# directed: selection between directed(TRUE) and undirected(FALSE) variant of mdc3 estimation
mdc3 <- function (data,channels_number,fmin,fmax,fstep,sampling_rate,trend,directed){
    library(utils)
    library(IRISSeismic)
    library(pracma)
    library(gsignal)
    library(matrixStats)
    if (!is.logical(directed)) {
        stop("directed must be a logical value")
    }
    if (dim(data)[1] != channels_number && dim(data)[2] != channels_number){
        stop('Channel number does not agree with data')
    }
    frequencies <- seq(fmin,fmax,fstep)
    scales <- round(sampling_rate/frequencies)
    scales <- unique(scales)
    frequencies <- sampling_rate/scales
    index_1 <- which(frequencies >= fmin)
    index_1 <- index_1[1]
    index_2 <- which(frequencies <= fmax)
    index_2 <- index_2[length(index_2)]
    frequencies <-  frequencies[index_1:index_2]
    scales <- round(sampling_rate/frequencies)
    scales_number <- length(scales)
    combinations <- combn(1:channels_number,2)
    sequence <- seq(1,channels_number^2)
    sequence <- matrix(sequence, ncol = channels_number, nrow = channels_number,byrow = TRUE)
    sequence <- sort(sequence[upper.tri(sequence)])
    if (dim(data)[1] == channels_number){
        data <- t(data)
    }
    tmp <- zeros(channels_number,channels_number)
    coefficient <- array(c(tmp, tmp), dim = c(channels_number, channels_number, scales_number))
    for (s in 1:scales_number){
        scale <- scales[s]
        finish <- 0
        covariance <- zeros(channels_number)
        variance <- zeros(1,channels_number)
        for (window in 1:floor(dim(data)[1]/scale)){
            start <- finish + 1
            finish <- start + scale-1
            tmp <- data[start:finish,]
            model <- lm(tmp ~ stats::poly(1:nrow(tmp), degree = trend, raw = TRUE))
            fitted_values <- predict(model)
            detrended <- tmp - fitted_values
            tmp <- detrended
            if (directed == FALSE){
                covariance_values <- cov(tmp)
            } else {
                covariance_values <- zeros(channels_number)
                xcov_results <- xcov(tmp,scale = 'biased')
                c <- xcov_results$C
                lags <- xcov_results$lags
                zero_lag_index <- which(lags == 0)
                leading_covariance <- c[1:zero_lag_index-1,]
                following_covariance <- c[(zero_lag_index+1):dim(c)[1],]
                max_leading_covariance <- colMaxs(leading_covariance)
                max_following_covariance <- colMaxs(following_covariance)
                min_leading_covariance <- colMins(leading_covariance)
                min_following_covariance <- colMins(following_covariance)
                for (combo in 1:dim(combinations)[2]){
                    sequence_index <- sequence[combo]
                    ch1 <- combinations[1,combo]
                    ch2 <- combinations[2,combo]
                    if (abs(max_leading_covariance[sequence_index]) > abs(min_leading_covariance[sequence_index])){
                        covariance_values[ch2,ch1] <- max_leading_covariance[sequence_index]
                    } else if (abs(max_leading_covariance[sequence_index]) < abs(min_leading_covariance[sequence_index])){
                        covariance_values[ch2,ch1] <- min_leading_covariance[sequence_index]
                    } else {
                        covariance_values[ch2,ch1] <- 0
                    }
                    if (abs(max_following_covariance[sequence_index]) > abs(min_following_covariance[sequence_index])){
                        covariance_values[ch1,ch2] <- max_following_covariance[sequence_index]
                    } else if (abs(max_following_covariance[sequence_index]) < abs(min_following_covariance[sequence_index])) {
                        covariance_values[ch1,ch2] <- min_following_covariance[sequence_index]
                    } else {
                        covariance_values[ch1,ch2] <- 0
                    }
                }
            }
            for (combo in 1:dim(combinations)[2]){
                ch1 <- combinations[1,combo]
                ch2 <- combinations[2,combo]
                covariance[ch1,ch2] <- covariance_values[ch1,ch2] + covariance[ch1,ch2]
                covariance[ch2,ch1] <- covariance_values[ch2,ch1] + covariance[ch2,ch1]
            }
            variance <- diag(cov(tmp)) + variance
        }
        covariance <- covariance/floor(dim(data)[1]/scale)
        variance <- variance/floor(dim(data)[1]/scale)
        for (combo in 1:dim(combinations)[2]){
            ch1 <- combinations[1,combo]
            ch2 <- combinations[2,combo]
            coefficient[ch1,ch2,s] <- covariance[ch1,ch2]/sqrt(variance[ch1]*variance[ch2])
            coefficient[ch2,ch1,s] <- covariance[ch2,ch1]/sqrt(variance[ch1]*variance[ch2])
        }
    }
    avg_coefficient <- zeros(channels_number)
    for (combo in 1:dim(combinations)[2]){
        ch1 <- combinations[1,combo]
        ch2 <- combinations[2,combo]
        signal1 <- data[,ch1]
        signal2 <- data[,ch2]
        fit_signal1 <- polyfit(1:length(signal1),signal1,trend)
        fit_signal1 <- polyval(fit_signal1,1:length(signal1))
        signal1 <- signal1 - fit_signal1
        fit_signal2 <- polyfit(1:length(signal2),signal2,trend)
        fit_signal2 <- polyval(fit_signal2,1:length(signal2))
        signal2 <- signal2 - fit_signal2
        signal1 <- ts(signal1,frequency=sampling_rate)
        signal2 <- ts(signal2,frequency=sampling_rate)
        signals <- ts.union(signal1,signal2)
        spectrum <- crossSpectrum(signals, detrend = FALSE)
        indexes <- zeros(length(frequencies),1)
        for (freq in 1:length(frequencies)){
            target_frequency <- frequencies[freq]
            index <- which(spectrum$freq == target_frequency)
            if (isempty(index)){
                index_1 <- which(spectrum$freq > target_frequency)
                index_1 <- index_1[1]
                diff_1 <- abs(spectrum$freq[index_1] - target_frequency)
                index_2 <- which(spectrum$freq < target_frequency)
                index_2 <- index_2[length(index_2)]
                diff_2 <- abs(spectrum$freq[index_2] - target_frequency)
                if (diff_1 < diff_2){
                    index <- index_1
                } else {
                    index <- index_2
                }
            }
            indexes[freq] <- index
        }
        spectrum <- spectrum$Pxy[indexes]
        spectrum <- abs(spectrum)
        weights <- spectrum/(sum(spectrum))
        if (directed == FALSE){
            connection_coefficients <- coefficient[ch1,ch2,]
            connection_coefficients <- atanh(connection_coefficients)
            avg_coefficient[ch1,ch2] <- tanh(sum(connection_coefficients*weights))
            avg_coefficient[ch2,ch1] <- avg_coefficient[ch1,ch2]
        } else {
            leading_connection_coefficients <- coefficient[ch2,ch1,]
            following_connection_coefficients <- coefficient[ch1,ch2,]
            leading_connection_coefficients <- atanh(leading_connection_coefficients)
            following_connection_coefficients <- atanh(following_connection_coefficients)
            avg_coefficient[ch2,ch1] <- tanh(sum(leading_connection_coefficients*weights))
            avg_coefficient[ch1,ch2] <- tanh(sum(following_connection_coefficients*weights))
        }
    }
    return(avg_coefficient)
}