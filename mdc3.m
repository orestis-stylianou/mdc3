% Multiscale Detrended Cross-Correlation Coefficient (based on 10.3389/fnins.2024.1422085)
% Author: Orestis Stylianou, orestisstylianou@rocketmail.com
% data: dataset for which mdc3 is to be estimated
% channels_number: number of channels/signals
% fmin: minimum frequency for mdc3 estimation
% fmax: maximum frequency for mdc3 estimation
% fstep: frequency increments from fmin to fmax for mdc3 estimation
% sampling_rate: sampling rate of the recorded signals
% trend: degree of the polynomial detrending
% directed: selection between directed(1) and undirected(0) variant of mdc3 estimation
function [avg_coefficient] = mdc3(data,channels_number,fmin,fmax,fstep,sampling_rate,trend,directed)
if ~islogical(directed)
     error('directed must be a logical value')
end
if size(data,1) ~= channels_number && size(data,2) ~= channels_number
    error('Channel number does not agree with data')
end
frequencies = fmin:fstep:fmax;
scales = round(sampling_rate./frequencies);
scales = unique(scales);
frequencies = sampling_rate./scales;
frequencies(frequencies < fmin |frequencies > fmax) = [];
scales = round(sampling_rate./frequencies);
scales_number = length(scales);
combinations = nchoosek(1:channels_number,2);
sequence = 1:channels_number^2;
sequence = reshape(sequence,[channels_number, channels_number])';
indexes = logical(triu(ones(size(sequence)), 1));
sequence = sort(sequence(indexes));
if size(data,1) == channels_number
        data = data';
end
coefficient = zeros(channels_number,channels_number,scales_number);
for s = 1:scales_number
    scale = scales(s);
    finish = 0;
    covariance = zeros(channels_number);
    variance = zeros(1,channels_number);
    for window = 1:floor(size(data,1)/scale)
        start = finish + 1; 
        finish = start + scale-1;
        tmp = data(start:finish,:);
        tmp = detrend(tmp,trend);
        if directed == 0
            covariance_values = cov(tmp);
        else
            covariance_values = zeros(channels_number);
            [c,lags] = xcov(tmp,'biased');
            zero_lag_index = find(lags == 0);
            leading_covariance = c(1:zero_lag_index-1,:);
            following_covariance = c(zero_lag_index+1:end,:);
            max_leading_covariance = max(leading_covariance);
            max_following_covariance = max(following_covariance);
            min_leading_covariance = min(leading_covariance);
            min_following_covariance = min(following_covariance);
            for combo = 1:size(combinations,1)
                 sequence_index = sequence(combo);
                 ch1 = combinations(combo,1);
                 ch2 = combinations(combo,2);
                 if abs(max_leading_covariance(sequence_index)) > abs(min_leading_covariance(sequence_index))
                    covariance_values(ch2,ch1) = max_leading_covariance(sequence_index);
                 elseif abs(max_leading_covariance(sequence_index)) < abs(min_leading_covariance(sequence_index))
                     covariance_values(ch2,ch1) = min_leading_covariance(sequence_index);
                 else
                    covariance_values(ch2,ch1) = 0; 
                 end
                 if abs(max_following_covariance(sequence_index)) > abs(min_following_covariance(sequence_index))
                     covariance_values(ch1,ch2) = max_following_covariance(sequence_index);
                 elseif abs(max_following_covariance(sequence_index)) < abs(min_following_covariance(sequence_index))
                     covariance_values(ch1,ch2) = min_following_covariance(sequence_index);
                 else
                     covariance_values(ch1,ch2) = 0;
                 end
            end
        end
        for combo = 1:size(combinations,1)
            ch1 = combinations(combo,1);
            ch2 = combinations(combo,2);
            covariance(ch1,ch2) = covariance_values(ch1,ch2) + covariance(ch1,ch2);
            covariance(ch2,ch1) = covariance_values(ch2,ch1) + covariance(ch2,ch1);
        end
        for channel = 1:channels_number
            variance(channel) = var(tmp(:,channel)) + variance(channel);
        end
    end
    covariance = covariance/floor(size(data,1)/scale);
    variance = variance/floor(size(data,1)/scale);
    for combo = 1:size(combinations,1)
        ch1 = combinations(combo,1);
        ch2 = combinations(combo,2);
        coefficient(ch1,ch2,s)  = covariance(ch1,ch2)/(sqrt(variance(ch1)*variance(ch2)));
        coefficient(ch2,ch1,s) = covariance(ch2,ch1)/(sqrt(variance(ch1)*variance(ch2)));
    end
end
avg_coefficient = zeros(channels_number);
for combo = 1:size(combinations,1)
    ch1 = combinations(combo,1);
    ch2 = combinations(combo,2);
    signal1 = detrend(data(:,ch1),trend);
    signal2 = detrend(data(:,ch2),trend);
    spectrum = cpsd(signal1,signal2,[],[],frequencies,sampling_rate);
    spectrum = abs(spectrum);
    weights = spectrum/(sum(spectrum));
    if directed == 0
        connection_coefficients = reshape(coefficient(ch1,ch2,:),[1 scales_number]);
        connection_coefficients = atanh(connection_coefficients);
        avg_coefficient(ch1,ch2) = tanh(sum(connection_coefficients.*weights));
        avg_coefficient(ch2,ch1) = tanh(sum(connection_coefficients.*weights));  
    else
        leading_connection_coefficients = reshape(coefficient(ch2,ch1,:),[1 scales_number]);
        following_connection_coefficients = reshape(coefficient(ch1,ch2,:),[1 scales_number]);
        leading_connection_coefficients = atanh(leading_connection_coefficients);
        following_connection_coefficients = atanh(following_connection_coefficients);
        avg_coefficient(ch2,ch1) = tanh(sum(leading_connection_coefficients.*weights));  
        avg_coefficient(ch1,ch2) = tanh(sum(following_connection_coefficients.*weights));
    end
end
end
