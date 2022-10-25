function [ amp_hf ] = get_amplitude(lfp, f_hf, pac_filter_bw_hf, ...
    pac_filter_type, pac_filter_order,fs )
%GET_AMPLITUDE_AND_PHASE Summary of this function goes here
%   Detailed explanation goes here

amp_hf = zeros(size(lfp));
f1_hz = (f_hf-pac_filter_bw_hf/2) ;
f2_hz = (f_hf+pac_filter_bw_hf/2) ;

if f1_hz < 1
    
    f1_hz = 1;
    pac_filter_bw_hf2 = f_hf - 1;
    f2_hz = (f_hf+pac_filter_bw_hf2/2);
end

if f2_hz >= fs/2;
    f2_hz = (fs/2)*0.95;
end

f1 = 2*f1_hz/fs;
f2 = 2*f2_hz/fs;

if strcmp(pac_filter_type, 'butter')
    
    [b_hf, a_hf] = butter(pac_filter_order,[f1 f2]);
    for k =1 : size(lfp,1)
        x_hf = filtfilt(b_hf,a_hf,lfp(k,:)); % filtering with phase compensation
        % extracting amplitude envelope out of x_hf
        amp_hf(k,:) = abs ( hilbert(x_hf) );
    end
end
if strcmp(pac_filter_type, 'fir1')
    for k =1 : size(lfp,1)
        x_hf=eegfilt(lfp(k,:),fs,f1_hz,f2_hz);
        % extracting amplitude envelope out of x_hf
        amp_hf(k,:) = abs ( hilbert(x_hf) );
    end
end




