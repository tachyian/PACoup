function [ phase_lf] = get_phase(lfp, f_lf, pac_filter_bw_lf, ...
    pac_filter_type, pac_filter_order,fs )
%GET_AMPLITUDE_AND_PHASE Summary of this function goes here
%   Detailed explanation goes here
phase_lf = zeros(size(lfp));
f1_hz = (f_lf-pac_filter_bw_lf/2) ;
f2_hz = (f_lf+pac_filter_bw_lf/2);

if f1_hz < 1
    f1_hz = 1;
    pac_filter_bw_lf2 = f_lf - 1;
    f2_hz = (f_lf+pac_filter_bw_lf2/2);
end
f1 = 2*f1_hz/fs;
f2 = 2*f2_hz/fs;

if strcmp(pac_filter_type, 'butter')
    [b_lf, a_lf]= butter(pac_filter_order, [f1 f2]);
    
    for k =1 : size(lfp,1)
        x_lf = filtfilt(b_lf,a_lf,lfp(k,:)); % filtering with phase compensation
        % extracting phase out of x_lf
        % phases in [-pi,pi] radians
        phase_lf(k,:) = angle ( hilbert(x_lf) );
    end
end
if strcmp(pac_filter_type, 'fir1')
    for k =1 : size(lfp,1)
        x_lf=eegfilt(lfp(k,:),fs,f1_hz,f2_hz);
        % extracting phase out of x_lf
        % phases in [-pi,pi] radians
        phase_lf(k,:) = angle ( hilbert(x_lf) );
    end
end




