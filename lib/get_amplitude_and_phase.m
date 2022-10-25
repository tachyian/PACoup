function [ phase_lf ,amp_hf ] = get_amplitude_and_phase(lfp, ...
    pac_f_lf0, pac_df_lf, pac_f_lff, pac_f_hf0, pac_df_hf, ...
    pac_f_hff, pac_filter_bw_lf, pac_filter_bw_hf, ...
    pac_filter_type, pac_filter_order,fs )
%GET_AMPLITUDE_AND_PHASE Summary of this function goes here
%   Detailed explanation goes here

phase_lf = [];
amp_hf = [];

% indexes to store modulation indexes in a matrix
f_lf_ind=1; %
f_hf_ind=1; %

%lfp_tmp(:) = lfp(:,k) ;
for f_lf=pac_f_lf0:pac_df_lf:pac_f_lff
    
    % low frequency
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
            phase_lf(f_lf_ind,:,k) = angle ( hilbert(x_lf) );
        end
    end
    if strcmp(pac_filter_type, 'fir1')
        for k =1 : size(lfp,1)
            x_lf=eegfilt(lfp(k,:),fs,f1_hz,f2_hz);
            % extracting phase out of x_lf
            
            % phases in [-pi,pi] radians
            phase_lf(f_lf_ind,:,k) = angle ( hilbert(x_lf) );
        end
    end

    
    f_lf_ind = f_lf_ind+1;
end

%%
% filtering signals across high frequency range
for f_hf=pac_f_hf0:pac_df_hf:pac_f_hff
    
    f1_hz = (f_hf-pac_filter_bw_hf/2) ;
    f2_hz = (f_hf+pac_filter_bw_hf/2) ;
    
    if f1_hz < 1
        f1_hz = 1;
        pac_filter_bw_hf2 = f_hf - 1;
        f2_hz = (f_hf+pac_filter_bw_hf2/2);
    end
    
    if f2_hz > fs/2;
        f2_hz = (fs/2)*0.95;
    end
    
    f1 = 2*f1_hz/fs;
    f2 = 2*f2_hz/fs;
    
    if strcmp(pac_filter_type, 'butter')
        [b_hf, a_hf] = butter(pac_filter_order,[f1 f2]);
        for k =1 : size(lfp,1)
            x_hf = filtfilt(b_hf,a_hf,lfp(k,:)); % filtering with phase compensation
            % extracting amplitude envelope out of x_hf
            amp_hf(f_hf_ind,:,k) = abs ( hilbert(x_hf) );
        end
    end
    if strcmp(pac_filter_type, 'fir1')
        for k =1 : size(lfp,1)
            x_hf=eegfilt(lfp(k,:),fs,f1_hz,f2_hz);
            % extracting amplitude envelope out of x_hf
            amp_hf(f_hf_ind,:,k) = abs ( hilbert(x_hf) );
        end
    end
    f_hf_ind = f_hf_ind+1;
end



