function [dataPac] = getPac(lfpIn, fs, config)
%GETPAC Returns data structure with phase amplitude coupling (PAC) measurements
%
%   Inputs
%       lfpIn: local field potential time series
%
%       fs: sampling frequency
%
%       config: data structure with PAC computation parameters
%           config.hfreq0: Initial frequency for amplitude (Hz)
%           config.hfreqf: Final frequency for amplitude (Hz)
%           config.lfreq0: Initial frequency for phase (Hz)
%           config.lfreqf: Final frequency for amplitude (Hz)
%           config.dhfreq: Resolution of frequency for amplitude (Hz) 
%           config.dlfreq: Resolution of frequency for phase (Hz)
%           config.method: 'entropy', 'direct_pac', or 'mean_vector'
%           config.filterType: 'butter' or 'fir1'
%           config.filterLfBw: Bandwidth used to extract low-frequency
%           oscillations associated with the phase of PAC (Hz)
%           config.filterHfBw: Bandwidth used to extract high-frequency
%           oscillations associated with the amplitude of PAC (Hz). Set to 
%           zero if you want the HF bandwidth be equal to filterHfBwOffset+2*Low-freq
%           config.filterHfBwOffset: Offset for bandwidth when config.filterHfBw
%           is set to zero. 
%           config.filterOrder: order of filters used to extract low- and
%           high-frequency oscillations
%           config.entropyNumBins: Number of discretized values of phase
%           used in the 'entropy' method. 
%
%   Outputs
%       dataPac: data structure with information on PAC measurements
%           dataPac.mi: matrix with modulation index values over
%           frequencies for phase and amplitude 
%           dataPac.phase: matrix with phase angles of PAC over
%           frequencies for phase and amplitude (rad)
%           dataPac.miMax: maximum modulation index across frequencies
%           dataPac.flMax: frequency for phase associated with maximum
%           modulation index (dataPac.flMax) (Hz)
%           dataPac.fhMax: frequency for amplitude associated with maximum
%           modulation index (dataPac.flMax) (Hz)
%           dataPac.phaseMax: phase angle associated with maximum
%           modulation index (dataPac.flMax) (rad)
%           dataPac.lfreqGrid: grid with frequencies for phase used in PAC
%           measurements (Hz)
%           dataPac.hfreqGrid: grid with frequencies for amplitude used in PAC
%           measurements (Hz)
%
%   Author: David Escobar / descobar@umn.edu

%%
lfp = reshape(lfpIn, 1 , []);
dt= 1/fs;
%%
[X,Y] = meshgrid(config.lfreq0:config.dlfreq:config.lfreqf, config.hfreq0:config.dhfreq:config.hfreqf);
%%
miMatrix = zeros ( size(X,1), size(X,2)) ;
phaseMatrix= zeros(size(miMatrix,1),size(miMatrix,2));
%%
f_lf_ind = 1;

if config.filterHfBw == 0
    for f_lf = config.lfreq0: config.dlfreq :config.lfreqf
        phase_lf = get_phase(lfp(1,:), f_lf, config.filterLfBw, ...
            config.filterType, config.filterOrder,fs );
        f_hf_ind = 1;
        for f_hf = config.hfreq0: config.dhfreq :config.hfreqf
            
            amp_hf = get_amplitude(lfp(1,:), f_hf, config.filterHfBwOffset+2*f_lf, ...
                config.filterType, config.filterOrder,fs );
            
            [miMatrix(f_hf_ind,f_lf_ind, 1), phaseMatrix(f_hf_ind,f_lf_ind, 1)] = get_mi( phase_lf, amp_hf, ...
                config.method , config.entropyNumBins);
            
            f_hf_ind = f_hf_ind+1;
        end
        f_lf_ind = f_lf_ind+1;
    end
else
    [phase_lf,amp_hf] = get_amplitude_and_phase(lfp, config.lfreq0, ...
        config.dlfreq, config.lfreqf, config.hfreq0, config.dhfreq, config.hfreqf, ...
        config.filterLfBw, config.filterHfBw,config.filterType, ...
        config.filterOrder, fs ) ;
    for f_lf_ind=1:size(phase_lf,1)
        for f_hf_ind=1:size(amp_hf,1)
            [ miMatrix(f_hf_ind,f_lf_ind), phaseMatrix(f_hf_ind,f_lf_ind)] = ...
                get_mi( phase_lf(f_lf_ind,:), amp_hf(f_hf_ind,:), ...
                config.method , config.entropyNumBins);
        end
    end 
end
%%
mi_max= max(max(miMatrix));

[f_hf_max,f_lf_max] =  find(miMatrix==max( miMatrix(:) ) ) ;
if isempty (f_hf_max)
    f_hf_max = 1;
end

if isempty (f_lf_max)
    f_lf_max = 1;
end

f_hf_max = round(mean(f_hf_max));
f_lf_max = round(mean(f_lf_max));

phase_max = phaseMatrix(f_hf_max, f_lf_max);

f_hf_max = config.hfreq0 + (f_hf_max-1)*config.dhfreq;
f_lf_max = config.lfreq0 + (f_lf_max-1)*config.dlfreq;

%%
dataPac.miMax = mi_max;
dataPac.flMax = f_lf_max;
dataPac.fhMax = f_hf_max;
dataPac.phaseMax = phase_max;
dataPac.mi = miMatrix;
dataPac.phase = phaseMatrix;
dataPac.lfreqGrid = X;
dataPac.hfreqGrid = Y;
end