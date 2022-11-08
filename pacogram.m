function [ miMat, phaseMat, time, lowFreq, phaseMax] = pacogram( lfp, twin, tstep, fs , config, plotType, range_mi)
%PACOGRAM Generates and/or plots the evolution of phase-amplitude coupling
%over time, in a plot referred to as pacogram. The horizontal axis is the time, the
%vertical axis is the low-frequency, and the color is the strength of PAC
%(modulation index)
%   
%   Inputs
%       lfp: signal to be analyzed
%
%       twin: time window to compute PAC (sec)
%
%       tstep: time increment to compute the pacogram
%
%       fs: sampling frequency
%
%       config: data structure with parameters to measure PAC. 
%
%       config: data structure with PAC computation parameters
%           config.hfreq0: Frequency for amplitude to be used (Hz)
%           config.hfreqf: Should be equal to config.hfreq0 (Hz)
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
%       plotType: Set to 0 if you do not want to plot the pacogram, to 1 if
%       you want to plot the pacogram, and to 2 if you want to plot the
%       pacogram and phase angle of PAC.
%
%   Outputs 
%       miMat: Matrix with modulation index values along time and
%       low-frequency samples
%
%       phaseMat: Matrix with phase angles along time and low-frequency
%       samples (rad)
%       
%       time: vector with values of time at which PAC was measured (sec)
%
%       lowFreq: vector vector with values of low-frequency at which PAC 
%       was measured (Hz)
%
%       phaseMax: vector with the maximum phase of PAC over time (rad)
%
%   Authors: David Escobar & Luke Johnson
%   Reviewed by Jaejin Lee
%   Added a line...

%%
fsize = 16;
miMat = [];
phaseMat = [];
phaseMax = [];
time = [];

dn = round(tstep*fs);
nwin = round(twin*fs);
count = 1;
time = [twin:dn/fs:length(lfp)/fs]-twin/2;
warningMsg =0;
for nn=nwin:dn:length(lfp)    
    
    config_local = config;
    config_local.hfreqf = config.hfreq0;
    
    dataPac = getPac( lfp(nn-nwin+1:nn), fs, config_local);
    mi_lfp = dataPac.mi;
    phase = dataPac.phase; 
    
    miMat(:, count) = mi_lfp;
    phaseMat(:, count) = phase;
    [~, ind_max] = max(mi_lfp); 
    phaseMax(count) =  phase(ind_max);
    count=count+1;
end
lowFreq = config.lfreq0:config.dlfreq:config.lfreqf;

if plotType >0 
     subplot(plotType, 1, 1); 
     [X1,Y1] = meshgrid(time,lowFreq); 
     n_interpx_pac = 20;
     n_interpy_pac = 20;
     stepx1 = diff(X1(1,1:2))/n_interpx_pac;
            stepy1 = diff(Y1(1:2,1))/n_interpy_pac;
            [plotx1, ploty1] = meshgrid(X1(1,1) : stepx1 : X1(1,end), Y1(1,1): stepy1 : Y1(end,1));                               
     plotz1 = interp2(X1,Y1,(miMat),plotx1,ploty1,'cubic');
            imagesc(plotx1(1,:),ploty1(:,1),plotz1); set(gca,'Ydir','normal')
     
     axis tight;
     contourcmap('jet');
     
     if sum(range_mi) == 0
         caxis('auto');
     else
        caxis(range_mi);
     end     
     %colormap('jet');
     ylabel('PAC freq. for phase (Hz)', 'fontsize',fsize);
     
     if plotType ==1
        xlabel('Time (sec)', 'fontsize',fsize);
     end
     set(gca, 'fontsize',fsize);
     
     hp = get(gca,'Position');
     cb = colorbar('Position', [hp(1)+hp(3)  hp(2)  0.015  hp(4)]);
     %cb = colorbar();
     ylabel(cb, 'Modulation Index (M.I.)', 'fontsize',fsize);
     
    if plotType ==2
        %axes( ha(plotType) ); 
        subplot(plotType, 1, 2); 
        plot(time, wrapTo2Pi(phaseMax)*180/pi, 'linewidth',3);
        xlabel('Time (sec)', 'fontsize',fsize);
        ylabel('Max M.I. phase (deg)', 'fontsize',fsize);
        set(gca, 'fontsize',fsize);
        axis tight;
    end 
end
end