%%
%   This example script presents an overview of the functions used to
%   aritificially generate phase-amplitude coupling (PAC) and measure it.
%
%   Author:     David Escobar
%               Department of Neurology
%               University of Minnesota
%%
close all;
clear all;
clc;
%%
addpath('../synthetic/');
addpath('lib/');

load ('Z:\descobar\data\pd\resting\offmed_offstim\Jo\Jo_20160104_2/Jo_20160104_2_all_sync');

m1array_to_remove_jo =[9,1; 4,2 ; 8,4; 4,8; 8,10];

fs = data.lfp_stn_fs(1);
dt=1/fs;
n0 = round(3*fs);

% High pass filter
fcut = 2*1/fs;
[bhp,ahp] = butter(2 , fcut , 'high');

%fs
%fcut = 2*700/fs;
%[blphf,alphf] = butter(2 , fcut );
%%
lfp_mean = get_lfp_mean( data.lfp_array, m1array_to_remove_jo);
lfp_m1 = lfp_mean(n0:end);                
xPac = filtfilt(bhp,ahp, lfp_m1);
            
figure;
plot(lfp_m1);

%
%%
% PAC computation parameters:
configPac.hfreq0 = 20; % initial frequency for amplitude
configPac.hfreqf = 300; % final frequency for amplitude
configPac.lfreq0 = 8; % initial frequency for phase8
configPac.lfreqf = 30; % final frequency for phase
configPac.dhfreq = 5; % increment in freq. for amplitude
configPac.dlfreq = 2; % increment in freq. for phase
configPac.method =  'entropy';
configPac.filterType= 'butter';
configPac.filterOrder= 2;
configPac.entropyNumBins=36;
configPac.filterLfBw=1;

% Set to zero if you want the HF bandwidth be selected automatically
% (filterHfBw = filterHfBwOffset+2*Low-freq.)
configPac.filterHfBw=0;

% This parameter is the offset used if the HF bandwdth is set automatically
configPac.filterHfBwOffset=10;

%%
dataPac = getPac( xPac, fs, configPac );
%%
contourPlot = 1;
titleStr = 'PAC synthetic';
rangeMi =0;
labelsOn = 1;
tickOn=1;
colorbarOn = 1;
figure;
comodulogram(dataPac.lfreqGrid, dataPac.hfreqGrid, dataPac.mi, contourPlot, rangeMi , labelsOn, tickOn,  colorbarOn );
%%
% Generating pacogram
twin = 10; %
tstep = 4; % start with a bigger step
plotType = 2;
rangeMi = [0 1.5e-3];
figure;
configPac.filterHfBw=0;
configPac.hfreq0 = 50; % freq. for amplitude - found in comodulogram
%configPac.hfreqf = f_a0;
[miMat, phaseMat, time, lowFreq, phaseMax] = pacogram(xPac, twin, tstep, fs, configPac, plotType, rangeMi);

%}