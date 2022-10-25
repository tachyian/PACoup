%%   
%   This example script presents an overview of the functions used to
%   aritificially generate phase-amplitude coupling (PAC) and measure it.
%
%   Author:     David Escobar
%               Department of Neurology
%               University of Minnesota
%%
close all;
clear var;
clc;
%%
addpath('../synthetic/');
addpath('lib/');
%%
% Defining parameters of the synthetic PAC signal

%sampling time
dt=0.0003;  % (sec)
fs= 1/dt;

%final time
tfinal = 30; % (sec)

% width: a value close to zero shrinks the region in which the high-frequency
% component has high amplitude oscillations. This width is the standard deviation of a
% Gaussean function. Use a width equal to zero to generate a PAC signal
% using trigonometric functions.  
width = 0.2; % (rad)

% phase: phase of the phase carrier signal at which we want coupling
phase = 0*pi/4; % (rad)

% f_p/Ap: frequency of the phase carrier signal (low frequency)
% The frequency increases linearly from fp0 to fpf
fp0 = 20; % (Hz)
fpf = 20.1; % (Hz)
Ap = 10;

% f_a/Aa: frequency/amplitude of the high frequency modulated component
% The frequency increases linearly from f_a0 to f_af
f_a0 = 150; % (Hz)
f_af = 150.1; % (Hz)
Aa =5;

% noise level
noiseLevel = 2;

%%
%   Calling the function that generates synthetic PAC signaals

t= 0:dt:tfinal;
N = length(t); 

fp0 = double(fp0);
fpf = double(fpf);

beta = (fpf-fp0)/tfinal; 
f_p = fp0 + beta*t/2; % time varying frequency

f_a0 = double(f_a0);
f_af = double(f_af);

beta = (f_af-f_a0)/tfinal; 
f_a = f_a0 + beta*t/2; % time varying frequency
[ xPac] = getSyntheticPac( t, 0*width, phase, f_p, f_a, Ap, Aa );
%%
% adding noise
fnoise = get1overfNoise( length(t), fs ); 
xPac = xPac + noiseLevel*fnoise;
%%
% plotting signal
figure;
plot(t,xPac, 'r', 'linewidth',2); hold on;
xlabel('Time (sec)','fontsize',16); 
ylabel('Amplitude','fontsize',16); 
legend('PAC synthetic signal');
%% 
% plotting frequencies for PAC
figure;
plot(t, f_p, 'r--', 'linewidth',3); hold on;
plot(t, f_a, 'b-', 'linewidth',3); hold on;
xlabel('Time (sec)','fontsize',16); 
legend('Phase freq', 'Amp. freq.'); 

%%
% PAC computation parameters:
configPac.hfreq0 = 20; % initial frequency for amplitude
configPac.hfreqf = 300; % final frequency for amplitude
configPac.lfreq0 = 4; % initial frequency for phase
configPac.lfreqf = 30; % final frequency for phase
configPac.dhfreq = 5; % increment in freq. for amplitude
configPac.dlfreq = 1; % increment in freq. for phase
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
twin = 15; % 
tstep = 2; % start with a bigger step 
plotType = 2;
rangeMi = [];
figure;
configPac.filterHfBw=0;
configPac.hfreq0 = f_a0; % freq. for amplitude - found in comodulogram
%configPac.hfreqf = f_a0; 
[miMat, phaseMat, time, lowFreq, phaseMax] = pacogram(xPac, twin, tstep, fs, configPac, plotType, rangeMi);