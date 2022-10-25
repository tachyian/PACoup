function [ xPac ] = getSyntheticPac(t, width, phaseShift, fp, fa, Ap, Aa )
%GETSYNTHETICPAC Returns a one dimensional signal with phase-amplitude
% coupling
%
%   Inputs
%       t: time (sec)
%       width: proportional to the duration of the high-freqeuncy
%       oscillations at a particular phaseShift of the low-frequency oscillation
%       phase: phase angle of coupling (rad)
%       fp: frequency for phase (Hz)
%       fa: frequency for amplitude (Hz)
%       Ap: amplitude of low-frequency oscillations
%       Aa: amplitude of high-frequency oscillations dimensional vector
%       with phase amplitude coupling.
%
%   Outputs:
%       xPac: time series with phase-amplitude coupling
%
%   Based on: Adriano B.L. and others, "Measuring Phase-Amplitude Coupling
%   Between Neuronal Oscillations of Different Frequencies," 2010
%
%   Author: David Escobar / descobar@umn.edu

%%

if width > 0
    s = 2*(fp.*t - floor(0.5+fp.*t) ) ;
    sigma= width;
    mu=0;
    % Gaussean function
    phi = (1/(sigma*sqrt(2*pi)))*exp(-(s - mu).^2/(2*sigma^2) );
    g= (phi-min(phi))/(max(phi)-min(phi));
    %A_fa= Aa*( (1-xi)*g+xi) ; % initial version from Tort
    % I think xi is not needed as sigma (above) modulates the dsitribution of
    % the modulated signal
    A_fa= Aa*g ;
    x_hf = sin(2*pi*fa.*t);
    
    xPac = A_fa.*x_hf + Ap*sin(2*pi*fp.*t+ phaseShift);
    %x_pac = Aa*sin(2*pi*fa*t ) + Ap*sin(2*pi*fp*t+ phaseShift) ;
else
    xPac = Aa*(1+cos(2*pi*fp.*t)).*sin(2*pi*fa.*t) + Ap*sin(2*pi*fp.*t+phaseShift);
end
end

