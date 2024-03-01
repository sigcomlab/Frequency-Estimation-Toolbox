%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% test_csfec is a main script to test the csfec function. The script is
% organized in two parts, providing two different examples:
%   1) a signal composed by 5 spectral tones (at different amplitudes) is
%   generated and analysed by the csfec function [1,2].
%   2) the file radar_unwrapped_phase.mat, containing real measurements of 
%   range profile from a single virtual antenna[1] when a single target is 
%   present within the field-of-view, is loaded. The csfec function is 
%   employed to precisely locate the target. The priori information of 
%   the target position is employed to limit the csfec search.
%   
% ** Authors**  
% Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it) 
% Pasquale Di Viesti (pasquale.diviesti@unimore.it)
% Giorgio Guerzoni (giorgio.guerzoni@unimore.it)
% Michele Mirabella (michele.mirabella@unimore.it)
% Elia Vignoli (elia.vignoli@unimore.it)
%
% ** References **
% [1] P. Di Viesti, A. Davoli, G. Guerzoni and G. M. Vitetta, "Novel 
% Deterministic Detection and Estimation Algorithms for Colocated 
% Multiple-Input Multiple-Output Radars," in IEEE Access, vol. 10, 
% pp. 2216-2255, 2022, doi: 10.1109/ACCESS.2021.3139200.
% [2] P. Di Viesti, A. Davoli, G. Guerzoni, et al. "Novel Methods for 
% Approximate Maximum Likelihood Estimation of Multiple Superimposed 
% Undamped Tones and Their Application to Radar Systems," TechRxiv, August 
% 03, 2021, doi: 10.36227/techrxiv.15054321.v2
% [3] G. Paterniani et al., "Radar-Based Monitoring of Vital Signs: A 
% Tutorial Overview," in Proceedings of the IEEE, vol. 111, no. 3, 
% pp. 277-317, March 2023, doi: 10.1109/JPROC.2023.3244362.
%
% v1.0 - march 2024
% Copyright (c) Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Website: https://www.sigcom.unimore.it/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% EXAMPLE 1
% generate a synthetic signal
params1.Ts = 0.001;
tstop = 20;
t = 0:params1.Ts:20;
shape = length(t);
f_comp = [3, 7.5, 8.99, 15, 42];
a_comp = [1.2, 1.9, 3, 1.15, 2.3];
params1.x0 = sum(a_comp.*exp(1i*2*pi*f_comp.*t.'),2);

% SFEC parameters
params1.M = 16;    % SFEC oversampling factor
params1.K = 5;     % four peaks to find
params1.NCFC = 2;  % SFEC refinement cycles
params1.f_min = 0;     % no lowerbound
params1.f_max = 0;   % no upperbound

[f, A, bin] = csfec(params1);

freq_hz = f / params1.Ts;  % equivalent to bin/(params.M*shape*params.Ts)
ampli = abs(A);

% plot the results
spectrum = abs(fft(params1.x0))/shape;
freq_axis = [0:shape-1]/(shape*params1.Ts);

figure
subplot(311)
plot(t, abs(params1.x0))     % time domain signal (absolute value)
xlabel('Time (s)')
ylabel('Absolute value')
title('EXAMPLE 1')
subplot(312)
plot(t, angle(params1.x0))     % time domain signal (phase)
xlabel('Time (s)')
ylabel('Phase (rad)')
subplot(313), hold on
plot(freq_axis, spectrum)   % frequency domain signal
scatter(freq_hz, abs(A))    % peaks found by SFEC
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% EXAMPLE 2
% load the range profile
frame = load ('radar_frame.mat');
params2.x0 = frame.radar_frame;
shape = length(params2.x0); 
d = 1.26;    % (m) target distance ground truth
min_r = 0.5;    % (m) range lowerbound (a priori information)
max_r = 2;    % (m) range upperbound (a priori information)
c = 3e8;

% radar parameter
params2.Ts = 1/1000000;    % (s) radar sampling time 
params2.kf = 7.8125e11;    % (Hz/s) radar chirp slope

% SFEC parameters
params2.M = 16;    % SFEC oversampling factor
params2.K = 1;     % a single peak to find (target)
params2.NCFC = 1;  % SFEC refinement cycles

frq_hz = [0:shape-1]./((shape-1)*params2.Ts); 
frq_r = frq_hz*c/(2*params2.kf);  % range bins

f_min = min_r*(2*params2.kf)/c;
f_max = max_r*(2*params2.kf)/c;

params2.f_min = f_min;     % (Hz) minimum range to search
params2.f_max = f_max;   % (Hz) maximum range to search

[f_br, A_br, bin_br] = csfec(params2);

d_est = f_br * c / (2*params2.Ts*params2.kf);  % equivalent to ???

% plot the results
N0 = shape * params2.M;
time_axis = [0:shape-1]*params2.Ts;     % (s) time axis
spectrum = abs(fft(params2.x0, N0))/shape;   % spectrum of params.x0
frq_hz = [0:shape*params2.M-1]./((shape*params2.M-1)*params2.Ts); 
frq_r = frq_hz*c/(2*params2.kf);  % range bin

figure()
subplot(211)
plot(time_axis, abs(params2.x0))     % time domain signal
xlabel('Time (s)')
ylabel('Displacement (rad)')
title('EXAMPLE 2')
subplot(212), hold on
plot(frq_r, spectrum)       % frequency domain signal
scatter(d_est, abs(A_br), 'red')    % target distance found by SFEC
xline(d, 'k')        % target distance ground truth as vertical lines
xlabel('Range (m)')
ylabel('Amplitude')
