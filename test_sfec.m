%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% test_sfec is a main script to test the sfec function. The script is
% organized in two parts, providing two different examples:
%   1) a signal composed by 5 spectral tones (at different amplitudes) is
%   generated and analysed by the sfec function [1,2].
%   2) the file radar_unwrapped_phase.mat, containing real measurements of 
%   the unwrapped phase signal (displacement) from contactless vital signs 
%   monitoring [1], is loaded. The sfec function is employed to precisely 
%   locate the breath rate (BR) and heart rate (HR).
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
params1.x0 = sum(a_comp.*sin(2*pi*f_comp.*t.'),2);

% SFEC parameters
params1.M = 16;    % SFEC oversampling factor
params1.K = 5;     % four peaks to find
params1.NCFC = 2;  % SFEC refinement cycles
params1.f_min = 0;     % no lowerbound
params1.f_max = 0;   % no upperbound

[f, A, bin] = sfec(params1);

freq_hz = f / params1.Ts;  % equivalent to bin/(params.M*shape*params.Ts)
ampli = abs(A);

% plot the results
spectrum = abs(fft(params1.x0))/shape;
freq_axis = [0:shape-1]/(shape*params1.Ts);

figure
subplot(211)
plot(t, params1.x0)     % time domain signal
xlabel('Time (s)')
ylabel('Amplitude')
title('EXAMPLE 1')
subplot(212), hold on
plot(freq_axis, spectrum)   % frequency domain signal
scatter(freq_hz, abs(A))    % peaks found by SFEC
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%% EXAMPLE 2
% load the unwrapped phase signal
displacement = load ('radar_unwrapped_phase.mat');
params2.x0 = displacement.radar_unwrapped_phase;
shape = length(params2.x0); 
br = 12.36;    % (acts/min) BR ground truth
hr = 83.29;    % (BPM) HR ground truth

% radar parameter
params2.Ts = 0.0632;    % (s) radar sampling time 

% A. first SFEC run to find BR (limited to breath frequency region)
% SFEC parameters
params2.M = 16;    % SFEC oversampling factor
params2.K = 1;     % a single peak to find (BR)
params2.NCFC = 1;  % SFEC refinement cycles
params2.f_min = 0.1;     % (Hz) minimum physiological frequency of BR
params2.f_max = 0.5;   % (Hz) maximum physiological frequency of BR

[f_br, A_br, bin_br] = sfec(params2);


freq_br_bpm = f_br / params2.Ts * 60;  % equivalent to bin_br/(params.M*shape*params.Ts)*60

% plot the results
N0 = shape * params2.M;
time_axis = [0:shape-1]*params2.Ts;     % (s) time axis
spectrum = abs(fft(params2.x0, N0))/shape;   % spectrum of params.x0
freq_axis = [0:N0-1]/(N0*params2.Ts)*60;     % frequency axis in beats (or acts) per minute (BPM)

figure()
subplot(211)
plot(time_axis, params2.x0)     % time domain signal
xlabel('Time (s)')
ylabel('Displacement (rad)')
title('EXAMPLE 2a - Breath rate estimation')
subplot(212), hold on
plot(freq_axis, spectrum)       % frequency domain signal
scatter(freq_br_bpm, abs(A_br), 'red')    % BR found by SFEC
xline(br, 'k')        % BR ground truth as vertical lines
xlabel('Frequency (BPM)')
ylabel('Amplitude')

% B. second SFEC run to find HR (limited to heart frequency region)
% SFEC parameters
params2.M = 16;    % SFEC oversampling factor
params2.K = 1;     % a single peak to find (HR)
params2.NCFC = 1;  % SFEC refinement cycles
params2.f_min = 0.833;     % (Hz) minimum physiological frequency of HR
params2.f_max = 2.667;   % (Hz) maximum physiological frequency of HR

[f_hr, A_hr, bin_hr] = sfec(params2);

freq_hr_bpm = f_hr / params2.Ts * 60;  % equivalent to bin_hr/(params.M*shape*params.Ts)*60

% plot the results
figure()
subplot(211)
plot(time_axis, params2.x0)     % time domain signal
xlabel('Time (s)')
ylabel('Displacement (rad)')
title('EXAMPLE 2b - Heart rate estimation')
subplot(212), hold on
plot(freq_axis, spectrum)       % frequency domain signal
scatter(freq_hr_bpm, abs(A_hr), 'red')    % HR found by SFEC
xline(hr, 'k')        % HR ground truth as vertical lines
xlabel('Frequency (BPM)')
ylabel('Amplitude')
