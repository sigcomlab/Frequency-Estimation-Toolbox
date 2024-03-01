%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% sfec function implements the single frequency estimation and 
% cancellation (SFEC) algorithm (version for real input signal) [1,2,3].
%
% ** Syntax**
% [f, A, bin] = sfec(params)
%
% ** Input arguments**
% - params: a class containing
%       1. params.x0: is a vector containing the (time-domain) signal to be
%       analysed.
%       2. params.Ts: is the sampling time (single value) of the signal 
%       params.X0.
%       3. params.M: is the oversampling factor (single value) that will be
%       employed for the FFTs.
%       4. params_SFEC.K: is the number of peaks to find (single value).
%       5. params_SFEC.NCFC: number of sfec refinement iterations 
%       (default=1).
%       6. params.f_min: is the lowerbound frequency (Hz) in the peaks 
%       search (none=0).
%       7. params.f_max: is the upperbound frequency (Hz) in the peaks 
%       search (none=0).          
%
% ** Output arguments**
% - f: is a vector of length params.K containing the normalized frequencies
% of the params.K strongest spectral tones. To retrieve the frequency
% values in Hz, the formula f/params.Ts should be applied.
% - A: is a vector of length params.K containing the complex amplitudes of
% the params.K strongest spectral tones.
% - bin: is a vector of length params.K containing the FFT bin indices 
% associated to the frequencies f.
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
% [3] L. Ferrari, G. M. Vitetta, et al., (10/10/2022), Method for 
% two-dimensional and three-dimensional imaging based on collocated 
% multiple-input multiple-output radars, (World Intellectual Property
% Organization N.WO2022233888A1)
%
% v1.0 - march 2024
% Copyright (c) Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Website: https://www.sigcom.unimore.it/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, A, bin] = sfec(params)

% make sure the input are correctly shaped
params = sfec_checkinput(params);

% Initalize variables
params.N = size(params.x0,1);  
params.N0 = params.N * params.M;     % DFT order
params.T0 = params.Ts;
params.F_trial = [0:1:params.N0-1] / (params.N0-1);
params.F_trial_min = params.F_trial / params.Ts;
params.iwr_freq_vect = params.F_trial./params.Ts;
params.v0 = (0:1:params.N0-1)/(params.N0*params.Ts);   % FFT frequencies
params.V0 = (0:1:params.N0-1)/(params.N0);             % FFT normalized frequencies

sn1 = (0:1:params.N-1);
sn2 = (0:1:params.N-1) .* sn1;
sn3 = (0:1:params.N-1) .* sn2;

params.C1 = fft(sn1, params.N0)/params.N;
params.C2 = fft(sn2, params.N0)/params.N;
params.C3 = fft(sn3, params.N0)/params.N;

params.iwr_freq_vect = (0:params.N0-1)/(params.N0*params.Ts);

params.seq = (0:1:params.N0-1);
params.np = [sn1; sn2; sn3];

% find the min and max bin for the peaks search
[~, params.RMinIdx] = min(abs(params.v0(1:round(end/2)) - params.f_min));
if params.f_max == 0    % if no upperbound is set 
    params.f_max = params.v0(round(end/2));
end
[~, params.RMaxIdx] = min(abs(params.v0(1:round(end/2)) - params.f_max));
        
params.x0 = params.x0;
params.x1 = params.x0.*(0:1:params.N-1).';
params.x2 = params.x1.*(0:1:params.N-1).';
params.x3 = params.x2.*(0:1:params.N-1).';

params.s0 = fft(params.x0, params.N0, 1)./ params.N;
params.s1 = fft(params.x1, params.N0, 1)./ params.N;
params.s2 = fft(params.x2, params.N0, 1)./ params.N;
params.s3 = fft(params.x3, params.N0, 1)./ params.N;

sig = [params.s0.'; params.s1.'; params.s2.'];
x = [params.x0.'; params.x1.'; params.x2.'];

f = zeros(params.K, 1);
A = zeros(params.K, 1);

sigt = sig(1, :);

RMinIdx = params.RMinIdx;
RMaxIdx = params.RMaxIdx;

for qq = 1 : params.NCFC

    for k = 1 : params.K
        if qq == 1
            % Run CSFE2
            % Estimation of the initial value of alpha through periodogram method
            [A(k), alph_i] = max(sigt(RMinIdx:RMaxIdx));
            alph = alph_i + RMinIdx - 1;
            f(k) = params.V0(alph);
            bin(k) = alph;
        end
        
        if k > 1
            ord = 1:k;
            for jj = ord
                idx = [1:1:jj-1, jj+1:1:k];
                fdiff = f(idx)-f(jj);
                fdiff_c = f(idx)+f(jj);
                to_lim = abs(fdiff)<eps;
                to_lim_c = abs(fdiff_c)<eps;
                q = exp(1j*2*pi*fdiff);               
                q_c = exp(-1j * 2 * pi * fdiff_c);
                [F0, F1, F2] = residualDFT_new(q, params, to_lim);
                [F0_c, F1_c, F2_c] = residualDFT_new(q_c, params, to_lim_c);
                
                AF0_local = sum(A(idx) .* F0 + conj(A(idx)).* F0_c, 1);
                AF1_local = sum(A(idx) .* F1 + conj(A(idx)).* F1_c, 1);
                AF2_local = sum(A(idx) .* F2 + conj(A(idx)).* F2_c, 1);
                
                F = f(jj);
                C = 1/params.N*(sum(x(1,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                    - AF0_local;
                X1a = 1/params.N*(sum(x(2,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                    - AF1_local;
                X2a = 1/params.N*(sum(x(3,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                    - AF2_local;
                
                C1_2l = params.C1(2*bin(jj)-1);
                C2_2l = params.C2(2*bin(jj)-1);
                
                w = exp(-1i*4*pi*f(jj));
                gf = 1/params.N0*(w^params.N0-1)/(w-1);

                Ar = (real(C)*(1-real(gf))-imag(C)*(imag(gf)))/(1-abs(gf)^2);
                Ai = (imag(C)*(1+real(gf))-real(C)*(imag(gf)))/(1-abs(gf)^2);
                A(jj) = Ar+1i*Ai;
                
                for it = 1:5
                    
                    b = -real(conj(A(jj))*X2a)+2*real(conj(A(jj))^2*C2_2l);
                    c = imag(conj(A(jj))*X1a)-imag(conj(A(jj))^2*C1_2l);
                    Delta = -c/b;
                    F = F + (Delta/(2*pi));
                    fdiff = f(idx)-F;
                    fdiff_c = f(idx)+F;
                    to_lim = abs(fdiff)<eps;
                    to_lim_c = abs(fdiff_c)<eps;
                    q = exp(1j*2*pi*fdiff);
                    q_c = exp(-1j*2*pi*fdiff_c);
                    [F0, F1, F2] = residualDFT_new(q, params, to_lim);
                    [F0_c, F1_c, F2_c] = residualDFT_new(q_c, params, to_lim_c);
                
                    AF0_local = sum(A(idx) .* F0 + conj(A(idx)).* F0_c, 1);
                    AF1_local = sum(A(idx) .* F1 + conj(A(idx)).* F1_c, 1);
                    AF2_local = sum(A(idx) .* F2 + conj(A(idx)).* F2_c, 1);

                    C = 1/params.N*(sum(x(1,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                        - AF0_local;
                    X1a = 1/params.N*(sum(x(2,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                        - AF1_local;
                    X2a = 1/params.N*(sum(x(3,:).*exp(-1j*2*pi*params.np(1,:)*F))) ...
                        - AF2_local;
                    
                    w = exp(-1i*4*pi*F);
                    gf = 1/params.N0*(w^params.N0-1)/(w-1);

                    Ar = (real(C)*(1-real(gf))-imag(C)*(imag(gf)))/(1-abs(gf)^2);
                    Ai = (imag(C)*(1+real(gf))-real(C)*(imag(gf)))/(1-abs(gf)^2);
                    A(jj) = Ar+1i*Ai;

                end      
            end
        end
        
        idx2 = [1:k];
        q = exp(1j * 2 * pi * (f(idx2) - params.seq / params.N0));
        q_c = exp(1j * 2 * pi * (-f(idx2) - params.seq / params.N0));
        to_lim = abs(f(idx2)-params.seq/params.N0)<eps;
        to_lim_c = abs(-f(idx2)-params.seq/params.N0)<eps;
        [F0] = residualDFT0 (q, params, to_lim);
        [F0_c] = residualDFT0 (q_c, params, to_lim_c);
        AF0_local = sum(A(idx2) .* F0 + conj(A(idx2)) .* F0_c, 1);
        sigt = sig(1,:) - AF0_local; 

%% decomment to see the effect of SFEC on the input signal

%         figure()
%         subplot(311)
%         plot(abs(sig(1,:)))
%         ylabel('Original Spectum')
%         subplot(312)
%         plot(abs(AF0_local))
%         ylabel('Cancelling Tone')
%         subplot(313)
%         plot(abs(sigt))
%         ylabel('Residual Spectrum')
%         xlabel('Frequency (Hz)')
        
    end
end

