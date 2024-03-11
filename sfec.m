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
end

%% function [F0] = residualDFT0(q, const, to_lim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F0] = residualDFT0(q, const, to_lim)

N = const.N;

qb = q.^(N);         

q_1 = q - 1;
den0 = q_1;

F0 = (qb - 1) ./ den0;

if any(any(to_lim))
    F0(to_lim) = const.N;
end

F0 = F0./ const.N;
end

%% function [F0, F1, F2] = residualDFT_new(q, const, to_lim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F0, F1, F2] = residualDFT_new(q, const, to_lim)

N = const.N;
N2 = N^2;
     
qb = q.^(N);           
qc = q.^(N + 1);         
qd = q.^(N + 2);         

q_1 = q - 1;
den0 = q_1;
den1 = q_1.^2;
den2 = q_1.^3;

F0 = (qb - 1) ./ den0;
F1 = ((N - 1) * qc - N * qb + q) ./ den1;
F2 = (((N - 1)^2) * qd + (-2*N2 + 2*N +1)*qc + N2*qb - q.^2 - q) ./ den2;

if any(any(to_lim))
    F0(to_lim) = const.N;
    F1(to_lim) = (const.N^2 - const.N)/2;
    F2(to_lim) = ((2*const.N-1)*(const.N-1)*const.N)/6;
end

F0 = F0./ const.N;
F1 = F1./ const.N;
F2 = F2./ const.N;
end

%% function [resh_input] = sfec_checkinput(input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% checkinput is a function to check the size, shape and values of the input
% parameters of sfec function.
%
% ** Syntax**
% [resh_input] = checkinput(input)
%
% ** Input arguments**
% - input: a class containing
%       1. input.x0: is a vector containing the (time-domain) signal to be
%       analysed.
%       2. input.Ts: is the sampling time (single value) of the signal 
%       params.X0.
%       3. input.M: is the oversampling factor (single value) that will be
%       employed for the FFTs.
%       4. input.K: is the number of peaks to find (single value).
%       5. input.NCFC: number of sfec refinement iterations 
%       (default=1).
%       6. input.f_min: is the lowerbound frequency (Hz) in the peaks 
%       search (none=0).
%       7. input.f_max: is the upperbound frequency (Hz) in the peaks 
%       search (none=0).          
%
% ** Output arguments**
% - resh_input: is a class containing the same attributes of input, but
% properly reshaped and checked.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resh_input] = sfec_checkinput(input)

    % check the input signal shape
    sz_x0 = size(input.x0);   
    if sz_x0(1)==1 % row vector
        resh_input.x0 = input.x0.';
    elseif sz_x0(2)==1 % column vector -> right shape
        resh_input.x0 = input.x0;
    else    % matrix and not vector
        disp('params.x0 should be a vector, got a matrix instead')
        disp('Press enter to go on')
        pause
    end
    % check the input signal type
    if ~isreal(resh_input.x0)
        disp('params.x0 should be real (double), got other type instead')
        disp('Press enter to go on')
        pause
    end

    
    % check the sampling time
    sz_Ts = size(input.Ts);    
    if sz_Ts(2)==1
        resh_input.Ts = input.Ts;
    else
        disp('params.Ts should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the oversampling factor
    sz_M = size(input.M);    
    if sz_M(2)==1
        resh_input.M = input.M;
    else
        disp('params.M should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the number of tones to search
    sz_K = size(input.K);    
    if sz_K(2)==1
        resh_input.K = input.K;
    else
        disp('params.K should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the number of refinement iterations
    sz_NCFC = size(input.NCFC);    
    if sz_NCFC(2)==1
        resh_input.NCFC = input.NCFC;
    else
        disp('params.NCFC should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the lowerbound for tones search
    sz_f_min = size(input.f_min);    
    if sz_f_min(2)==1
        resh_input.f_min = input.f_min;
    else
        disp('params.f_min should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the upperbound for tones search
    sz_f_max = size(input.f_max);    
    if sz_f_max(2)==1
        resh_input.f_max = input.f_max;
    else
        disp('params.f_max should be a value, got a vector or matrix instead')
        disp('Press enter to go on')
        pause
    end

    % check the values of upperbound and lowerbound
    if input.f_min > input.f_max
        disp('Error: params.f_max should be >= than params.f_min')
        disp('Press enter to go on')
        pause
    end
end

