%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% DRAEC3D function implements the reduced dimension version of the Doppler range angle estimator
% with successive compensation algorithm [1].
%
% ** Syntax**
% T = DRAEC3D(H, params)
%
% ** SIGNAL MODEL: H[m,n,nt] = \sum_{l=0}^{L-1} A_l exp(1i*2*pi*m*fD_l)
% exp(-1i*2*pi*n*fr_l) exp(-1i*2*pi*nt*f_theta_l) + w[m,n,nt]
%
% ** Input arguments**
% - H : M-by-N-by-Nant 3D measurement tensor/matrix whose element (m,n,nt)
%       obeys to the model described above
% - params: set of parameters initialized by the function
%       'initialize_draec3D_par'
%
% ** Output argument**
% - T: a class containing
%   1. fr: is a vector of length params.K containing the normalized
%   frequencies associated to the 1st dimension of the input matrix H
%   2. fD: is a vector of length params.K containing the normalized
%   frequencies associated to the 2nd dimension of the input matrix H
%   3. f_theta: is a vector of length params.K containing the normalized
%   frequencies associated to the 3rd dimension of the input matrix H
%   3. A: is a vector of length params.K containing the complex amplitudes of
%   the params.K strongest spectral 3D tones.
%   4. time: scalar value related to the execution time of the algorithm
%
% ** Authors**
% Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Pasquale Di Viesti (pasquale.diviesti@unimore.it)
% Giorgio Guerzoni (giorgio.guerzoni@unimore.it)
% Michele Mirabella (michele.mirabella@unimore.it)
% Elia Vignoli (elia.vignoli@unimore.it)
%
% ** References **
% [1] M. Mirabella, P. Di Viesti and G. M. Vitetta, "Deterministic Algorithms
% for Four-Dimensional Imaging in Colocated MIMO OFDM-Based Radar Systems," in
% IEEE Open Journal of the Communications Society, vol. 4, pp. 1516-1543,
% 2023, doi: 10.1109/OJCOMS.2023.3292796.
%
% v1.0 - april 2024
% Copyright (c) Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Website: https://www.sigcom.unimore.it/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function T = DRAEC3D(H,params)

T.f_theta = NaN(params.csfdec_settings.K,params.csfec_settings.K);
T.fD = NaN(params.csfdec_settings.K,params.csfec_settings.K);
T.fr = NaN(params.csfdec_settings.K,params.csfec_settings.K);
T.A = NaN(params.csfdec_settings.K,params.csfec_settings.K);
T.cl_freq = 0;
Tstart = tic();
Hr = H(:,:,params.ref_ant);                                         % fix a reference sample in angular domain
T1 = CSFDEC(Hr, params.csfdec_settings);                            % run the CSFDEC algorithm
AA = abs(T1.A).^2;
TPi = AA/max(AA) > params.det_th2D;
T1.fD = T1.fD(TPi,1);
T1.fr = T1.fr(TPi,1);

for ii=1:sum(TPi)
    Af2 = (exp(1i*2*pi*T1.fD(ii,1)*params.csfdec_settings.mp(:,1).')).'*exp(-1i*2*pi*T1.fr(ii,1)*params.csfdec_settings.np(1,:)); % compensating factor (MxN)
    Hf2 = squeeze(mean(H.*conj(Af2),[1 2])).';
    [ft, A] = CSFECn(Hf2, params.csfec_settings);
    AA = abs(A).^2;
    APi = AA/max(AA) > params.det_th1D;
    T.A(ii,APi) = A(APi);
    T.f_theta(ii,APi) = ft(APi);
end

T.f_theta = T.f_theta(~isnan(T.f_theta));
NT = numel(T.f_theta);

Af1 = permute(exp(-1i*2*pi*T.f_theta.*params.csfec_settings.np(1,:)),[3 4 2 1]);  % 1-by-1-Nant-by-K compensating factor
for jj = 1:NT
    Hf1 = squeeze(mean(H.*conj(Af1(:,:,:,jj)),3));                                  % sum along the antenna elements
    T2 = CSFDEC(Hf1, params.csfdec_settings);
    AA = abs(T2.A).^2;
    TPi2 = AA/max(AA) > params.det_th2D;
    N_found = sum(TPi2);
    T.fD(1:N_found,jj) = T2.fD(TPi2,1);
    T.fr(1:N_found,jj) = T2.fr(TPi2,1);
    T.A(1:N_found,jj) = T2.A(TPi2,1);
end
T.time = toc(Tstart);
% sort and return
index_duplicates = ~isnan(T.fr);
if NT == sum(index_duplicates,"all")
    T.fr = T.fr(index_duplicates);
    T.fD = T.fD(index_duplicates);
    T.A = T.A(index_duplicates);
    [T.fr,idx_sort] = sort(T.fr,'ascend');
    T.fD = T.fD(idx_sort);
    T.A = T.A(idx_sort);
    T.f_theta = T.f_theta(idx_sort);
else
    disp('Close frequencies detected');
    T.cl_freq = 1;
end
end


% SUPPORTING FUNCTIONS

function [f, A] = CSFECn(H, const)

% CSFECp - Computation of a set of delays and complex amplitudes given noisy samples
%
% Syntax:  [f, A] = CSFECp(sig, x, it_residual, re_est, const)
%
% Inputs:
%    sig - vector of samples of the IDFT of the baseband signal
%    x - 3xN matrix: [x; n*x; n^2*x];
%    const - set of constants and parameters
%
% Outputs:
%    f - vector of delays
%    A - vector of corresponding complex amplitudes
%
% Example:
%    [f, A] = CSFEC(sig, x, 20, 5, const)


% Initalize variables
f = zeros(const.K,1);
A = zeros(const.K,1);
sig = ifft(H,const.N0,2).*const.LN;
x = [H;const.np(1,:).*H;const.np(2,:).*H];

sigt = sig;

for qq = 1 : const.re_est
    % Coarse estimation of f
    for k = 1 : const.K
        if qq == 1
            % Run CSFE
            % Estimation of the initial value of alpha through periodogram method
            [A(k), alph] = max(sigt);
            f(k) = const.V0(alph);
        end
        F_num = nnz(f);
        notk = [1:1:k-1,k+1:1:F_num]; % get all other indexes
        F = f(k);

        for it = 1:const.it_residual
            if isempty(notk)
                AF0_local = 0;
                AF1_local = 0;
                AF2_local = 0;
            else
                fdiff = f(notk) - F;
                q = exp(-1j*2*pi*fdiff);
                to_lim = abs(fdiff)<eps;
                [F0, F1, F2] = residualDFT (q, const, to_lim);
                AF0_local = sum(A(notk) .* F0, 1);
                AF1_local = sum(A(notk) .* F1, 1);
                AF2_local = sum(A(notk) .* F2, 1);
            end
            C = 1/const.N*(sum(x(1,:).*exp(1j*2*pi*const.np(1,:)*F))) - AF0_local;
            X1a = 1/const.N*(sum(x(2,:).*exp(1j*2*pi*const.np(1,:)*F))) - AF1_local;
            X2a = 1/const.N*(sum(x(3,:).*exp(1j*2*pi*const.np(1,:)*F))) - AF2_local;
            % Compute leakage terms

            % Computation of a, b
            b = -real(conj(C)*X2a);
            c = -imag(conj(C)*X1a);
            Delta = -c/b;
            F = F + (Delta/(2*pi)); % frequency estimate (coarse + fine)
        end

        f(k) = F;
        A(k) = C;
        % cancellation
        idx2 = 1:k;
        q = exp(-1j*2*pi*(f(idx2)-const.V0));
        to_lim = abs(f(idx2)-const.V0)<eps;
        [F0, ~, ~] = residualDFT (q, const, to_lim);
        AF0_local = sum(A(idx2) .* F0, 1);
        sigt = sig - AF0_local;

    end
end
end

% -------------------------------------------------------------------------

function [F0, F1, F2] = residualDFT(q, const, to_lim)
N = const.N;
N2 = N^2;

% prepare num powers for range and Doppler
qrb = q.^(N);         % qr^(N)
qrc = q.^(N + 1);     % qr^(N+1)
qrd = q.^(N + 2);     % qr^(N+2)

% prepare den powers for range
qr_1 = q - 1;
den0r = qr_1;
den1r = qr_1.^2;
den2r = qr_1.^3;

% Compute F0, F1, F2

% along range
F0 = (qrb - 1) ./ den0r;
F1 = ((N - 1) * qrc - N * qrb + q) ./ den1r;
F2 = (((N - 1)^2) * qrd + (-2*N2 + 2*N +1)*qrc + N2*qrb - q.^2 - q) ./ den2r;

% check limits along range
if any(any(to_lim))
    F0(to_lim) = N;
    F1(to_lim) = (N2 - N)/2;
    F2(to_lim) = ((2*N-1)*(N-1)*N)/6;
end

F0 = (F0)./ N;
F1 = (F1)./ N;
F2 = (F2)./ N;
end

% -------------------------------------------------------------------------