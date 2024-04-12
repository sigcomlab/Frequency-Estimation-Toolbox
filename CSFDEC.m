%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% CSFDEC function implements the complex single frequency-delay estimation and 
% cancellation (CSFDEC) algorithm [1].
%
% ** Syntax**
% T = CSFDEC(H, params)
%
% ** SIGNAL MODEL: H[m,n] = \sum_{l=0}^{L-1} A_l exp(1i*2*pi*m*fD_l)
% exp(-1i*2*pi*n*fr_l) + w[m,n]
%
% ** Input arguments**
% - H               : MxN measurement matrix
% - params          : a class containing
%   1. M            : first dimension of the matrix H (number of rows)
%   2. N            : second dimension of the matrix H (number of columns)
%   3. K            : number of 2D tones to search
%   4. LM           : oversampling factor along rows of H
%   5. LN           : oversampling factor along columns of H
%   6. M0           : params.M*params.LM;
%   7. N0           : params.N*params.LN;
%   8. mp           : params.mp = [m0 m0.^2 m0.^3]; where m0 = (0:1:par.M-1).';
%   9. np           : params.np = [n0; n0.^2; n0.^3]; where n0 = (0:1:par.N-1);
%   10. seq_m       : par.seq_m = 0:1:par.M0-1;
%   11. seq_n       : par.seq_n = 0:1:par.N0-1;
%   12. Vm0         : par.Vm0 = par.seq_m/par.M0;
%   13. Vn0         : par.Vn0 = par.seq_n/par.N0;
%   14. it_residual : number of iterations to compute frequency residuals
%   15. re_est      : number of re-estimations to be carried out
%     
%
% ** Output argument**
% - T: a class containing
%   1. fr: is a vector of length params.K containing the normalized
%   frequencies associated to the columns of the input matrix H
%   2. fD: is a vector of length params.K containing the normalized
%   frequencies associated to the rows of the input matrix H
%   3. A: is a vector of length params.K containing the complex amplitudes of
%   the params.K strongest spectral 2D tones.
%   4. time: scalar value related to the execution time of the algorithm
%   5. Ys: residual 2D spectrum at the end of each of the params.K+1
%   cancellation phases (dimension --> (params.M0)x(params.N0)x(params.K+1))
%
% ** Authors**  
% Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it) 
% Pasquale Di Viesti (pasquale.diviesti@unimore.it)
% Giorgio Guerzoni (giorgio.guerzoni@unimore.it)
% Michele Mirabella (michele.mirabella@unimore.it)
% Elia Vignoli (elia.vignoli@unimore.it)
%
% ** References **
% [1] M. Mirabella, P. D. Viesti, A. Davoli and G. M. Vitetta, 
% "An Approximate Maximum Likelihood Method for the Joint Estimation of 
% Range and Doppler of Multiple Targets in OFDM-Based Radar Systems," 
% in IEEE Transactions on Communications, vol. 71, no. 8, pp. 4862-4876, 
% Aug. 2023, doi: 10.1109/TCOMM.2023.3280562.
%
% v1.0- april 2024
% Copyright (c) Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Website: https://www.sigcom.unimore.it/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





function T = CSFDEC(H, params)

%Prepare estimation vectors
fr = zeros(params.K,1);
fd = zeros(params.K,1);
A = zeros(params.K,1);

% Useful definitions
M = params.M;
N = params.N;
M0 = params.M0;
N0 = params.N0;
im0 = 1:M;
im1 = (M +1): 2*M;
im2 = (2*M +1): 3*M;
im3 = (3*M +1): 4*M;
in0 = 1:N;
in1 = (N +1): 2*N;
in2 = (2*N +1): 3*N;
in3 = (3*N +1): 4*N;
normF = 1/(M*N);
n0 = params.np(1,:);
n1 = params.np(2,:);
n2 = params.np(3,:);
m0 = params.mp(:,1);
m1 = params.mp(:,2);
m2 = params.mp(:,3);
% START ALGORITHM
Tstart = tic();
Y = ifft(fft(H,M0,1),N0,2)*N0/(M*N);

x = [H H.*n0 H.*n1 H.*n2;
    H.*m0 H.*m0.*n0 H.*m0.*n1 H.*m0.*n2;
    H.*m1 H.*m1.*n0 H.*m1.*n1 H.*m1.*n2;
    H.*m2 H.*m2.*n0 H.*m2.*n1 H.*m2.*n2];

%initialize the variables

Yt = Y;
T.Ys = zeros(M0,N0,params.K+1);
T.Ys(:,:,1) = Y;
for qq = 1 : params.re_est
for k = 1 : params.K
    if qq==1
    % Initialization of A f and fd through periodogram method
    A(k) = max(Yt,[],[1 2]);
    [i_m, i_n] = find(Yt == A(k));
    fr(k) = params.Vn0(i_n);
    fd(k) = params.Vm0(i_m);
    end

    delta = 0;      % initial residual for range
    F_num = nnz(fr);
    notk = [1:1:k-1,k+1:1:F_num]; % get all other indexes
    Fr = fr(k);
    Fd = fd(k);
    for it = 1 : params.it_residual
        if isempty(notk)
            AF00_local = 0;
            AF01_local = 0;
            AF02_local = 0;
            AF11_local = 0;
            AF12_local = 0;
            AF13_local = 0;
            AF10_local = 0;
            AF20_local = 0;
            AF21_local = 0;
            AF22_local = 0;
            AF23_local = 0;
            AF31_local = 0;
            AF32_local = 0;
        else
            fr_diff = fr(notk)-Fr;
            fd_diff = fd(notk)-Fd;
            q.rY = exp(-1i*2*pi*fr_diff);
            q.dY = exp(1i*2*pi*fd_diff);
            to_lim.rY = abs(fr_diff)<eps;
            to_lim.dY = abs(fd_diff)<eps;
            [F0, F1, F2, F3] = residual2DFT (q, params, to_lim);

            AF00_local = getsum(A,notk,F0.d,F0.r);
            AF01_local = getsum(A,notk,F0.d,F1.r);
            AF02_local = getsum(A,notk,F0.d,F2.r);
            AF11_local = getsum(A,notk,F1.d,F1.r);
            AF12_local = getsum(A,notk,F1.d,F2.r);
            AF13_local = getsum(A,notk,F1.d,F3.r);
            AF10_local = getsum(A,notk,F1.d,F0.r);
            AF20_local = getsum(A,notk,F2.d,F0.r);
            AF21_local = getsum(A,notk,F2.d,F1.r);
            AF22_local = getsum(A,notk,F2.d,F2.r);
            AF23_local = getsum(A,notk,F2.d,F3.r);
            AF31_local = getsum(A,notk,F3.d,F1.r);
            AF32_local = getsum(A,notk,F3.d,F2.r);
        end


        C = normF *sum(sum(x(im0,in0).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF00_local;
        Y01a = normF *sum(sum(x(im0,in1).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF01_local;
        Y02a = normF *sum(sum(x(im0,in2).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF02_local;
        Y11a = normF *sum(sum(x(im1,in1).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF11_local;
        Y12a = normF *sum(sum(x(im1,in2).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF12_local;
        Y21a = normF *sum(sum(x(im2,in1).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF21_local;
        Y22a = normF *sum(sum(x(im2,in2).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF22_local;
        Y31a = normF *sum(sum(x(im3,in1).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF31_local;
        Y32a = normF *sum(sum(x(im3,in2).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF32_local;
        Y10a = normF *sum(sum(x(im1,in0).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF10_local;
        Y20a = normF *sum(sum(x(im2,in0).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF20_local;
        Y13a = normF *sum(sum(x(im1,in3).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF13_local;
        Y23a = normF *sum(sum(x(im2,in3).*exp(-1j*2*pi*params.mp(:,1)*Fd),1).*exp(1j*2*pi*params.np(1,:)*Fr),2) - AF23_local;

        d3 = (delta^3)/6;
        d2 = (delta^2)/2;

        bd = d3*imag(conj(C)*Y23a) - d2*real(conj(C)*Y22a) - delta*imag(conj(C)*Y21a) + real(conj(C)*Y20a);
        cd = d3*real(conj(C)*Y13a) + d2*imag(conj(C)*Y12a) - delta*real(conj(C)*Y11a) - imag(conj(C)*Y10a);

        omega = -cd/bd;

        o3 = (omega^3)/6;
        o2 = (omega^2)/2;

        br = -o3*imag(conj(C)*Y32a) - o2*real(conj(C)*Y22a) + omega*imag(conj(C)*Y12a) + real(conj(C)*Y02a);
        cr = o3*real(conj(C)*Y31a) - o2*imag(conj(C)*Y21a) - omega*real(conj(C)*Y11a) + imag(conj(C)*Y01a);

        delta = -cr/br;

        Fr = Fr + (delta/(2*pi)); % Range (coarse + fine)
        Fd = Fd + (omega/(2*pi)); % Doppler (coarse + fine)

    end % end for it = 1 : params.it_residual
    fr(k) = Fr;
    fd(k) = Fd;
    A(k) = C;

    % Cancellation

    idx2 = 1:k;

    q.rY = exp(-1i*2*pi*(fr(idx2)-params.Vn0));
    q.dY = exp(1i*2*pi*(fd(idx2)-params.Vm0));
    to_lim.rY = abs(fr(idx2)-params.Vn0)<eps;
    to_lim.dY = abs((fd(idx2)-params.Vm0))<eps;
    [F0, ~, ~, ~] = residual2DFT (q, params, to_lim);
    AF00_local = getsum(A,idx2,F0.d,F0.r);
    Yt = Y - AF00_local;
    if qq == params.re_est
        T.Ys(:,:,k+1) = Yt;
    end
end
end
T.fr = fr;
T.fD = fd;
T.A = A;
T.time = toc(Tstart);
end

%% Other functions

function AF = getsum(A,A_index,Fd,Fr)
jj = 1:numel(A_index);
AF = A(A_index).'.*Fd(jj,:).'*Fr(jj,:);
end

% ----------------------------------------------------------------

function [F0, F1, F2, F3] = residual2DFT(q, const, to_lim)

N = const.N;
M = const.M;
N2 = N^2;
M2 = M^2;
N3 = N^3;
M3 = M^3;

% prepare num powers for range and Doppler
qr = q.rY;             % qr
qrb = qr.^(N);         % qr^(N)
qrc = qr.^(N + 1);     % qr^(N+1)
qrd = qr.^(N + 2);     % qr^(N+2)
qre = qr.^(N + 3);     % qr^(N+3)

qd = q.dY;             % qr
qdb = qd.^(M);         % qr^(M)
qdc = qd.^(M + 1);     % qr^(M+1)
qdd = qd.^(M + 2);     % qr^(M+2)
qde = qd.^(M + 3);     % qr^(M+3)

% prepare den powers for range and Doppler
qr_1 = qr - 1;
den0r = qr_1;
den1r = qr_1.^2;
den2r = qr_1.^3;
den3r = qr_1.^4;

qd_1 = qd - 1;
den0d = qd_1;
den1d = qd_1.^2;
den2d = qd_1.^3;
den3d = qd_1.^4;

% Compute F0, F1, F2

% along range
F0.r = (qrb - 1) ./ den0r;
F1.r = ((N - 1) * qrc - N * qrb + qr) ./ den1r;
F2.r = (((N - 1)^2) * qrd + (-2*N2 + 2*N +1)*qrc + N2*qrb - qr.^2 - qr) ./ den2r;
F3.r = (((N - 1)^3) * qre + (-3*N3 + 6*N2 -4)*qrd +(3*N3-3*N2-3*N-1)*qrc -N3*qrb +qr.^3+4*qr.^2+qr) ./ den3r;

% along Doppler
F0.d = (qdb - 1) ./ den0d;
F1.d = ((M - 1) * qdc - M * qdb + qd) ./ den1d;
F2.d = (((M - 1)^2) * qdd + (-2*M2 + 2*M +1)*qdc + M2*qdb - qd.^2 - qd) ./ den2d;
F3.d = (((M - 1)^3) * qde + (-3*M3 + 6*M2 -4)*qdd +(3*M3-3*M2-3*M-1)*qdc -M3*qdb +qd.^3+4*qd.^2+qd) ./ den3d;

% check limits along range and Doppler
if any(any(to_lim.rY))
    F0.r(to_lim.rY) = const.N;
    F1.r(to_lim.rY) = (const.N^2 - const.N)/2;
    F2.r(to_lim.rY) = ((2*const.N-1)*(const.N-1)*const.N)/6;
    F3.r(to_lim.rY) = ((const.N-1)^2*const.N^2)/4;
end

if any(any(to_lim.dY))
    F0.d(to_lim.dY) = const.M;
    F1.d(to_lim.dY) = (const.M^2 - const.M)/2;
    F2.d(to_lim.dY) = ((2*const.M-1)*(const.M-1)*const.M)/6;
    F3.d(to_lim.dY) = ((const.M-1)^2*const.M^2)/4;
end


% final results These are Wk with k=0,1,2,3 (for range and Doppler)

F0.r = (F0.r)./ const.N;
F1.r = (F1.r)./ const.N;
F2.r = (F2.r)./ const.N;
F3.r = (F3.r)./ const.N;

F0.d = (F0.d)./ const.M;
F1.d = (F1.d)./ const.M;
F2.d = (F2.d)./ const.M;
F3.d = (F3.d)./ const.N;

end

% ----------------------------------------------------------------