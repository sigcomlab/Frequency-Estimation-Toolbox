%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% test_csfdec is a main script to test the csfdec function. The script is
% organized in two parts, providing two different examples:
%   1) a signal composed by 4 spectral 2D tones (at different amplitudes) is
%   generated and analysed by the csfdec function [1].
%   2) a signal is generated in a Montecarlo simulation to verify the root
%   mean square error (RMSE) performance of the CSFDEC when 2 close 2D tones are
%   present. The RMSEs are compared against the Cramèr-Rao lower bound
%   (CRLB) for the considered scenario

%% EXAMPLE 1: 2D HARMONIC RETRIEVAL PROBLEM

clc, clear, close all
% signal generation --> SIGNAL MODEL: H[m,n] = \sum_{l=0}^{L-1} A_l exp(1i*2*pi*m*fD_l)
% exp(-1i*2*pi*n*fr_l) + w[m,n]
M = 256;
N = 128;
K = 4;

m = 0:M-1;
n = 0:N-1;
A = [1 0.85 0.75 0.6];
Fm_spacing = 3;                                                             % normalized Fm spacing
Fn_spacing = 2.5;                                                           % normalized Fn spacing
Fm = 2.3/M + Fm_spacing*(0:K-1)/M;                                          % vector Fm
Fn = 1.8/N + Fn_spacing*(0:K-1)/N;                                          % vector Fn
Hi = (A.*exp(1i*2*pi*m.'.*Fm))*exp(-1i*2*pi*n.*Fn.');                       % noiseless signal (MxN matrix)
H = awgn(Hi,20,"measured");                                                 % noisy signal (e.g., SNR = 20 dB)

% initialize CSFDEC parameters
LM = 4;                                                                     % oversampling along rows of H
LN = 4;                                                                     % oversampling along columns of H
it_residual = 10;                                                           % number of refinement steps for each 2D tone
re_est = 2;                                                                 % number of re-estimation steps
par = initialize_csfdec_par(M,N,K,LM,LN,it_residual,re_est);                % initialization of the parameters

% run the CSFDEC
T = CSFDEC(H, par);
disp(strcat("Estimation ended in ",num2str(T.time)," s"));

% PLOTS
[fn_est,idx_sort] = sort(T.fr,'ascend');                                    % sort frequencies Fn
fm_est = T.fD(idx_sort);                                                    % sort Fm vector according to Fn sorting
A_est = T.A(idx_sort);                                                      % sort A vector according to Fn sorting

% plot 1) show estimated 2D tones
figure()
hold on
scatter(Fm,Fn,25*A,'filled',"red",'DisplayName','Ground Truth')
scatter(fm_est,fn_est,25*abs(A_est),"cyan",'DisplayName','Estimates')
title('Estimation accuracy')
xlabel('F_M [.]')
ylabel('F_N [.]')
legend('Location','northwest')

% plot 2) residual spectral terms
for ii = 1:size(T.Ys,3)
    figure()
    mesh(abs(T.Ys(:,:,ii)))
    zlim([0 1.25*max(abs(A))])
    xlabel('m')
    ylabel('n')
    title(strcat("2D residual spectrum at step ",num2str(ii)))
end





%% SUPPORTING FUNCTIONS

function par = initialize_csfdec_par(M,N,K,LM,LN,it_residual,re_est)
% check input parameters before creating the par structure
par.M = M;
par.N = N;
par.K = K;
par.LM = LM;
par.LN = LN;
par.M0 = par.M*par.LM;
par.N0 = par.N*par.LN;
m0 = (0:1:par.M-1).';
par.mp = [m0 m0.^2 m0.^3];
n0 = (0:1:par.N-1);
par.np = [n0; n0.^2; n0.^3];
par.seq_m = 0:1:par.M0-1;
par.seq_n = 0:1:par.N0-1;
par.Vm0 = par.seq_m/par.M0;
par.Vn0 = par.seq_n/par.N0;
par.it_residual = it_residual;
par.re_est = re_est;
end




