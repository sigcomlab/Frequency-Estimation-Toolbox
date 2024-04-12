%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ** Description **
% test_draec3D is a main script to test the DRAEC3D function. The script is
% organized in two parts, providing two different examples:
%   1) a signal composed by 4 spectral 3D tones (at different amplitudes) is
%   generated and analysed by the DRAEC3D function [1],[2],[3].
%
% NOTES : the DRAEC3D is capable of detecting similar frequencies within a
% specific dimension of the signal. For the sake of simplicity, the
% provided example consider distinc frequencies only.
%
% REFERENCES :
%
% [1] M. Mirabella, P. Di Viesti and G. M. Vitetta, "Deterministic Algorithms
% for Four-Dimensional Imaging in Colocated MIMO OFDM-Based Radar Systems," in
% IEEE Open Journal of the Communications Society, vol. 4, pp. 1516-1543,
% 2023, doi: 10.1109/OJCOMS.2023.3292796.
%
% [2] M. Mirabella, P. D. Viesti, A. Davoli and G. M. Vitetta, 
% "An Approximate Maximum Likelihood Method for the Joint Estimation of 
% Range and Doppler of Multiple Targets in OFDM-Based Radar Systems," 
% in IEEE Transactions on Communications, vol. 71, no. 8, pp. 4862-4876, 
% Aug. 2023, doi: 10.1109/TCOMM.2023.3280562.
%
% [3] P. Di Viesti, A. Davoli, G. Guerzoni, et al. "Novel Methods for 
% Approximate Maximum Likelihood Estimation of Multiple Superimposed 
% Undamped Tones and Their Application to Radar Systems," TechRxiv, August 
% 03, 2021, doi: 10.36227/techrxiv.15054321.v2
%
% REVIEW VERSION : V1.0
%
% last review: 12/04/2024

%% EXAMPLE 1: 3D HARMONIC RETRIEVAL PROBLEM

clc, clear, close all
% signal generation --> SIGNAL MODEL: H[m,n,q] = \sum_{l=0}^{L-1} A_l exp(1i*2*pi*m*fD_l)
% exp(-1i*2*pi*n*fr_l) exp(-1i*2*pi*n*f_theta_l) + w[m,n,q]
M = 64;
N = 128;
Nant = 16;
K = 4;

m = 0:M-1;
n = 0:N-1;
nt = 0:Nant-1;
A = [1;0.85;0.75;0.6];
Fm_spacing = 3;                                                                      % normalized Fm spacing
Fn_spacing = 2.5;                                                                    % normalized Fn spacing
Fm = 2.3/M + Fm_spacing*(0:K-1).'/M;                                                 % vector Fm
Fn = 1.8/N + Fn_spacing*(0:K-1).'/N;                                                 % vector Fn
F_theta = [1.5/Nant; 3.2/Nant; 4/Nant; 7.4/Nant];                                    % vector F_theta
Hi = squeeze(sum(A.*exp(1i*2*pi*m.*Fm).*exp(-1i*2*pi*permute(n,[1 3 2]).*Fn) ...
    .*exp(-1i*2*pi*permute(nt,[1 4 3 2]).*F_theta),1));                              % noiseless signal (MxNxNant matrix)
H = awgn(Hi,20,"measured");                                                          % noisy signal (e.g., SNR = 20 dB)

% initialize DRAEC3D parameters
LM = 4;                                                                              % oversampling along 1st dimension of H
LN = 4;                                                                              % oversampling along 2nd dimension of H
LNant = 8;                                                                           % oversampling along 3rd dimension of H
it_r1 = 10;                                                                          % number of refinement steps for each 2D tone
it_r2 = 10;                                                                          % number of refinement steps for each single tone
re_est1 = 2;                                                                         % number of 2D tone re-estimation steps
re_est2 = 2;                                                                         % number of single tone re-estimation steps
par = initialize_draec3D_par(M,N,Nant,K,LM,LN,LNant,it_r1,re_est1,it_r2,re_est2);    % initialization of the parameters

% run the DRAEC 3D
T = DRAEC3D(H, par);
disp(strcat("Estimation ended in ",num2str(T.time)," s"));


if T.cl_freq == 0
    % PLOTS
    figure()
    hold on
    scatter3(Fm,Fn,F_theta,25*A,'filled',"red",'DisplayName','Ground Truth')
    scatter3(T.fD,T.fr,T.f_theta,25*abs(T.A),"cyan",'DisplayName','Estimates')
    title('Estimation accuracy')
    xlabel('$F_M$ [.]','Interpreter','latex')
    ylabel('$F_N$ [.]','Interpreter','latex')
    zlabel('$F_{\theta}$ [.]','Interpreter','latex')
    legend('Location','northwest')

end


%% SUPPORTING FUNCTIONS

function par = initialize_draec3D_par(M,N,Nant,K,LM,LN,LNant,it_r1,re_est1,it_r2,re_est2)
% check input parameters before creating the par structure
par.ref_ant = floor(Nant/2);
par.det_th2D = 0.15;
par.det_th1D = 0.75;
par.csfdec_settings = initialize_csfdec_par(M,N,K,LM,LN,it_r2,re_est2);
par.csfec_settings = initialize_csfec_par(Nant,K,LNant,it_r1,re_est1);
end

% -------------------------------------------------------------------------

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

% -------------------------------------------------------------------------

function par = initialize_csfec_par(N,K,LN,it_residual,re_est)
% check input parameters before creating the par structure
par.N = N;
par.K = K;
par.LN = LN;
par.N0 = par.N*par.LN;
n0 = (0:1:par.N-1);
par.np = [n0; n0.^2; n0.^3];
par.seq_n = 0:1:par.N0-1;
par.V0 = par.seq_n/par.N0;
par.it_residual = it_residual;
par.re_est = re_est;
end

% -------------------------------------------------------------------------




