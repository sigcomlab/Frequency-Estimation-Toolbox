%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Tis function is required by the sfe/csfe algorithm
%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F0, F1, F2] = residualDFT_new(q, const, to_lim)

N = const.N;
N2 = N^2;
% N3 = N^3;

% qa = q.^(N - 1);        %q^(N-1)
qb = q.^(N);           %q^(N)
qc = q.^(N + 1);           %q^(N+1)
qd = q.^(N + 2);           %q^(N+2)
% qe = q.^(N + 3);           %q^(N+2)

q_1 = q - 1;
den0 = q_1;
den1 = q_1.^2;
den2 = q_1.^3;
% den3 = q_1.^4;

F0 = (qb - 1) ./ den0;
F1 = ((N - 1) * qc - N * qb + q) ./ den1;
F2 = (((N - 1)^2) * qd + (-2*N2 + 2*N +1)*qc + N2*qb - q.^2 - q) ./ den2;
% F3 = (((N - 1)^3) * qe + (-3*N3 + 6*N2 -4)*qd + ...
%     (3*N3-3*N2-3*N-1)*qc -N3*qb +q.^3+4*q.^2+q) ./ den3;

if any(any(to_lim))
    F0(to_lim) = const.N;
    F1(to_lim) = (const.N^2 - const.N)/2;
    F2(to_lim) = ((2*const.N-1)*(const.N-1)*const.N)/6;
%     F3(to_lim) = ((const.N-1)^2*const.N^2)/4;
end

F0 = F0./ const.N;
F1 = F1./ const.N;
F2 = F2./ const.N;
% F3 = F3./ const.N;
end

