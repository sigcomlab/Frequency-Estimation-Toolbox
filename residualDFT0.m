%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This algorithm is required by the sfe/csfe
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

