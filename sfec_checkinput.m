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
%
% ** Authors**  
% Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it) 
% Pasquale Di Viesti (pasquale.diviesti@unimore.it)
% Giorgio Guerzoni (giorgio.guerzoni@unimore.it)
% Michele Mirabella (michele.mirabella@unimore.it)
% Elia Vignoli (elia.vignoli@unimore.it)
%
% v1.0 - march 2024
% Copyright (c) Prof. Giorgio Matteo Vitetta (giorgiomatteo.vitetta@unimore.it)
% Website: https://www.sigcom.unimore.it/
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
