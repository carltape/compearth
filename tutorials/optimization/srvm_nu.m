function nu = srvm_nu(a,b)
%SRVM_NU square-root variable metric: compute scalar nu from scalars a and b
%
% see formulas of Williamson (1975) and Hull and Tapley (1977).
%
% The objective is to determine a value for nu that keeps the
% pre-conditioning matrix symmetric positive definite.
%
% Carl Tape, 10-June-2007
%

ratio = b/a;

if ratio < 1
    % Note that two solutions are possible, but the negative sign is
    % preferred, because it makes nu continuous in the region of b/a = 0
    % (Hull and Tapley, 1977).
    %nu = (1 + sqrt(1-ratio)) / ratio;
    nu = (1 - sqrt(1-ratio)) / ratio;

elseif ratio == 1
    nu = 0.9;    % "a value of nu such that nu < 1"

else
    nu = 1.0;
end

%==========================================================================