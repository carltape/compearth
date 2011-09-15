%
% function x = randomvec(xmin, xmax, n)
%
% This generates a vector of n random numbers in the range [xmin,xmax].
%
% Is there a matlab built-in function that does this?
% There is randi for integers, but what about for non-integers?
%

function x = randomvec(xmin0, xmax0, n)

% flip input order if needed
xmin = min([xmin0 xmax0]);
xmax = max([xmin0 xmax0]);

x = (xmax - xmin)*rand(n,1) + xmin;

%================================================
