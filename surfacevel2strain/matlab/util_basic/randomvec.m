function x = randomvec(xmin0,xmax0,n)
%RANDOMVEC generate a vector of n random numbers in the range [xmin,xmax]
%
% I'm not aware of a built-in matlab function that does this.
% There is randi for integers, but what about for non-integers?

% flip input order if needed
xmin = min([xmin0 xmax0]);
xmax = max([xmin0 xmax0]);

x = (xmax - xmin)*rand(round(n),1) + xmin;