function x = randomint(xmin, xmax, n)
% generates a random integer between xmin and xmax
%
% EXAMPLE:
%   x = randomint(-3, 21, 1000); figure; hist(x); 
%   title(['min = ' num2str(min(x)) ', max = ' num2str(max(x))]);
%
% =======OBSOLETE -- JUST USE MATLAB'S randi=============
%

x0 = rand(n,1);
x = round((xmax - xmin)*x0 + xmin);
