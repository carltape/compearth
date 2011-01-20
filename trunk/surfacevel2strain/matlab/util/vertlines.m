%
% function [xmat, ymat] = vertlines(xvec,ymin,ymax)
%
% This function inputs a vector of x-values and outputs two matrices that
% will plot vertical lines at the x-values specified in xvec.
%
% EXAMPLE:
%   [xmat,ymat] = vertlines(linspace(0,10,6),-1,4); figure; plot(xmat,ymat,'r');
%
% calls xxx
% called by xxx
% 

function [xmat, ymat] = vertlines(xvec,ymin,ymax)

n = length(xvec);
xvec = xvec(:)';        % xvec must be a row
if length([ymin(:) ; ymax(:) ]) ~= 2
    error('ymin and ymax must be single numbers');
end

% switch if incorrectLy entered
if ymin > ymax
    temp = ymin;
    ymin = ymax;
    ymax = temp;
end

xmat = [xvec; xvec];
ymat = [ymin*ones(1,n); ymax*ones(1,n)];

%=====================================================
