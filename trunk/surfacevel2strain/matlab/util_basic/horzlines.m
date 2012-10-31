function [xmat, ymat] = horzlines(yvec,xmin,xmax)
%HORZLINES create matrices for easy plotting of horizontal lines
%
% function [xmat, ymat] = horzlines(yvec,xmin,xmax)
%
% This function inputs a vector of y-values and outputs
% two matrices that will plot horizontal lines
% at the y-values specified in xvec.
%
% EXAMPLE:
%   [xmat,ymat] = horzlines(linspace(0,10,6),-1,4); figure; plot(xmat,ymat,'r');
% 

n = length(yvec);
yvec = yvec(:)';        % must be a row

xmin = xmin(:)';
xmax = xmax(:)';

if length(xmin) == 1, xmin = xmin*ones(1,n); end
if length(xmax) == 1, xmax = xmax*ones(1,n); end

% switch if incorrectLy entered
iflip = find( xmin > xmax ); nflip = length(iflip);
temp = xmin(iflip);
xmin(iflip) = xmax(iflip);
xmax(iflip) = temp;
if nflip > 0, disp('  '); disp([ num2str(nflip) ' horzlines entered in reverse order']); end

xmat = [xmin; xmax];
ymat = [yvec; yvec];

%==========================================================================
