%
% function [xvec,yvec] = gridvec(xmin,xmax,numx,ymin,ymax)
%
% This function inputs specifications for creating a grid
% of uniformly spaced points, reshaped into column vectors
% of the x- and y-coordinates.  Note that dx = dy.
%
% calls xxx
% called by xxx
%

function [xvec,yvec] = gridvec(xmin,xmax,numx,ymin,ymax)

xvec0 = linspace(xmin,xmax,numx);
dx = xvec0(2) - xvec0(1);
yvec0 = [ymin : dx : ymax];

[X,Y] = meshgrid(xvec0,yvec0);
[a,b] = size(X);
xvec = reshape(X,a*b,1);
yvec = reshape(Y,a*b,1);

%====================================================