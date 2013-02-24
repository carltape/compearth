function [xvec,yvec,numy,a,b,dx] = gridvec(xmin,xmax,numx,ymin,ymax)
%GRIDVEC returns two vectors representing a uniform 2D grid of points
%
% This function inputs specifications for creating a grid
% of uniformly spaced points, reshaped into column vectors
% of the x- and y-coordinates.  Note that dx = dy.

xvec0 = linspace(xmin,xmax,numx);
dx = xvec0(2) - xvec0(1);
yvec0 = [ymin : dx : ymax];

[X,Y] = meshgrid(xvec0,yvec0);
[a,b] = size(X);
xvec = reshape(X,a*b,1);
yvec = reshape(Y,a*b,1);
numy = length(yvec0);