%
% function [xvec,yvec] = gridvec(xmin,xmax,numx,ymin,ymax)
% CARL TAPE, 23-March-2005
% printed xxx
%
% This function inputs specifications for creating a grid
% of uniformly spaced points, reshaped into column vectors
% of the x- and y-coordinates.  Note that dx = dy.
%
% See griddataXB.m
%
% calls xxx
% called by wave2d_basis_fun.m, spline_wang_D.m
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