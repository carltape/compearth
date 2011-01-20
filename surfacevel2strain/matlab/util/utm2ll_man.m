%
% function [x,y] = utm2ll_man(xi,yi,i_zone,i_type) 
%
% "Manual" version of utm2ll, where the source code is explicit.
% Compare with the Matlab blackbox version in utm2ll.m.
%
% calls utm2ll_src.m
% called by xxx
%

function [x,y] = utm2ll_man(xi,yi,i_zone,i_type) 

[a,b] = size(xi);
num = a*b;
xvec0 = reshape(xi,num,1);
yvec0 = reshape(yi,num,1);

for ii=1:length(xi)
    if mod(ii,1000) == 0, disp(sprintf('%i out of %i',ii,length(xi))); end
    [xvec(ii),yvec(ii)] = utm2ll_src(xvec0(ii),yvec0(ii),i_zone,i_type);
end
x = reshape(xvec,a,b);
y = reshape(yvec,a,b);
        
%==========================================================================
