function [inds,xval,yval] = getsubset(x,y,ax0)
%GETSUBSET get a subset of x-y points within ax0 = [xmin xmax ymin ymax]
%
% input x and y can be vectors or matrices

inds = find( and( and( x >= ax0(1), x <= ax0(2)), and(y >= ax0(3), y <= ax0(4))) );
if nargout==3
    xval = x(inds);
    yval = y(inds);
end

disp(sprintf('getsubset.m: %i points in the subset out of %i',length(inds(:)),length(x(:))));

% % correct if xmin and xmax are reversed
% if xmin > xmax
%     xtemp = xmax;
%     xmax  = xmin;
%     xmin  = xtemp;
% end
% 
% if xmin <= min(x)
%     iq1 = 1;
% else
%     itemp = find(abs(x - xmin) < trsh);
%     iqt1 = round(length(itemp)/2);
%     iq1 = itemp(iqt1);
% end
% 
% if xmax >= max(x)
%     iq2 = length(x);
% else    
%     itemp = find(abs(x - xmax) < trsh);     % find nearest x-values
%     iqt2 = round(length(itemp)/2);          % take the middle index
%     iq2 = itemp(iqt2);
% end
% 
% iran = [iq1 : iq2];
