%
% function mm = counter2(nr,nc)
%
% This program is used for subplotting indexing.  Given the subplot
% dimension (nr,nc), this gives a vector that orders the subplots counting
% BY COLUMN, rather than by row.
%
% For example, if you have a 4 x 2 plotting figure, then you would index
% mm = counter2(4,2) = [1 3 5 7 2 4 6 8]
%   |  1  2  |
%   |  3  4  |
%   |  5  6  |
%   |  7  8  |
%
% calls xxx
% called by xxx
%

function mm = counter2(nr,nc)

kk = 1;
for jj=1:nc
    for ii=1:nr
        mm(kk) = nc*(ii-1) + jj;
        %disp([ ii jj mm(kk)]);
        kk = kk+1;
    end
end

%===================================================
