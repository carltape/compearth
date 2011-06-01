%
% function theta = fangle(va,vb)
% Carl Tape, 12-April-2011
%
% INPUT
%   va      3 x n set of initial vectors
%   vb      3 x n set of rotated vectors
%
% OUTPUT
%   theta   n x 1 set of rotation angles
% 
%
% calls flength.m
% called by xxx
%

function theta = fangle(va,vb)

% check dimensions
[n1a,n2a] = size(va);
[n1b,n2b] = size(vb);
%if any([n1a n1b] ~= 3), error('va and vb must be 3 x n'); end
%if n2a ~= n2b
%    error('va and vb must be 3 x n');
%else
%    n = n2a;
%end
if or(n2a ~= n2b, n1a ~= n1b)
    error('va and vb must be same dimension');
else
    n = n2a;
end

ma = flength(va);
mb = flength(vb);
vadvb = va(1,:).*vb(1,:) + va(2,:).*vb(2,:) + va(3,:).*vb(3,:);
theta = 180/pi * acos( vadvb ./ (ma .* mb ) );

%--------------
% EXAMPLES

if 0==1
    va = [1 0 0]'; vb = [0 0 3]'; theta = fangle(va,vb)
    va = [1 0 0 0 0 0]'; vb = [0 0 0 0 0 3]'; theta = fangle(va,vb)
end

%===========================================================
