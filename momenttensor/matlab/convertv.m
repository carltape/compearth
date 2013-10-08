function [vout,T] = convertv(i1,i2,v)
%CONVERTV convert vectors among different bases
%
% INPUT
%   i1      index of input moment tensor basis (see convert_MT.m)
%   i2      index of output moment tensor basis (see convert_MT.m)
%   V       3 x n set of vectors (or 3 x 3 x n set of bases)
%
% OUTPUT
%   vout    3 x n set of vectors in basis of i2 (or 3 x 3 x n set of bases in i2)
%   T       transformation matrix to change basis of v from i1 to i2:
%              vout = T*v
%
% See also convert_MT.m
%
% Carl Tape, 10/2013
%

% get transformation matrix
T = convert_MT(i1,i2);

[a,b,c] = size(v);
if c==1
    if a~=3, v = v'; n = b; disp('convertv.m: taking transpose of input v'); else n = a; end
    vout = T*v;
else
    if and(a==3,b==3)
        vout = NaN(a,b,c); 
        for ii=1:c
           vout(:,:,ii) = T*v(:,:,ii);
        end
    else
        whos v
        disp('check dimensions of input v'); 
    end
end
    
%==========================================================================
% EXAMPLES

if 0==1
    % input standard basis for up-south-east
    % in this case we convert from up-south-east to north-east-down
    % (1,0,0) points up    for v and is (0,0,-1) for vout
    % (0,1,0) points south for v and is (-1,0,0) for vout
    % (0,0,1) points east  for v and is (0,1,0) for vout
    i1 = 1; i2 = 2;
    v = eye(3,3)
    [vout,T] = convertv(i1,i2,v)
    
    % what about the other way?
    % in this case we convert from north-east-down to up-south-east
    % (1,0,0) points north for v and is (0,-1,0) for vout
    % (0,1,0) points east  for v and is (0,0,1) for vout
    % (0,0,1) points down  for v and is (-1,0,0) for vout
    i1 = 2; i2 = 1;
    v = eye(3,3)
    [vout,T] = convertv(i1,i2,v)
    
    % check that converting the basis vectors is the same as converting
    % the moment tensor
    i1 = 1; i2 = 2;
    M = rand(6,1);
    Mmat = Mvec2Mmat(M,1);
    [U,D] = eig(Mmat);
    % convert the basis vectors
    [Uout,T] = convertv(i1,i2,U)
    Mcheck = Uout*D*Uout'
    % convert the moment tensor
    Mout = convert_MT(i1,i2,M);
    Mvec2Mmat(Mout,1)
end


%==========================================================================

