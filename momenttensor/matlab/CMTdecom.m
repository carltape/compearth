function [lam,U] = CMTdecom(M,isort)
%CMTDECOM decompose a set of moment tensors into eigenvalues + basis
%
% INPUT
%   M       6 x n moment tensors with some unspecified basis (e.g., up-south-east)
%           M = [M11 M22 M33 M12 M13 M23]
%   isort   optional: sorting of eigenvalues (default=1)
%
% OUTPUT
%   lam     3 x n set of eigenvalues
%   U       3 x 3 x n set of bases (SAME BASIS AS THE INPUT M)
%
% Inverse program to CMTrecom.m
%
% calls Mdim.m
%
% Carl Tape, 01-Nov-2010

% get deviatoric part
%[Miso,Mdev] = CMTdecom_iso(M);

% make sure M is 6 x n
[M,n] = Mdim(M);

M11 = M(1,:); M22 = M(2,:); M33 = M(3,:);
M12 = M(4,:); M13 = M(5,:); M23 = M(6,:);

% compute eigenvalues and orthonormal eigenvectors
% NOTE: lams(M) = lams(Mdev) + lams(Miso)

lam = zeros(3,n);
U = zeros(3,3,n);

% sorting of eigenvalues
% 1: highest to lowest, algebraic: lam1 >= lam2 >= lam3
% 2: lowest to highest, algebraic: lam1 <= lam2 <= lam3
% 3: highest to lowest, absolute : | lam1 | >= | lam2 | >= | lam3 |
% 4: lowest to highest, absolute : | lam1 | <= | lam2 | <= | lam3 |
if nargin==1
    isort = 1;
else
    if ~any(isort==[1:4]), error(sprintf('isort (%i) must be 1,2,3 or 4',isort)); end
end
slabs = {'lam1 >= lam2 >= lam3','lam1 <= lam2 <= lam3',...
    '| lam1 | >= | lam2 | >= | lam3 |','| lam1 | <= | lam2 | <= | lam3 |'};
disp(sprintf('isort = %i: eigenvalues/eigenvectors sorted by %s',isort,slabs{isort}));

for ii = 1:n
    % moment tensor
    Mx = [ M11(ii) M12(ii) M13(ii) ;
             M12(ii) M22(ii) M23(ii) ;
             M13(ii) M23(ii) M33(ii) ];

    % Matlab eigenvalue ordering is lowest to highest in ALGEBRAIC sense
    [V,D] = eig(Mx);
    lams = diag(D);
    %(inv(V)*Mcmt*V - D) / norm(Mcmt)

    % sort eigenvalues -- IMPORTANT for moment tensor conventions,
    % such as the formula for epsilon and M0, or for computing
    % strike, dip, and rake angles (for double couples)
    isign = sign(lams);
    if isort==1
        [lsort, inds] = sort(lams,'descend');
    elseif isort==2
        [lsort, inds] = sort(lams,'ascend');
    elseif isort==3
        [lsort, inds] = sort(abs(lams),'descend');
        lsort = lsort.*isign(inds);
    elseif isort==4
        [lsort, inds] = sort(abs(lams),'ascend');
        lsort = lsort.*isign(inds);
    end
    Vsort = V(:,inds);      % sort eigenvectors
    %(inv(Vsort)*Mcmt*Vsort - diag(lsort)) / norm(Mcmt)  % check

    lam(:,ii) = lsort;
    U(:,:,ii) = Vsort;
end

% ensure that U is right-handed
U = Udetcheck(U);

%==========================================================================

if 0==1
    clear, close all, clc
    % generate a set of random moment tensors (but all same magnitude)
    npt = 1e3;      % number of randomly generated points (MTs)
    M0ref = 1e16;   % constant M0 for all events
    iarg = [M0ref npt];
    M = CMTspace(iarg,1);
    [lam,U] = CMTdecom(M);
end
    
%==========================================================================