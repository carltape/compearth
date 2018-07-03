function F = vm_F(F0,icprior,u,v)
%VM_F construct F by repeatedly calling vm_F_chi.m
%
% See vm_F_chi.m for details.
%
% INPUT:
%    F0      nparm x nparm initial preconditioner
%    icprior nparm x nparm prior covariance matrix
%    u       nparm x niter matrix of nparm x 1 vectors
%    v       niter x 1 vector of stored values
%
% OUTPUT:
%    F      nparm x nparm preconditioner that estimates the curvature H^(-1)
%
% calls vm_F_chi.m
%
% Carl Tape, 05-June-2007
%

[~,nparm] = size(F0);

F = zeros(nparm,nparm);
for ii = 1:nparm
    chi = zeros(nparm,1);
    chi(ii) = 1;
    Fchi = vm_F_chi(chi,F0,icprior,u,v);
    F(:,ii) = Fchi;
end

%=======================================================================