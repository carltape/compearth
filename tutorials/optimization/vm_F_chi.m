function Fchi = vm_F_chi(chi,F0,icprior,u,v)
%VM_F_CHI variable metric algorithm without matrices
%
% This function is based on Tarantola (2005), Section 6.22, Eq. 6.347.
% It computes the operation of F on an arbitrary vector chi by storing a
% set of vectors (u_k) and scalars (v_k).
%
% EVERYTHING HERE IS ASSUMED TO BE IN THE NONHAT-NOTATION.
%
% INPUT:
%    chi     nparm x 1 vector
%    F0      nparm x nparm initial preconditioner
%    icprior nparm x nparm prior covariance matrix
%    u       nparm x niter matrix of nparm x 1 vectors
%    v       niter x 1 vector of stored values
%
% OUTPUT:
%    Fchi   nparm x 1 vector of F*chi
%
% CARL TAPE, 05-June-2007
%

[~,niter] = size(u);

% compute F*chi
Fchi = F0*chi;
for jj = 1:niter
    vtmp = chi'*icprior*u(:,jj);
    Fchi = Fchi + vtmp/v(jj) * u(:,jj);
end

%==========================================================================