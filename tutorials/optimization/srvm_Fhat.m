function Fhat = srvm_Fhat(Shat0,niter,nu,a,w)
%SRVM_FHAT construct matrix Fhat by repeatedly calling srvm_Shat_chi.m
%
% See srvm_Shat_chi.m for details.
%
% INPUT:
%   S0hat   nparm x nparm initial matrix
%   nu, a   niter x 1 vectors of stored scalar values
%   w       nparm x niter matrix of stored vectors
%
% OUTPUT:
%   Fhat    nparm x nparm matrix ( F = S * Shat' * icprior )
%
% Carl Tape, 05-June-2007
%

[~,nparm] = size(Shat0);

Fhat = zeros(nparm,nparm);
for ii = 1:nparm
    chi = zeros(nparm,1);
    chi(ii) = 1;
    ShatT_chi  = srvm_Shat_chi(chi,niter,Shat0,nu,a,w,1);
    Fhat_chi   = srvm_Shat_chi(ShatT_chi,niter,Shat0,nu,a,w,0);
    Fhat(:,ii) = Fhat_chi;
end

%==========================================================================