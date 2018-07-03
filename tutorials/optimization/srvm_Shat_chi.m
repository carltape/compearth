function Shat_chi = srvm_Shat_chi(chi,niter,S0hat,nu,a,w,itranspose)
%SRVM_SHAT_CHI compute Shat*chi or Shat'*chi from a set of w vectors and nu/a scalars
%
% Square-root variable metric method
%
% INPUT:
%   chi         nparm x 1 input vector
%   S0hat       nparm x nparm initial matrix
%   nu, a       niter x 1 vectors of stored scalar values
%   w           nparm x niter matrix of stored vectors
%   itranspose  =1 for Shat'*chi; =0 for Shat*chi
%
% OUTPUT:
%   Shat_chi    nparm x 1 vector given by Shat*chi or Shat'*chi
%
% called by srvm_Fhat.m
%
% Carl Tape, 05-June-2007
%

%niter = length(nu);

if itranspose == 1
    % operate Shat'*chi using vector operations
    
    Shat_chi = S0hat'*chi;    % initialize
    
    for jj = 1 : 1 : niter
        wtemp = w(:,jj);                                 % vector
        xtemp = wtemp'*Shat_chi;                         % scalar
        Shat_chi = Shat_chi - xtemp*nu(jj)/a(jj)*wtemp;  % vector
    end

elseif itranspose == 0
    % operate Shat*chi using vector operations
    
    Shat_chi = chi;    % initialize
    
    for jj = niter : -1 : 1
        %irev = niter-jj+1;
        wtemp = w(:,jj);
        xtemp = wtemp'*Shat_chi;
        Shat_chi = Shat_chi - xtemp*nu(jj)/a(jj)*wtemp;
    end
    Shat_chi = S0hat * Shat_chi;
    
else
    error('Error: choose itranspose = 1 or 0');  
end

% function Fhat_chi = srvm_Fhat_chi(chi,S0hat,nu,a,w)
% 
% niter = length(nu);
% 
% % initialize Shat_chi
% Shat_chi = S0hat'*chi;
% 
% % operate Shat'*chi using vector operations
% for jj = 1:niter
%     wtemp = w(:,jj);                                 % vector
%     xtemp = wtemp'*Shat_chi;                         % scalar
%     Shat_chi = Shat_chi - xtemp*nu(jj)/a(jj)*wtemp;  % vector
% end
% 
% % operate S*Shat'*chi using vector operations
% Fhat_chi = Shat_chi;
% for jj = 1:niter
%     irev = niter-jj+1;
%     wtemp = w(:,irev);
%     xtemp = wtemp'*Fhat_chi;
%     Fhat_chi = Fhat_chi - xtemp*nu(irev)/a(irev)*wtemp;
% end
% 
% Fhat_chi = S0hat * Fhat_chi;

%==========================================================================