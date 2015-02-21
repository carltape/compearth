function [omega,xi0,U] = CMT2omega_xi0(X1,X2,iorthoU,idisplay)
%CMT2OMEGA_XI0 compute the omega and xi0 angles between two moment tensors
%
% INPUT
%   X1,X2     can be either U1,U2 or M1,M2 (see below)
%   iorthoU   =0 to NOT orthogonalize; >1 to orthogonalize (see Uorth.m) 
%             if X1=M1 and X2=M2, then iorthoU=0 should always be okay
%   idisplay  optional: =1 to display details
%
%   M1,M2     6 x n moment tensors: M = [M11 M22 M33 M12 M13 M23]
%   U1,U2     3 x 3 x n bases
%
% OUTPUT
%   omega     9-angle between closest double couples
%   xi0       minimum rotation angle between principal axes
%   U         3 x 3 x n rotation matrix U = U1' * U2
%
% For details, see TapeTape2012 "Angle between principal axis triples".
%
% EXAMPLES: see below and also TT2012kaganAppE.m
% Set bfigure=true to plot histograms of the distributions.
% 
% Carl Tape 8/11/2012

bfigure = true;

% default: no information displayed
if ~exist('idisplay','var'), idisplay = 0; end 

if and(length(size(X1))==2,numel(X1)~=9)
    disp('CMT2omega_xi0.m: input arrays are moment tensors');
    M1 = X1; M2 = X2;
    
    % check that M is 6 x n
    [M1,n1] = Mdim(M1);
    [M2,n2] = Mdim(M2);
    if n1~=n2, error('M1 and M2 must be equal in dimension'); end
    n = n1;

    % get basis
    % note 1: by design, M represents a truly symmetric matrix;
    %         therefore the basis will be exactly orthogonal
    % note 2: no need for eigenvalues, though it is key that the
    %         columns of U are sorted in order lam1 >= lam2 >= lam3
    [~,U1] = CMTdecom(M1);
    [~,U2] = CMTdecom(M2);

else
    disp('CMT2omega_xi0.m: input arrays are bases');
    U1 = X1; U2 = X2;
    % check dimensions
    [~,~,n1] = size(U1);
    [~,~,n2] = size(U2);
    if n1~=n2, error('U1 and U2 must be equal in dimension'); end
    n = n1;
end

% ensure that U are rotation matrices: det U = 1
% note: this may already be done in CMTdecom.m
U1 = Udetcheck(U1);
U2 = Udetcheck(U2);

if iorthoU > 0
    disp('orthogonalizing U');
    U1 = Uorth(U1,iorthoU);
    U2 = Uorth(U2,iorthoU);
end
  
% OMEGA
% find angle between closest double couples
lam0 = repmat([1 0 -1]',1,n);
MDC1 = CMTrecom(lam0,U1);
MDC2 = CMTrecom(lam0,U2);
omega = CMT2omega(MDC1,MDC2);

% XI
% compute U = U1'*U2
U = UiU(U1,U2);
xi0 = U2xi0(U,0,idisplay);

disp(sprintf('%i/%i xi0 values are imaginary',n-sum(~imag(xi0)),n));
disp(sprintf('%i/%i omega values are imaginary',n-sum(~imag(omega)),n));

if idisplay==1
    disp('CMT2omega_xi0.m all omega and xi0:');
   for ii=1:n
      disp(sprintf('%12i/%6i: omega = %8.4f, xi0 = %8.4f',ii,n,omega(ii),xi0(ii))); 
   end
end

if and(bfigure, n>1)
   OMAX = 0.07;  % will depend on bin size
   figure; nr=2; nc=1;
   subplot(nr,nc,1); hold on; plot_histo(omega,[0:5:180]);
   plot([90 90],[0 OMAX],'r','linewidth',2);
   set(gca,'xtick',0:10:180); %axis([0 120 0 OMAX]);
   xlabel('omega angle describing difference in orientation');
   title(sprintf('OMEGA: min = %.2f, max = %.2f',min(omega),max(omega)));
   
   PMAX = 0.12;  % will depend on bin size
   subplot(nr,nc,2); hold on; plot_histo(xi0,[0:5:120]);
   plot([90 90],[0 PMAX],'r','linewidth',2);
   set(gca,'xtick',0:10:120); %axis([0 120 0 PMAX]);
   xlabel('xi0 angle describing difference in orientation');
   title(sprintf('XI: min = %.2f, max = %.2f',min(xi0),max(xi0)));
end

%==========================================================================
% EXAMPLES

if 0==1
    % Appendix E of TapeTape2012 "Angle betweeen principal axis triples"
    % xi0 = 102.5 deg as in Eq E1
    % see TT2012kaganAppE.m for details
    M1 = [-0.4120    0.0840    0.3280    0.3980   -1.2390    1.0580]'*1e19;
    M2 = [ 5.0540   -2.2350   -2.8190   -0.4760    5.4200    5.5940]'*1e18;
    [omega,xi0] = CMT2omega_xi0(M1,M2,0,1); omega, xi0
end

%==========================================================================