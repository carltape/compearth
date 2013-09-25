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
% EXAMPLES: see below.
% Set ifigure=1 to plot histograms of the distributions.
% 
% Carl Tape 8/11/2012

% default: no information displayed
if ~exist('idisplay','var'), idisplay = 0; end 

if and(length(size(X1))==2,numel(X1)~=9)
    disp('input arrays are moment tensors');
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
    disp('input arrays are bases');
    U1 = X1; U2 = X2;
    % check dimensions
    [~,~,n1] = size(U1);
    [~,~,n2] = size(U2);
    if n1~=n2, error('U1 and U2 must be equal in dimension'); end
    n = n1;
end

% ensure that U are rotation matrices: det U = 1
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
% convert to matrix
M1mat = Mvec2Mmat(MDC1,1);
M2mat = Mvec2Mmat(MDC2,1);
% calculate omega
cosom = zeros(n,1);
for ii=1:n
   M1x =  M1mat(:,:,ii);
   M2x =  M2mat(:,:,ii);
   cosom(ii) = dot(M1x(:),M2x(:));
end
cosom = cosom / 2;  % since |Lam0| = sqrt(2)
% correction for numerical errors
ipos = cosom > 1;
ineg = cosom < -1;
disp(sprintf('%i/%i dot products > 1',sum(ipos),n));
disp(sprintf('%i/%i dot products < 1',sum(ineg),n));
% this correction is needed for comparing U's that are very close to each other
cosom(ipos) = 1;
cosom(ineg) = -1;
omega = acos(cosom) * 180/pi;

% XI
% compute U = U1' * U2
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

ifigure = 1;
if and(ifigure==1, n>1)
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
    clear, close all, clc
    % Appendix E of TapeTape2012 "Angle betweeen principal axis triples"
    % Kagan (1991) comparison events from GCMT catalog
    eid1 = 'C010677A'; eid2 = 'C092680B';   % New Guinea
    fac1 = 1e19; fac2 = 1e18;
    if 0==1     % if you have access to the full GCMT catalog
        [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid,elabel,...
            str1,dip1,rk1,str2,dip2,rk2,lams,pl1,az1,pl2,az2,pl3,az3] = readCMT;
        %eid1 = 'B060586B'; eid2 = 'B062486F';  % south Pacific
        ieid1 = find(strcmp(eid1,eid)==1);
        ieid2 = find(strcmp(eid2,eid)==1);
        % event 1
        disp(sprintf('%s %s',eid{ieid1},datestr(otime(ieid1),29)));
        disp(sprintf('str/dip/rk (1) = %i/%i/%i',str1(ieid1),dip1(ieid1),rk1(ieid1)));
        disp(sprintf('str/dip/rk (2) = %i/%i/%i',str2(ieid1),dip2(ieid1),rk2(ieid1)));
        disp(sprintf('pl1/az1 = %i/%i, pl2/az2 = %i/%i, pl3/ax3 = %i/%i',...
            pl1(ieid1),az1(ieid1),pl2(ieid1),az2(ieid1),pl3(ieid1),az3(ieid1)));
        M1 = M(:,ieid1);
        % event 2
        disp(sprintf('%s %s',eid{ieid2},datestr(otime(ieid2),29)));
        disp(sprintf('str/dip/rk (1) = %i/%i/%i',str1(ieid2),dip1(ieid2),rk1(ieid2)));
        disp(sprintf('str/dip/rk (2) = %i/%i/%i',str2(ieid2),dip2(ieid2),rk2(ieid2)));
        disp(sprintf('pl1/az1 = %i/%i, pl2/az2 = %i/%i, pl3/ax3 = %i/%i',...
            pl1(ieid2),az1(ieid2),pl2(ieid2),az2(ieid2),pl3(ieid2),az3(ieid2)));
        M2 = M(:,ieid2);
    else
        % New Guinea
        M1 = [-0.4120    0.0840    0.3280    0.3980   -1.2390    1.0580]'*1e19;
        M2 = [ 5.0540   -2.2350   -2.8190   -0.4760    5.4200    5.5940]'*1e18;
    end
    % event 1
    disp('///////////////// EVENT 1 /////////////////');
    displayCMTshort(M1/fac1,'%16.3f');
    disp(sprintf('multiplication factor is %.0e N-m',fac1));
    [lam,U] = CMTdecom(M1); U = Udetcheck(U); U14 = Ufour(U);
    % note that FRAME 4 will match the U1 listed in TapeTape Appendix E
    for kk=1:4
        disp(sprintf('========= FRAME %i ===========',kk));
        U1 = U14(:,:,kk);
        U1, det(U1), diag(lam), U1*diag(lam)*U1', Mvec2Mmat(M1,1)
    end
    % event 2
    disp('///////////////// EVENT 2 /////////////////');
    displayCMTshort(M2/fac2,'%16.3f');
    disp(sprintf('multiplication factor is %.0e N-m',fac2));
    [lam,U] = CMTdecom(M2); U = Udetcheck(U); U24 = Ufour(U);
    for kk=1:4
        disp(sprintf('========= FRAME %i ===========',kk));
        U2 = U24(:,:,kk);
        U2, det(U2), diag(lam), U2*diag(lam)*U2', Mvec2Mmat(M2,1)
    end
    
    % now consider the 'difference' matrix U = U1^-1 U2
    U12 = U14(:,:,4)'*U24(:,:,4)   % to match results in paper
    [xi0,ixi0,q] = U2xi0(U12,0,1);     % q1 matches the results in paper
    [omega,xi0] = CMT2omega_xi0(M1,M2,0,1); omega, xi0
end

%==========================================================================