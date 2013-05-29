function [T,k] = lam2Tk(lam)
%LAM2TK input eigenvalues of moment tensor, output T-k values of Hudson et al. (1989)
%
% INPUT
%   lam     3 x n set of eigenvalues
% OUTPUT
%   T       CLVD measure of Hudson et al. (1989)
%   k       ISO measure of Hudson et al. (1989)
%
% NOTE: The input eigenvalues do NOT need to be sorted in any order.
%
% Carl Tape, 12/2012

if numel(lam)==3, n=1; else [~,n] = size(lam); end
lamI = 1/3*sum(lam);

% deviatoric eigenvalues
lamd = lam - repmat(lamI,3,1);

% T: CLVD component
% Julian et al. (1998), Eq. (18)
% (see also CMT2epsilon.m)
[~,imat] = sort(abs(lamd),'descend');
epsilon = zeros(1,n);
for ii=1:n
    i1 = imat(3,ii);   % smallest in absolute sense
    i2 = imat(1,ii);   % largest in absolute sense
    if all(lamd(:,ii)==0)
        disp('WARNING: input is isotropic -- we set epsilon=0');
    else
        epsilon(ii) = -lamd(i1,ii) / abs( lamd(i2,ii) );
    end
end
T = -2*epsilon;

% k: isotropic component
% Ford et al. (2009), Baig and Urbancic (2010)
lammaxd = max(abs(lamd));

% % Julian et al. (1998), Eq. (19) -- THIS EQUATION IS AMBIGUOUS
% % SHOULD WE GET M3, THEN TAKE THE DEVIATORIC PART (AS HERE),
% % OR SHOULD WE TAKE THE DEV PART, THEN GET MDEV3?
% % maximum (in the absolute sense) eigenvalue -- can be positive or negative
% lammax = zeros(1,n);
% for ii=1:n
%    [~,imax] = max(abs(lam(:,ii)));
%    lammax(ii) = lam(imax,ii);
% end
% lammaxd = lammax - lamI;

% other unclear equations for k:
% Hudson et al. (1989), Foulger et al. (2004)

k = lamI ./ ( abs(lamI) + abs(lammaxd) );

%==========================================================================

if 0==1
    % get a grid of regular cartesian points
    dg = 1;
    gvec = [-30:dg:30];
    bvec = [-90+dg:dg:90-dg];
    [G,B] = meshgrid(gvec,bvec);
    gamma = G(:);
    delta = B(:);
    n = length(gamma);
    M0 = 1/sqrt(2)*ones(n,1);
    lam = lune2lam(gamma,delta,M0);
    
    [T,k] = lam2Tk(lam);
    figure; nr=1; nc=2; msize = 4^2;
    subplot(nr,nc,1); scatter(gamma,delta,msize,T,'filled'); colorbar, axis equal, axis tight;
    xlabel('gamma'); ylabel('delta'); title('T from T-k plots');
    subplot(nr,nc,2); scatter(gamma,delta,msize,k,'filled'); colorbar, axis equal, axis tight;
    xlabel('gamma'); ylabel('delta'); title('k from T-k plots');
    
    % simple test
    T = 0.5, k = -0.3,
    lam = Tk2lam(T,k)
    [Tcheck,kcheck] = lam2Tk(lam)
end

%==========================================================================
