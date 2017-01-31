function [omegasort,isort] = lam2omegacrit(lam,bdisplay)
%LAM2OMEGACRIT compute angles on the eigensphere to critical points
%
% Here it is enough to only consider the eigenvalues of moment .
% In other words, we consider a moment tensor in its standard basis.
% By permuting the eigenvalues, we keep the same basis U = (e1,e2,e3),
% so that the arc-distance on the eigensphere is the same as the omega
% angle between two moment tensors.
%
% Excluding the critical angle w0 = 0, there are, in general, four unique
% critical angles for each set of eigenvalues. Here are the exceptions:
%    3 unique angles:   lune longitude gamma = 0
%    1 unique angle:    lune boundary
%    0 angle:           isotropic points
%
% See also omegacritmod.m
%
% INPUT
%   lam         3 x n matrix of eigenvalues
%   bdisplay    =true to display details
%
% OUTPUT        
%   omega       n x 5 matrix of sorted omega angles to the 5 other critical
%                   points on the eigensphere
%                   NOTE: omega(:,3) = omega(:,4) ALWAYS
%   isort       n x 5 matrix of indices that show how omega was sorted from
%                   a predetermined order of permuted eigenvalues
%
% Carl Tape, 14-Feb-2015
%

if nargin==1, bdisplay=false; end

[~,n] = size(lam);

% put lambda in the fundamental lune
lam = sort(lam,'descend');

% normalize eigenvalues
lamlen = sqrt( lam(1,:).^2 + lam(2,:).^2 + lam(3,:).^2 );
lam = lam ./ repmat(lamlen,3,1);

% permute eigenvalues
lam1 = lam;                 % fundamental lune
lam2 = lam([1 3 2],:);
lam3 = lam([2 3 1],:);
lam4 = lam([2 1 3],:);
lam5 = lam([3 1 2],:);
lam6 = lam([3 2 1],:);

% calculate omega (arc distances) to critical points
omega = NaN(n,5);
omega(:,1) = getomega(lam1,lam2);
omega(:,2) = getomega(lam1,lam3);
omega(:,3) = getomega(lam1,lam4);
omega(:,4) = getomega(lam1,lam5);
omega(:,5) = getomega(lam1,lam6);

% sort critical angles
[omegasort,isort] = sort(omega,2,'ascend');
%omega, omegasort, isort

if bdisplay
   disp('note that these will be ordered DIFFERENTLY from the omega returned by lam2omegacrit.m:');
   for ii=1:n
       disp(sprintf('===== lune point %i / %i ========',ii,n));
       disp(sprintf('lam1 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam1(:,ii),0));
       disp(sprintf('lam2 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam2(:,ii),omega(ii,1)));
       disp(sprintf('lam3 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam3(:,ii),omega(ii,2)));
       disp(sprintf('lam4 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam4(:,ii),omega(ii,3)));
       disp(sprintf('lam5 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam5(:,ii),omega(ii,4)));
       disp(sprintf('lam6 = (%6.3f, %6.3f, %6.3f) = omega = %8.2f',lam6(:,ii),omega(ii,5)));
   end
   [gamma,delta] = lam2lune(lam);
   disp('sorted omega angles to critical points (2,3,4,5,6) for each input lune point');
   for ii=1:n
       disp(sprintf('lune point %3i: (%6.3f, %6.3f, %6.3f) (%5.1f, %5.1f) : omegas = %6.2f %6.2f %6.2f %6.2f %6.2f',...
           ii,lam1(:,ii),gamma(ii),delta(ii),omegasort(ii,:)));
   end
end

%--------------------------------------

function w = getomega(lam1,lam2)

[~,n] = size(lam1);

% note: lam are already normalized
cosom = zeros(n,1);
for ii=1:n
   cosom(ii) = dot(lam1(:,ii),lam2(:,ii));
end
% correction for possible numerical errors
% this correction is needed for comparing U's that are very close to each other
ipos = cosom > 1;
ineg = cosom < -1;
if sum(ipos), disp(sprintf('%i/%i dot products > 1 (max %.3e)',sum(ipos),n,max(cosom(ipos))-1 )); end
if sum(ineg), disp(sprintf('%i/%i dot products < 1 (max %.3e)',sum(ineg),n,min(cosom(ipos))+1 )); end
cosom(ipos) = 1;
cosom(ineg) = -1;
w = acos(cosom) * 180/pi;

%==========================================================================
% EXAMPLE

if 0==1
    % EXAMPLE 1: simple case
    clear, close all, clc
    n = 7;
    lam = randi([-9 9],3,n)
    [omega,isort] = lam2omegacrit(lam,true);

    % EXAMPLE 2: deviatoric points
    gamma = [-29 -25:5:0];
    delta = zeros(size(gamma));
    lam = lune2lam(gamma,delta);
    [omega,isort] = lam2omegacrit(lam,true);
    
    %% EXAMPLE 2b: fixed lune latitude
    dgam = 10;
    gamma = [-30:dgam:0];
    n = length(gamma);
    delta = 0*ones(n,1);    % fixed
    lam = lune2lam(gamma,delta);
    omega = lam2omegacrit(lam);
    ncrit = NaN(n,1);
    for ii=1:n
        omegamod = omegacritmod(omega(ii,:))
        ncrit(ii) = length(omegamod);
    end
    
    %% EXAMPLE 2c: fixed lune longitude
    ddelta = 30;
    delta = [0:ddelta:90];
    n = length(delta);
    gamma = 0*ones(n,1);   % fixed
    lam = lune2lam(gamma,delta);
    omega = lam2omegacrit(lam);
    ncrit = NaN(n,1);
    for ii=1:n
        omegamod = omegacritmod(omega(ii,:))
        ncrit(ii) = length(omegamod);
    end
    
    %% EXAMPLE 3: omega for each critical point on the lune
    % get a grid of regular cartesian points
    dg = 1;
    gvec = [-30:dg:30];
    %bvec = [-90+dg:dg:90-dg];
    bvec = [-90:dg:90];
    [G,B] = meshgrid(gvec,bvec);
    gamma = G(:);
    delta = B(:);
    lam = lune2lam(gamma,delta);
    n = length(gamma);
    omega = lam2omegacrit(lam);

    % add a set of zeros to represent the first critical point
    omega = [zeros(n,1) omega];

    % plot
    figure; nr=2; nc=3; msize=3^2; cinc=15;
    for kk=1:6
        W = reshape(omega(:,kk),size(G));
        subplot(nr,nc,kk); hold on;
        %scatter(gamma,delta,msize,omega(:,kk),'filled');
        pcolor(G,B,W); shading flat;
        contour(G,B,W,[cinc:cinc:(180-cinc)]','k','linewidth',1);
        axis equal, axis([-30 30 -90 90]); caxis([0 180]);
        title({sprintf('critical point #%i',kk),
            sprintf('\\omega = [%.0f to %.0f]',min(W(:)),max(W(:)))});
    end
    orient tall; wysiwyg;
    %print(gcf,'-dpng','/home/carltape/lam2omegacrit');
end

%==========================================================================
