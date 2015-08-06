function plotMT_dc(kap1,theta1,sig1,kap2,theta2,sig2)
%PLOTMT_DC plot histograms of strike-dip-rake angles for both fault planes
%
% INPUT
%    kap1       strike angle, deg (0,360)
%    theta1     dip angle, deg (0,-90)
%    sig1       rake angle, deg (-180,180)
%
% The second set of fault angles is optional.
% You can also use dcfaultpar2nodal.m to get the complementary set of angles.
%
% See CMT2dcfaultpar.m and CMT2TT.m
% 

deg = 180/pi;
figure; nr=3;
if nargin==3
    nc = 1;
else
    nc = 2;
end

% number of bins in each histogram
nbin = 36;
nedge = nbin + 1;
hflat = 1/nbin;

kapbin = linspace(0,360,nedge);
thetabin = linspace(0,90,nedge);
sigbin = linspace(-180,180,nedge);

ktick = [0:60:360];
dtick = [0:15:90];
ltick = [-180:60:180];

if nc==2
    subplot(nr,nc,1); plot_histo(kap1,kapbin); xlabel('strike-1, deg'); set(gca,'xtick',ktick);
    subplot(nr,nc,2); plot_histo(kap2,kapbin); xlabel('strike-2, deg'); set(gca,'xtick',ktick);
    subplot(nr,nc,3); plot_histo(theta1,thetabin); xlabel('dip-1, deg'); set(gca,'xtick',dtick);
    subplot(nr,nc,4); plot_histo(theta2,thetabin); xlabel('dip-2, deg'); set(gca,'xtick',dtick);
    subplot(nr,nc,5); plot_histo(sig1,sigbin); xlabel('rake-1, deg'); set(gca,'xtick',ltick);
    subplot(nr,nc,6); plot_histo(sig2,sigbin); xlabel('rake-2, deg'); set(gca,'xtick',ltick);
else
    subplot(nr,nc,1); hold on; plot_histo(kap1,kapbin); xlabel('strike, deg'); set(gca,'xtick',ktick);
    plot([0 360],hflat*[1 1],'r--','linewidth',2);
    
    % the dip is NOT flat for uniform distribution
    subplot(nr,nc,2); hold on; plot_histo(theta1,thetabin); xlabel('dip, deg'); set(gca,'xtick',dtick);
    acont = linspace(0,90);
    nmax = hflat*pi/2;
    plot(acont,nmax*sin(acont/deg),'r--','linewidth',2);
    
    subplot(nr,nc,3); hold on; plot_histo(sig1,sigbin); xlabel('rake, deg'); set(gca,'xtick',ltick);   
    plot([-180 180],hflat*[1 1],'r--','linewidth',2);
end

%==========================================================================
