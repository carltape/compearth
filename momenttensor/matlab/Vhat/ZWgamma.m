function [Z,W,a,b,c] = ZWgamma(phi,sigma,gamma)
%ZWGAMMA
%
% INPUT
%   mu      in radians
%   nu      in radians
%   gamma   in radians
%

[a,b,c] = abcgamma(phi,sigma,gamma);

Z = -b ./ (2*a);
W = -(b.^2 - 4*a.*c)./(4*a);

if 0==1
   % test calculation
   format long
   deg = 180/pi;
   [Z,W,a,b,c] = ZWgamma(10/deg,15/deg,5/deg)
    
   % plot variations in the functions
   close all
   nphi = 50;
   nsigma = 2*nphi;
   deg = 180/pi;
   phivec = linspace(-pi/2,pi/2,nphi);
   sigmavec = linspace(-pi/2,pi/2,nsigma);
   [Pu,Su] = meshgrid(phivec,sigmavec);
   
   ax0 = [-1 1 -1 1]*90; clims = [0 2]; climsB = [-1 1];
   dticks = [-90:45:90];
   
   gvec = [-25:5:25]/deg;
   for ii = 1:length(gvec);
       gamma = gvec(ii);
       [Zgam,Wgam,Agam,Bgam,Cgam] = ZWgamma(Pu,Su,gamma);
       Pup = Pu*deg;
       Sup = Su*deg;

       figure; nr=3; nc=2;
       subplot(nr,nc,1); pcolor(Pup,Sup,Agam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\phi'); ylabel('\sigma'); title(sprintf('a_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,2); pcolor(Pup,Sup,Bgam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\phi'); ylabel('\sigma'); title(sprintf('b_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,3); pcolor(Pup,Sup,Cgam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\phi'); ylabel('\sigma'); title(sprintf('c_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       
       subplot(nr,nc,5); pcolor(Pup,Sup,Zgam); shading flat; axis equal, axis(ax0); caxis(climsB); colorbar;
       xlabel('\phi'); ylabel('\sigma'); title(sprintf('Z_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,6); pcolor(Pup,Sup,Wgam); shading flat; axis equal, axis(ax0); caxis(climsB); colorbar;
       xlabel('\phi'); ylabel('\sigma'); title(sprintf('W_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       
       %orient tall; print(gcf,'-dpng',sprintf('ZWgamma_igamma%2.2i',ii));
   end
end
