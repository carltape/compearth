function [S,T,a,b,c] = STgamma(mu,nu,gamma)
% 
%
% INPUT
%   mu      in radians
%   nu      in radians
%   gamma   in radians
%

[a,b,c] = abcgamma(mu,nu,gamma);

S = -b ./ (2*a);
T = -(b.^2 - 4*a.*c)./(4*a);

if 0==1
   %%
   close all
   nmu = 50;
   nnu = 2*nmu;
   deg = 180/pi;
   muvec = linspace(-pi,pi,nmu);
   nuvec = linspace(-pi,pi,nnu);
   [Mu,Nu] = meshgrid(muvec,nuvec);
   
   ax0 = [-180 180 -180 180]; clims = [0 2]; climsB = [-1 1];
   dticks = [-180:90:180];
   
   gvec = [-25:5:25]/deg;
   for ii = 1:length(gvec);
       gamma = gvec(ii);
       [Sgam,Tgam,Agam,Bgam,Cgam] = STgamma(Mu,Nu,gamma);
       Mup = Mu*deg;
       Nup = Nu*deg;

       figure; nr=3; nc=2; ax0 = [-180 180 -180 180]; clims = [0 2];
       subplot(nr,nc,1); pcolor(Mup,Nup,Agam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\mu'); ylabel('\nu'); title(sprintf('a_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,2); pcolor(Mup,Nup,Bgam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\mu'); ylabel('\nu'); title(sprintf('b_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,3); pcolor(Mup,Nup,Cgam); shading flat; axis equal, axis(ax0); caxis(clims); colorbar;
       xlabel('\mu'); ylabel('\nu'); title(sprintf('c_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       
       subplot(nr,nc,5); pcolor(Mup,Nup,Sgam); shading flat; axis equal, axis(ax0); caxis(climsB); colorbar;
       xlabel('\mu'); ylabel('\nu'); title(sprintf('S_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
       subplot(nr,nc,6); pcolor(Mup,Nup,Tgam); shading flat; axis equal, axis(ax0); caxis(climsB); colorbar;
       xlabel('\mu'); ylabel('\nu'); title(sprintf('T_\\gamma for \\gamma = %.1f',gamma));
       set(gca,'xtick',dticks,'ytick',dticks);
   end
end
