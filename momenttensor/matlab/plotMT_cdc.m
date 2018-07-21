function plotMT_cdc(nu,alpha,phi,zeta)
%PLOTMT_CDC plot histograms of the entries of Mij

buse_greek = false;

%edges_nu = linspace(-1.001,0.501,16);    % to keep 0.50 values in the plot
edges_nu = linspace(-5,5,21);
%edges_nu = [-1:0.1:0.5];
%edges_nu = [-15:0.5:2.5];               % to show all nu values

edges_alpha = [0:5:180];
edges_phi = [-180:10:180];
edges_zeta = [0:5:90];

fsizex = 14;
itype = 1;

figure; nr=2; nc=2;
subplot(nr,nc,1); hold on; plot_histo(nu,edges_nu,itype);
if buse_greek
    xlabel('Poisson ratio \nu','fontsize',fsizex);
else
    xlabel('Poisson ratio nu','fontsize',fsizex);
end
box on;    % note sure why this is needed, but it is
%title('(a)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

subplot(nr,nc,2); plot_histo(alpha,edges_alpha,itype);
if buse_greek
    xlabel('\alpha = \angle (N, S)','fontsize',fsizex);
else
    xlabel('alpha = angle(N, S)','fontsize',fsizex);
end
set(gca,'xtick',[0:30:180]);
%title('(b)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

subplot(nr,nc,3); plot_histo(phi,edges_phi,itype);
if buse_greek
    xlabel('Azimuth on lune \phi','fontsize',fsizex);
else
    xlabel('Azimuth on lune phi','fontsize',fsizex);
end
set(gca,'xtick',[-180:60:180]);
%title('(c)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

subplot(nr,nc,4); plot_histo(zeta,edges_zeta,itype);
if buse_greek
    xlabel('Crack fraction \zeta','fontsize',fsizex);
else
    xlabel('Crack fraction zeta','fontsize',fsizex);
end
%title('(d)','Units','normalized','Position',[-0.1 1.05],'HorizontalAlignment','left','fontsize',fsizex-2);

%==========================================================================
