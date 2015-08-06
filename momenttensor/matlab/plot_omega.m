function [fres,N,Nplot,centers] = plot_omega(omega,nbin)
%PLOT_OMEGA plot a distribution of omega angles (compared with FMT distribution)
%
% INPUT
%   omega   angle between moment tensors, degrees
%
% see also plot_xi0.m, plot_omegadc.m
%

deg = 180/pi;
% default: 2-degree bins
if nargin==1, nbin=90; end

edges = linspace(0,180,nbin+1)'/deg;
dbin = edges(2)-edges(1);

figure; nr=4; nc=1;

% histogram of omega
subplot(nr,nc,1); hold on;
[N,Nplot,centers] = plot_histo(omega/deg,edges,3);    % samples

disp(sprintf('first bin is %.1f deg wide and contains %i samples',dbin*deg,N(1)));
disp(sprintf('last bin is %.1f deg wide and contains %i samples',dbin*deg,N(end)));
disp(sprintf('integration check: sum(Nplot)*dbin = %16.12f',sum(Nplot)*dbin));

[p,t] = omegapdf(nbin);     % homogeneous PDF

% histogram of deviations from best fit
% stairs assumes the evaluations are at the bin EDGES
%stairs(edges,omegadcpdf(edges),'r','linewidth',2);
plot(t,p,'r-','linewidth',2)
xlabel('\omega, radians');
legend('sampled PDF','homogeneous PDF','location','northwest');
title(sprintf('\\omega: min = %.3f (%.1f), max = %.3f (%.1f), n = %.0f',...
    min(omega)/deg,min(omega),max(omega)/deg,max(omega),length(omega)));

res = Nplot./p;
subplot(nr,nc,2); hold on; xmax = 180/deg; ymax = 2;
bar(centers,res,1,'c');
plot([0 pi],[1 1],'r--','linewidth',2);
xlabel('\omega, radians'); axis([0 xmax 0 ymax]);
ylabel('post / homo');

res = Nplot - p;
fres = mean(abs(res));
subplot(nr,nc,3); hold on; ymax = 0.05;
bar(centers,res,1,'c');
plot([0 pi],[0 0],'r--','linewidth',2);
xlabel('\omega, radians'); axis([0 xmax -ymax ymax]);
ylabel('post - homo');
title(sprintf('mean(abs(post - homo)) = %.4f',fres));

res = log(Nplot./p);
subplot(nr,nc,4); hold on; ymax = 0.2;
bar(centers,res,1,'c');
plot([0 pi],[0 0],'r--','linewidth',2);
xlabel('\omega, radians'); axis([0 xmax -ymax ymax]);
ylabel('ln( post / homo )');
    
%==========================================================================
