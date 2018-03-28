function plotMT_Mij(M,edges,Mlab)
%PLOTMT_MIJ plot histograms of the entries of Mij
%
% M      6 x n moment tensors, M = [M11 M22 M33 M12 M13 M23]
% edges  Ne-dimensional vector for bin edges (or 6 x Ne sets of edges)
% Mlab   optional: labels for each subplot
%

PLOT_UNIFORM_CURVES = true;
PLOT_AS_PDFS = true;

if PLOT_AS_PDFS, itype = 3; else itype = 2; end

% make sure M is 6 x n
[M,n] = Mdim(M);

if nargin~=3
    Mlab = {'M11','M22','M33','M12','M13','M23'};
end

% edges can be specified for all six entries
% (for now we require the same number of bins)
USE_DIFFERENT_EDGES = false;
if ~any(size(edges)==1)
    USE_DIFFERENT_EDGES = true;
    [n1,n2] = size(edges);
    if n2==6, edges=edges'; end     % default 6 x Ne
end

mnorm = norm_MT(M);
if PLOT_UNIFORM_CURVES
    mfac = mean(mnorm);
    if std(mnorm) > 0.01*mfac
        warning('input moment tensors do NOT have the same norm');
        disp('--> setting PLOT_UNIFORM_CURVES = false');
        PLOT_UNIFORM_CURVES = false;
    else
        disp('input moment tensors have the same norm');
    end
end

figure; nr=3; nc=2;
pvec = [1 3 5 2 4 6];
for ii=1:6
   subplot(nr,nc,pvec(ii)); hold on;
   if USE_DIFFERENT_EDGES, xedges = edges(ii,:); else xedges = edges; end
   plot_histo(M(ii,:),xedges,itype);
   if PLOT_UNIFORM_CURVES
       dbin = xedges(2) - xedges(1);
       %xplot = linspace(min(xedges),max(xedges),100);
       if ii <= 3   % M11, M22, M33
           xplot = linspace(-mfac,mfac,100);
           uplot = 8/(3*pi) / mfac * (1 - xplot.^2 / mfac^2).^(3/2);
       else         % M12, M13, M23
           xplot = linspace(-mfac/sqrt(2),mfac/sqrt(2),100);
           uplot = 8*sqrt(2)/(3*pi) / mfac * (1 -2*xplot.^2 / mfac^2).^(3/2);
       end
       if ~PLOT_AS_PDFS, uplot = uplot*dbin; end
       plot(xplot,uplot,'r','linewidth',2);
   end
   xlabel(Mlab{ii});
end

Mdiag = M(1:3,:);
Moff  = M(4:6,:);

disp(sprintf('plotMT_Mij.m (%i moment tensors):',n));
disp(sprintf('Mnorm min/mean/max  : %.3e/%.3e/%.3e',min(mnorm),mean(mnorm),max(mnorm)));
disp(sprintf('M0    min/mean/max  : %.3e/%.3e/%.3e',min(mnorm)/sqrt(2),mean(mnorm)/sqrt(2),max(mnorm)/sqrt(2)));
disp(sprintf('Diagonal entries    : min(Mij) = %e, max(Mij) = %e',min(Mdiag(:)),max(Mdiag(:))));
disp(sprintf('Off-diagonal entries: min(Mij) = %e, max(Mij) = %e',min(Moff(:)),max(Moff(:))));

if min(M(:)) < edges(1),
    disp(sprintf('plotMT_Mij.m WARNING: min(M) (%e) < edges(1) (%e)',min(M(:)),edges(1)));
end
if max(M(:)) > edges(end),
    disp(sprintf('plotMT_Mij.m WARNING: max(M) (%e) < edges(2) (%e)',max(M(:)),edges(end)));
end

%==========================================================================
