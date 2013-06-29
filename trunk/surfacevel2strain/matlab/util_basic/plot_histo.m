function [N,Nplot] = plot_histo(hdat,edges,itype)
%PLOT_HISTO plot a histogram with cyan bars and black boundaries

if nargin==2, itype=2; end

dbin = edges(2) - edges(1);
Ntotal = length(hdat);
[N,bin] = histc(hdat,edges);
switch itype
    case 1, Nplot = N; xlab = 'Count';
    case 2, Nplot = N/Ntotal; xlab = 'Fraction';
    case 3, Nplot = N/Ntotal/dbin; xlab = 'PDF';
end
bar(edges,Nplot,'histc');
xlim([min(edges) max(edges)]);
ylabel(sprintf('%s (N=%i)',xlab,Ntotal));

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 1 1],'EdgeColor','k');
