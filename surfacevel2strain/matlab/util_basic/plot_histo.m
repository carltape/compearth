function N = plot_histo(hdat,edges)
%PLOT_HISTO plot a histogram with cyan bars and black boundaries

Ntotal = length(hdat);
[N,bin] = histc(hdat,edges);
bar(edges,N/Ntotal,'histc');
xlim([min(edges) max(edges)]);
ylabel(sprintf('p(x) (N=%i)',Ntotal));
%ylabel(['Fraction (N = ' num2str(Ntotal) ')']);

h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 1 1],'EdgeColor','k');