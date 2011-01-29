%
% function N = plot_histo(hdat,edges)
%
% This function uses "histc" and "bar" to plot
% a histogram with cyan bars and black boundaries.
%
% calls xxx
% called by xxx
%

function N = plot_histo(hdat,edges)

Ntotal = length(hdat);
[N,bin] = histc(hdat,edges);
bar(edges,N/Ntotal,'histc');
xlim([min(edges) max(edges)]);
ylabel([' Fraction (N = ' num2str(Ntotal) ')']);

h = findobj(gca,'Type','patch'); set(h,'FaceColor',[0 1 1],'EdgeColor','k');

%=======================================================================