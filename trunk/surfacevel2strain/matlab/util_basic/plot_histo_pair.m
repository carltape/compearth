function plot_histo_pair(hdat1,hdat2,st1,st2,edges)
%PLOT_HISTO_PAIR plots two dets of data as a superimposed histogram

hold on;

Ntotal = length(hdat1);
[N,bin] = histc(hdat1,edges);
bar(edges,N/Ntotal,'histc');
h = findobj(gca,'Type','patch'); set(h,'FaceColor',[0 1 1],'EdgeColor','k');

Ntotal = length(hdat2);
[N,bin] = histc(hdat2,edges);
j = bar(edges,N/Ntotal,'histc');
set(j,'FaceColor','none','EdgeColor','r','linewidth',3);
xlim([min(edges) max(edges)]);
ylabel('Fraction of total');
legend({sprintf('%s (N = %i)',st1,length(hdat1)),...
        sprintf('%s (N = %i)',st2,length(hdat2))},...
        'interpreter','none');
% note: need to add xlabel
