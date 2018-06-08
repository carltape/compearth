function [N,Nplot,centers] = plot_histo(hdat,edges,itype,make_plot)
%PLOT_HISTO plot a histogram with cyan bars and black boundaries
%
% INPUT
%   hdat        input data to bin
%   edges       vector defining the edges of the bins (for hdat)
%   itype       optional: type of histogram (=1,2,3) [default = 2]
%   make_plot   optional: plot histogram [default = true]
%
% This uses Matlab's functions histc and bar.
% Carl Tape, 1/1/2008

hdat = hdat(:);
barcolor = [1 1 1]*0.8;

% default is to plot the histogram
if nargin==3, make_plot=true; end

% default is to plot the fraction of the total for each bin
if nargin==2, itype=2; make_plot=true; end

% bin width (only relevant if bins are the same width)
dbin = edges(2) - edges(1);

Ntotal = length(hdat);
[N,bin] = histc(hdat,edges);
switch itype
    case 1, Nplot = N; xlab = 'Count';
    case 2, Nplot = N/Ntotal; xlab = 'Fraction';
    case 3, Nplot = N/Ntotal/dbin; xlab = 'PDF';
        %if length(unique(edges)) > 1
        if std(diff(edges))/mean(diff(edges)) > 1e-4       % ad hoc criterion
            unique(diff(edges))
            warning('PDF is not implemented to allow bins with varying widths');
        end
    otherwise, error('itype = %i -- it must be 1,2, or 3',itype); 
end

if make_plot
    bar(edges,Nplot,'histc');
    xlim([min(edges) max(edges)]);
    ylabel(sprintf('%s (N=%i)',xlab,Ntotal));

    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',barcolor,'EdgeColor','k');

    if length(hdat) ~= sum(N)
       warning('(plot_histo.m): You may want to extend the histogram edges -->');
       disp(sprintf(' there are %i/%i input that are outside the specified range',...
           length(hdat)-sum(N),length(hdat)));
       %disp(sprintf(' the number of input (%i) does not equal the sum of bin counts (%i).',length(hdat),sum(N)));
    end
    
    % if the bin widths are different sizes, then matlab plots asterix symbols
    % at the base of the bins (WHY?), so here we delete them
    if length(unique(edges)) > 1
        sh=findall(gcf,'marker','*');
        delete(sh);
    end
end

centers = edges + dbin/2;
centers(end) = [];
% the last bin of N (and Nplot) will count any values that match EDGES(end)
% (see histc), so we cut them
N(end) = [];
Nplot(end) = [];
