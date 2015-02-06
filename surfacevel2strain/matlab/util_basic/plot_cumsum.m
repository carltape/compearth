function dcum = plot_cumsum(d,p)
%PLOT_CUMSUM plot normalized cumulative sum for a set of input points
%
% INPUT
%   d       vector or array of input values
%   p       optional: return the d value that marks p fraction of the cumulative distribution
%
% OUTPUT
%   dcum    d marking the p fraction of area under the cumulative distribution
%
% see also plot_histo.m
%

busep = true;
if nargin==1
    busep = false; if nargout==2, error('only 1 output argument allowed'); end
end

d = d(:);
dmin = min(d);
dmax = max(d);

% cumulative sum
dsort = sort(d);
dcum = cumsum(dsort);
dcumn = dcum/max(dcum);

xran = [dmin dmax];
yran = [0 1] + 0.01*[-1 1];
hold on;
plot(dsort,dcumn,'linewidth',2);
axis([xran yran]);
xlabel('x');
ylabel('p(x), normalized cumulative sum');
set(gca,'ytick',[0:0.1:1]);
%grid on;

if busep
    % from cumulative sum
    dcum = interp1(dcumn,dsort,p);    
    % percentile from matlab (this is very close to xcum)
    %xper = prctile(d,p*100);
    
    plot([0 dcum],[p p],'r--',[dcum dcum],[0 p],'r--');
    title(sprintf('p(%.2f) = %.2f',dcum,p));
end

%==========================================================================
% EXAMPLE

if 0==1
    %%
    clear, close all
    % EXAMPLE 1: simple Gaussian
    n = 1e5; a0 = 7;
    d = a0+randn(n,1);
    p = 0.85;       % percentage
    
    figure; plot_cumsum(d);
    figure; plot_cumsum(d,p);
    
    % EXAMPLE 2: non-Gaussian distribution
    b0 = 10;
    a = a0+randn(n/2,1); b = b0 + 0.5*randn(n/2,1);
    d = [a ; b];
    d(find(and(d > 6, d < 7)))=[];
    
    figure; plot_cumsum(d,p);
    
    % comparison plot with a histogram
    edges = [1:0.25:14];
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot_histo(d,edges); xlim([edges(1) edges(end)]);
    subplot(nr,nc,2); plot_cumsum(d,p); xlim([edges(1) edges(end)]);

    %title(sprintf('%.2f percentile = %.2f, cumsum = %.2f',p,px,xcum));
    %plot(px*[1 1],[0 0.1],'r--');
end

%==========================================================================
