function plotMT_lam(lam,edges)
%PLOTMT_LAM plot histograms of three eigenvalues

if nargin==1
    nbin = 18;
    nedge = nbin + 1;
    lmin = min(lam(:));
    lmax = max(lam(:));
    edges = linspace(lmin,lmax,nedge);
end

[n1,n2] = size(lam);
if n1~=3, error('lam must be 3 x n'); end

lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);

figure; nr=3; nc=1;
subplot(nr,nc,1); hold on; plot_histo(lam1,edges); 
xlabel('fist eigenvalue');

subplot(nr,nc,2); hold on; plot_histo(lam2,edges);
xlabel('second eigenvalue');

subplot(nr,nc,3); hold on; plot_histo(lam3,edges);
xlabel('third eigenvalue');
