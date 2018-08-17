function plotMT_lam(lam)
%PLOTMT_LAM plot histograms of three eigenvalues

PLOT_UNIFORM_CURVES = true;
PLOT_AS_PDFS = true;

if PLOT_AS_PDFS, itype = 3; else itype = 2; end

[n1,n2] = size(lam);
if n1~=3, error('lam must be 3 x n'); end

lam1 = lam(1,:);
lam2 = lam(2,:);
lam3 = lam(3,:);
mnorm = sqrt( lam1.^2 + lam2.^2 + lam3.^2 );

if PLOT_UNIFORM_CURVES
    mfac = mean(mnorm);
    if std(mnorm) > 0.01*mfac
        warning('input moment tensors do NOT have the same norm');
        disp('--> setting PLOT_UNIFORM_CURVES = false');
        PLOT_UNIFORM_CURVES = false;
    else
        disp(sprintf('input moment tensors have the same norm (%.3e)',mfac));
    end
end

nbin = 36;
nedge = nbin + 1;
lmax = max(abs(lam(:)));
edges = linspace(-lmax,lmax,nedge);

figure; nr=3; nc=1;
subplot(nr,nc,1); hold on;
plot_histo(lam1,edges,itype); 
if PLOT_UNIFORM_CURVES
    npt = 500;
    xplot = linspace(-lmax,lmax,npt);
    uplot = plam(xplot,1,mfac);
    plot(xplot,uplot,'r','linewidth',2);
end
xlabel('fist eigenvalue');
xlim(1.05*lmax*[-1 1]);

subplot(nr,nc,2); hold on;
plot_histo(lam2,edges,itype);
if PLOT_UNIFORM_CURVES
    uplot = plam(xplot,2,mfac);
    plot(xplot,uplot,'r','linewidth',2);
end
xlabel('second eigenvalue');
xlim(1.05*lmax*[-1 1]);

subplot(nr,nc,3); hold on;
plot_histo(lam3,edges,itype);
if PLOT_UNIFORM_CURVES
    uplot = plam(xplot,3,mfac);
    plot(xplot,uplot,'r','linewidth',2);
end
xlabel('third eigenvalue');
xlim(1.05*lmax*[-1 1]);

%==========================================================================

function p = plam(x,ilam,mfac)

p = zeros(size(x));

if ilam==1, x = -x; end

z1 = -1 * mfac;
z2 = -1/sqrt(2) * mfac;
z3 = -1/sqrt(3) * mfac;
z4 =  1/sqrt(3) * mfac;

i1 = find(x < z1);
i2 = find(and(x >= z1, x < z2));
i3 = find(and(x >= z2, x < z3));
i4 = find(and(x >= z3, x <= z4));
i5 = find(x > z4);

% if ~isempty(i1)
%     ilam, x(i1), warning('input lambdas outside left endpoint');
% end
% if ~isempty(i5)
%     ilam, x(i5), warning('input lambdas outside right endpoint');
% end

if ilam==2
    % lambda2
    inds = find(and(x >= z2, x <= -z2));
    xx = x(inds)/mfac;
    p(inds) = 4/(3*pi) / mfac * (2 - 4*xx.^2).^(3/2);
        
else
    % equation is for lamda3(x)
    x2 = x(i2)/mfac;
    x3 = x(i3)/mfac;
    x4 = x(i4)/mfac;
    p(i2) = 8/(3*pi) / mfac * (7*x2.^2 - 1) .* (1 - x2.^2).^(1/2);
    p(i3) = 8/(3*pi) / mfac * ( (7*x3.^2 - 1) .* (1 - x3.^2).^(1/2) ...
                            + sqrt(2)*(1 - 2*x3.^2).^(3/2) );
    p(i4) = sqrt(8)/(3*pi) / mfac * (sqrt(2)* (7*x4.^2 - 1) .* (1 - x4.^2).^(1/2) ...
                + 2*(1 - 2*x4.^2).^(3/2) - x4.*(x4.^2 + 3) );
        
end

%==========================================================================
