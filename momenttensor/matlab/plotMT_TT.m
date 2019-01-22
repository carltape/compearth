function plotMT_TT(gamma,delta,M0,kappa,theta,sigma,bvw)
%PLOTMT_TT plot histograms of the parameters of TapeTape2012
%
% INPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%   M0          seismic moment, N-m
%   kappa       strike angle, degrees: [0,360]
%   theta       dip angle, degrees: [0,90]
%   sigma       slip (or rake) angle, degrees: [-90,90]
%   bvw         OPTIONAL: =true if using vw coordinates instead of lune coordinates
%
% bvw = false (DEFAULT)
%   gamma is lune longitude
%   delta is lune latitude
%
% bvw = true
%   gamma is v (TapeTape2015)
%   delta is w = 3*pi/8 - u (TapeTape2015)
%
% Note that h = cos(dip) is the parameter needed for uniformity of
% orientations. Instead of plotting the histograms for dip, we could plot
% histograms for h.
% 
% See example in run_uniformMT.m
%

if nargin==6, bvw = false; end

% USER PARAMETERS
bfixbinwidthdip = true;     % fix bin width or scale to make uniform distributions flat
MAGPERT = 0.5;              % the percent perturbation is approximately 100*MAGPERT
nbin = 18;                  % number of bins in each histogram
ns = 100;                   % number of points on smoothed curves

deg = 180/pi;
iplothisto = 3;

theta = theta/deg;          % convert dip to radians
if ~bvw
    % use lune longitude and lune latitude of Tape and Tape (2012)
    g1 = -30/deg; g2 = 30/deg;
    d1 = -90/deg; d2 = 90/deg;
    glab = 'lune longitude, deg';
    dlab = 'lune latitude, deg';
    stag = 'lune';
    
    % it is safer to handle everything in radians, since we are dealing
    % with PDFs; here we will convert only the lune coordinates
    gamma = gamma/deg;
    delta = delta/deg;
    
    % see Tape and Tape (2015)
    % lune longitude
    gsmooth  = linspace(g1,g2,ns);
    xgsmooth = cos(3*gsmooth) / (2/3);
    % lune latitude
    dsmooth  = linspace(d1,d2,ns);
    xdsmooth = cos(dsmooth).^4 / (3*pi/8);
else
    % use v-w coordinates of Tape and Tape (2015)
    g1 = -1/3; g2 = 1/3;
    d1 = -3*pi/8; d2 = 3*pi/8;
    glab = 'v (lune longitude)';
    dlab = 'w (lune latitude)';
    stag = 'vw';
    
    gsmooth  = linspace(g1,g2,ns);
    xgsmooth = ones(ns,1) / (2/3);
    dsmooth  = linspace(d1,d2,ns);
    xdsmooth = ones(ns,1) / (3*pi/4);
end

disp(sprintf('plotMT_TT.m: %i moment tensors, %s coordinates',length(gamma),stag));
if bfixbinwidthdip
    disp('bfixbinwidthdip = true: dip bins have fixed widths');
else
    disp('bfixbinwidthdip = false: dip bin widths vary (h = cos(theta))');
end

% we plot the magnitudes on a log scale w.r.t. median value
M0norm = 10^median(log10(M0));
stM0 = sprintf('ln( M0 / %.2e )',M0norm);
pM0 = log(M0/M0norm);

nedge = nbin + 1;
hflat = 1/nbin;

gammabin = linspace(g1,g2,nedge);           % radians
deltabin = linspace(d1,d2,nedge);           % radians
magbin   = MAGPERT * linspace(-1,1,nedge);
kappabin = linspace(0,360,nedge);           % degrees
sigmabin = linspace(-90,90,nedge);          % degrees
if bfixbinwidthdip
    thetabin = linspace(0,90,nedge)/deg;    % radians
else
    thetabin = acos(linspace(1,0,nedge));   % radians
end

figure; nr=3; nc=2;

% gamma
subplot(nr,nc,1); hold on;
plot_histo(gamma,gammabin,iplothisto);
plot(gsmooth,xgsmooth,'r','linewidth',2);
if ~bvw
    % WARNING: angles are in radians, and the pdf is based on calculations in
    % radians, but here we label the bin values as degrees
    xlabs = [-30:6:30];
    set(gca,'xtick',xlabs/deg,'xticklabel',numvec2cell(xlabs));
end
xlabel(glab);
%title(sprintf('median = %.2f deg, mean = %.2f deg',median(gamma)*deg,mean(gamma)*deg));
%plot([-30 30], hflat*[1 1],'r','linewidth',2);

% delta
subplot(nr,nc,2); hold on;
plot_histo(delta,deltabin,iplothisto);
plot(dsmooth,xdsmooth,'r','linewidth',2);
if ~bvw
    % WARNING: angles are in radians, and the pdf is based on calculations in
    % radians, but here we label the bin values as degrees
    xlabs = [-90:30:90];
    set(gca,'xtick',xlabs/deg,'xticklabel',numvec2cell(xlabs));
end
xlabel(dlab);
%plot([-90 90], hflat*[1 1],'r','linewidth',2);

% magnitude
subplot(nr,nc,3);
plot_histo(pM0,magbin); xlabel(stM0); 

% strike angle
subplot(nr,nc,4); hold on;
plot_histo(kappa,kappabin); xlabel('strike, deg'); set(gca,'xtick',0:60:360);
plot([0 360], hflat*[1 1],'r','linewidth',2);

% dip angle
tsmooth = linspace(0,90,ns)/deg;
fsmooth = sin(tsmooth);
subplot(nr,nc,5); hold on;
if bfixbinwidthdip
    plot_histo(theta,thetabin,iplothisto);
    plot(tsmooth,fsmooth,'r','linewidth',2);
else
    plot_histo(theta,thetabin,2);
    plot([0 90]/deg, hflat*[1 1],'r','linewidth',2);
end
% WARNING: angles are in radians, and the pdf is based on calculations in
% radians, but here we label the bin values as degrees
xlabs = [0:30:90];
set(gca,'xtick',xlabs/deg,'xticklabel',numvec2cell(xlabs));
xlabel('dip, deg');

% slip angle (rake)
subplot(nr,nc,6); hold on;
plot_histo(sigma,sigmabin); xlabel('rake, deg'); set(gca,'xtick',-180:30:180);
plot([-90 90], hflat*[1 1],'r','linewidth',2);

%--------------------------------------------------------------------------

if bvw
    % the five histograms of v, w, kappa, h, sigma should all be flat for a
    % uniform distribution of moment tensors
    figure; nr=3; nc=2;

    % v
    subplot(nr,nc,1); hold on;
    plot_histo(gamma,gammabin,iplothisto);
    plot(gsmooth,xgsmooth,'r','linewidth',2);
    xlabel('v');

    % w
    subplot(nr,nc,2); hold on;
    plot_histo(delta,deltabin,iplothisto);
    plot(dsmooth,xdsmooth,'r','linewidth',2);
    xlabel('w');
    %plot([-90 90], hflat*[1 1],'r','linewidth',2);

    % magnitude
    subplot(nr,nc,3);
    plot_histo(pM0,magbin); xlabel(stM0); 

    % strike angle
    subplot(nr,nc,4); hold on;
    plot_histo(kappa/deg,kappabin/deg); xlabel('strike, radians');
    plot([0 360]/deg, hflat*[1 1],'r','linewidth',2);

    % cosine of dip angle
    thetabin = linspace(0,1,nedge);
    subplot(nr,nc,5); hold on;
    plot_histo(cos(theta),thetabin,2); xlabel('h = cos(dip)');
    plot([0 90]/deg, hflat*[1 1],'r','linewidth',2);

    % slip angle (rake)
    subplot(nr,nc,6); hold on;
    plot_histo(sigma/deg,sigmabin/deg); xlabel('rake, radians');
    plot([-90 90]/deg, hflat*[1 1],'r','linewidth',2);
end

%==========================================================================
