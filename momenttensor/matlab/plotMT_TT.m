function plotMT_TT(gamma,delta,M0,kappa,theta,sigma)
%PLOTMT_TT plot histograms of the parameters of TapeTape2012
%
% INPUT
%   gamma       angle from DC meridian to MT point (-30 <= gamma <= 30)
%   delta       angle from deviatoric plane to MT point (-90 <= delta <= 90)
%   M0          seismic moment, N-m
%   kappa       strike angle, degrees: [0,360]
%   theta       dip angle, degrees: [0,90]
%   sigma       slip (or rake) angle, degrees: [-90,90]
% 
% See example below.
%

% USER PARAMETERS
bfixbinwidthdip = true;    % fix bin width or scale to make uniform distributions flat
MAGPERT = 0.5;      % the percent perturbation is approximately 100*MAGPERT
nbin = 18;          % number of bins in each histogram

deg = 180/pi;
iplothisto = 3;

% it is safer to handle everything in radians, sincer we are dealing with PDFs
% here we will convert only the lune coordinates and dip
gamma = gamma/deg;
delta = delta/deg;
theta = theta/deg;

% we plot the magnitudes on a log scale w.r.t. median value
M0norm = 10^median(log10(M0));
stM0 = sprintf('ln( M0 / %.2e )',M0norm);
pM0 = log(M0/M0norm);

nedge = nbin + 1;
hflat = 1/nbin;

gammabin = linspace(-30,30,nedge)/deg;      % radians
deltabin = linspace(-90,90,nedge)/deg;      % radians
magbin   = MAGPERT * linspace(-1,1,nedge);
kappabin = linspace(0,360,nedge);
thetabin = acos(linspace(1,0,nedge));       % radians
sigmabin = linspace(-90,90,nedge);
if bfixbinwidthdip
    thetabin = linspace(0,90,nedge)/deg;    % radians
end

figure; nr=3; nc=2;

% gamma
gsmooth = linspace(-30,30,100)/deg;
fsmooth = cos(3*gsmooth) / (2/3);       % Tape and Tape (2015)
subplot(nr,nc,1); hold on;
plot_histo(gamma,gammabin,iplothisto);
plot(gsmooth,fsmooth,'r','linewidth',2);
% WARNING: angles are in radians, and the pdf is based on calculations in
% radians, but here we label the bin values as degrees
xlabs = [-30:6:30];
set(gca,'xtick',xlabs/deg,'xticklabel',numvec2cell(xlabs));
xlabel('lune longitude, deg');
%title(sprintf('median = %.2f deg, mean = %.2f deg',median(gamma)*deg,mean(gamma)*deg));
%plot([-30 30], hflat*[1 1],'r','linewidth',2);

% delta
dsmooth = linspace(-90,90,100)/deg;
fsmooth = cos(dsmooth).^4 / (3*pi/8);   % Tape and Tape (2015)
subplot(nr,nc,2); hold on;
plot_histo(delta,deltabin,iplothisto);
plot(dsmooth,fsmooth,'r','linewidth',2);
% WARNING: angles are in radians, and the pdf is based on calculations in
% radians, but here we label the bin values as degrees
xlabs = [-90:30:90];
set(gca,'xtick',xlabs/deg,'xticklabel',numvec2cell(xlabs));
xlabel('lune latitude, deg');
%plot([-90 90], hflat*[1 1],'r','linewidth',2);

% magnitude
subplot(nr,nc,3);
plot_histo(pM0,magbin); xlabel(stM0); 

% strike angle
subplot(nr,nc,4); hold on;
plot_histo(kappa,kappabin); xlabel('strike, deg'); set(gca,'xtick',0:60:360);
plot([0 360], hflat*[1 1],'r','linewidth',2);

% dip angle
tsmooth = linspace(0,90,100)/deg;
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

%==========================================================================
