function plotMT_lune(gamma,delta,gedges,dedges)
%PLOTMT_LUNE plot distributions of gamma and delta
%
% INPUT
%   gamma   lune longitudes, in degrees
%   delta   lune latitudes, in degrees
%   gedges  optional: vector of edges of bins for gamma, in degrees
%   dedges  optional: vector of edges of bins for delta, in degrees
%

deg = 180/pi;
fsize = 8;
iplothisto = 3;     % =3: plot histograms as probability density function

% everything needs to be in radians, sincer we are dealing with PDFs
gamma = gamma/deg;
delta = delta/deg;

if nargin==2
    dinc = 5;
    gedges = [-30:dinc:30]/deg;
    dedges = [-90:dinc:90]/deg;
end

n = length(gamma);
ndelta = length(dedges)-1;
ngamma = length(gedges)-1;

% analytical solution (Tape and Tape, 2015)
K = 2/3;        % int[-pi/6,pi/6]( cos(3*gamma) )
dsmooth = linspace(-30,30,100)/deg;
csmooth = (1/K) *cos(3*dsmooth);

% plot histograms of LUNE LONGITUDE for each lune latitude band
nr=6; nc=3;
for ii=1:ndelta
    % subset of lune points in this latitude band
    d1 = dedges(ii);
    d2 = dedges(ii+1);
    ikeep = find(and( delta > d1, delta <= d2 ));
    
    kk = mod(ii-1,nr*nc);
    if kk==0, figure; end
    subplot(nr,nc,kk+1); hold on;
    plot_histo(gamma(ikeep),gedges,iplothisto,true); ylim([0 2]);
    plot(dsmooth,csmooth,'r','linewidth',2);
    title(sprintf('\\delta = [%.1f, %.1f] (%.3f)',d1*deg,d2*deg,length(ikeep)/n));
    fontsize(fsize); 
end

% analytical solution (Tape and Tape, 2015)
K = 3*pi/8;         % int[0,pi]( sin(beta).^4 ) = int[-pi/2,pi/2]( cos(delta).^4 )
dsmooth = linspace(-90,90,100)/deg;
csmooth = (1/K) * cos(dsmooth).^4;

% plot histograms of LUNE LATITUDE for each lune longitude band
nr=4; nc=3;
for ii=1:ngamma
    % subset of lune points in this longitude band
    g1 = gedges(ii);
    g2 = gedges(ii+1);
    ikeep = find(and( gamma > g1, gamma <= g2 ));
    
    kk = mod(ii-1,nr*nc);
    if kk==0, figure; end
    subplot(nr,nc,kk+1); hold on;
    plot_histo(delta(ikeep),dedges,iplothisto,true); ylim([0 1]);
    plot(dsmooth,csmooth,'r','linewidth',2);
    title(sprintf('\\gamma = [%.1f, %.1f] (%.3f)',g1*deg,g2*deg,length(ikeep)/n));
    fontsize(fsize);
end
