function [p,t] = kaganpdf(x,bfigure)
%KAGANPDF probability density function for random rotations
%
% INPUT
%   x   A: an integer number of points for discretization OR
%       B: an input vector of points with for which to evaluate the PDF
%           INPUT IS IN RADIANS, NOT DEGREES
%
% OUTPUT
%   p   PDF discretized based on input choice
%   t   discretization of xi0 minimum rotation angles (IN RADIANS)
%
% See Tape and Tape (2012), "Angle between principal axis triples"
%
% EXAMPLES (see longer example below):
%   t = linspace(0,2*pi/3,100); p = kaganpdf(t,true);
%   [p,t] = kaganpdf(10,true);
%

if nargin==1, bfigure=false; end

% key xi0 points on curve
t0 = 0;
t1 = pi/2;
t2 = acos(-1/3);
t3 = 2*pi/3;

if length(x)==1
   n = round(x); 
   dt = t3/n;
   t = [dt/2 : dt : t3-dt/2];   % centers of bins
   tx1 = t - dt/2;              % left edge of bin
   tx2 = t + dt/2;              % right edge of bin
else
   t = x;
   n = length(t);
   dt = t(2)-t(1);
   tx1 = [];
   tx2 = [];
end

% evaluate the PDF at the bin centers
p = kaganfun(t);

if length(x)==1
    % numerically integrate the PDF for each bin
    pint = zeros(n,1);
    xtol = 1e-12;
    for ii=1:n
        pint(ii) = integral(@(x) kaganfun(x), tx1(ii), tx2(ii),'AbsTol',xtol,'RelTol',xtol);
    end
    
    % integration check
    disp(sprintf('integration check using %i bars: %.14f',n,sum(p)*dt));
    disp(sprintf('integration check using matlab: %.14f',sum(pint)));

    p = pint/dt;
end

% % cumulative density function
% cint = cumsum(pint);    % pint will be very accurate
% %cint = zeros(n,1);
% %for ii=1:n
% %    cint(ii) = integral(fun, t0, tx2(ii),'AbsTol',xtol,'RelTol',xtol);
% %end

if bfigure
   tsmooth = linspace(t0,t3,400);
   psmooth = kaganfun(tsmooth);
   pt1 = 4/pi;
   pt2 = 4/pi*(-8/3 + sqrt(8));
   
   figure; nr=2; nc=1; 
   
   subplot(nr,nc,1); hold on;
   plot(tsmooth,psmooth,'b');
   plot(t,p,'r.--');
   plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   plot([t0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([0 t3 0 1.5]);
   xlabel('t, radians');
   ylabel('p_{\xi_0}(t)');
   title(sprintf('PDF for random orientations (n = %i)',n));
   
   subplot(nr,nc,2); hold on;
   bar(t,p,1,'c');
   plot(t,p,'r.');
   plot(tsmooth,psmooth,'b');
   plot([t1 t1],[0 pt1],'k--',[t2 t2],[0 pt2],'k--');
   plot([0 t1 t2 t3],[0 pt1 pt2 0],'bo','markersize',14,'markerfacecolor','r','linewidth',1);
   axis equal; axis([t0 t3 0 1.5]);
   xlabel('t, radians');
   ylabel('p_{\xi_0}(t)');
   title(sprintf('PDF for random orientations (n = %i)',n));
end

%--------------------------------------------------------------------------

function p = kaganfun(t)

t0 = 0;
t1 = pi/2;
t2 = acos(-1/3);
t3 = 2*pi/3;

i1 = find(and(t >= t0, t <= t1));
i2 = find(and(t > t1, t <= t2));
i3 = find(and(t > t2, t <= t3));

% Tape and Tape, 2012, Eq 59 amd B4
% first interval
p(i1) = 4/pi*(1-cos(t(i1)));
% second interval
p(i2) = 4/pi*(-2 + 2*cos(t(i2)) + 3*sin(t(i2)));
% third interval
th0 = acos( cos(t(i3)/2) ./ sqrt(-cos(t(i3))) );
fA = th0  - tan(t(i3)/2) .* asin(sin(th0)/sqrt(2));
fB = pi/4 - tan(t(i3)/2) .* asin(sin(pi/4)/sqrt(2));
p(i3) = 48*sin(t(i3))/pi^2 .* (fB - fA);

%==========================================================================
% EXAMPLE

if 0==1
    % plot the PDF
    t = linspace(0,2*pi/3,121); p = kaganpdf(t,true);
    [p,t] = kaganpdf(10,true);
    
    % see also plot_xi0.m
    n = 1e6;    % try 1e7
    U = randomU(n);         % random rotation matrices
    xi0 = U2xi0(U,1,0);
    nbin = 120;
    plot_xi0(xi0,nbin);
end
    
%==========================================================================
