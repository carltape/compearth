%
% function [glon,glat,gq,nvec,axmat] = getspheregrid(ax0,qmin,qmax)
% Carl Tape, 01-Jan-2011
%
% Extract spherical-triangular gridpoints within a lat-lon box for grid
% orders qmin to qmax. This is adapted from a Fortran77 code based on the
% grids presented in 
%   Z. Wang and F. A. Dahlen, "Spherical-spline parameterization of
%   three-dimensional Earth models," Geophysical Research Letters, 1995
%
% These grids are part of the software package surfacevel2strain and
% were used in
%   Tape, Muse, Simons, Dong, Webb, "Multiscale Estimation of GPS velocity
%   fields," Geophysicsl Journal International, 2009.
%
% NOTE: The lon-lat box assumes no crossing of lon = 0 or lon = 180. We
% need to generalize this to allow for corrections.
%
% calls xxx
% called by surfacevel2strain.m
%

function [glon,glat,gq,nvec,axmat] = getspheregrid(ax0,qmin,qmax)

QMAXMAX = 14;
if qmax > QMAXMAX
    error(sprintf('qmax (%i) exceeds QMAXMAX (%i)',qmax,QMAXMAX));
end

% check if input longitudes are [0,360] or [-180,180]
if any(ax0(1:2) > 180)
    ilon360 = 1;
    axg = [0 360 -90 90];
else
    ilon360 = 0;
    axg = [-180 180 -90 90];
end
ax1 = ax0;

% UGLY FIX to split any box crossing lon=180 into two boxes -- this is no
% longer needed.
%if or(any(ax0(1:2) > 180),any(ax0(1:2) < -180))
%     ilon360 = 1;
%     ax1([1 2]) = wrapTo180(ax0([1 2]));
%     if diff(ax1([1 2])) ~= diff(ax0([1 2]))
%         DEPS = 1e-4;  % avoid repeating points on lon = +/- 180
%         ax0w = [ax1(1) 180 ax0(3:4)];
%         ax0e = [-180+DEPS ax1(2) ax0(3:4)];
%         disp('input axes box spans lon = -180 -- divide into two regions:');
%         disp(sprintf('west region: [%.6f %.6f %.1f %.1f]',ax0w));
%         disp(sprintf('east region: [%.6f %.6f %.1f %.1f]',ax0e));
%         
%         figure; nr=2; nc=1;
%         subplot(nr,nc,1);
%         plot(ax0([1 2 2 1 1]),ax0([3 3 4 4 3]),'k');
%         axis equal; axis([0 360 -90 90]); grid on;
%         subplot(nr,nc,2); hold on;
%         plot(ax0w([1 2 2 1 1]),ax0w([3 3 4 4 3]),'r');
%         plot(ax0e([1 2 2 1 1]),ax0e([3 3 4 4 3]),'b');
%         axis equal; axis([-180 180 -90 90]); grid on;
%         
%         itwobox = 1;
%     end
% else
%     ilon360 = 0;
% end

% parameters
rad = pi/180;
c72 = cos(rad*72);
base = acos(c72/(1-c72));

% determines how many gridpoints to take outside the data region
dfac = 3;

lonmin0 = ax1(1);
lonmax0 = ax1(2);
latmin0 = ax1(3);
latmax0 = ax1(4);
% thmin = (90 - latmax)*rad;
% thmax = (90 - latmin)*rad;
% phmin = lonmin*rad;
% phmax = lonmax*rad;

% if the patch does not cover the sphere, this indicates that a subset
% lon-lat patch of points is requested
Apatch = areaquad(ax1(3),ax1(1),ax1(4),ax1(2));
if abs(Apatch - 1) < 1e-4
    disp('getspheregrid.m: full sphere');
    isub = 0;
else
    disp('getspheregrid.m: lon-lat subregion of sphere');
    isub = 1;
end

%-------------------------------

gph = [];
gth = [];
gq = [];

qvec = [qmin:qmax];
nq = length(qvec);
nvec = zeros(nq,1);
axmat = zeros(nq,4);

% loop over grids
for ii = 1:nq
    q = qvec(ii);
    
   % length of triangular patch (deg), which is also the approximate length
   % scale of a basis function associated with the gridpoint
   dbase_deg = base / 2^q / rad;
    
   % If a subregion is wanted, then this is a crude way to ensure that some
   % gridpoints OUTSIDE the target region will be allowed, which is
   % important for capturing long-scalelength effects near the boundary.
   % NOTE 1: a better way would be to compute the distance from a potential
   %         gridpoint to the closest point on the lon-lat box boundary
   % NOTE 2: most of these outside gridpoints will be thresholded later
   if isub==1
        dlon = dfac*dbase_deg;
        dlat = dfac*dbase_deg;
        lonmin = lonmin0 - dlon;
        lonmax = lonmax0 + dlon;
        latmin = latmin0 - dlat;
        latmax = latmax0 + dlat;
   else
       lonmin = axg(1); lonmax = axg(2); latmin = axg(3); latmax = axg(4);
   end
    
    % NOTE: if a target region is CLOSE to lon=180, then these expanded
    %       boxes may cross lon=180; for now we do not allow the expanded
    %       axes to cross; note that most of these outside gridpoints
    %       will be thresholded later
    %if (lonmin < -180) lonmin = -180; end
    %if (lonmax >  180) lonmax =  180; end
    if (lonmin < axg(1)) lonmin = axg(1); end
    if (lonmax > axg(2)) lonmax = axg(2); end
    if (latmin < -90) latmin = -90; end
    if (latmax > 90) latmax = 90; end
    thmin = (90 - latmax)*rad;
    thmax = (90 - latmin)*rad;
    phmin = lonmin*rad;
    phmax = lonmax*rad;
    
    % display information
    disp('  ');
    disp(sprintf('GRID ORDER q = %i',q));
    disp(sprintf('  dbase = %.3e deg',dbase_deg));
    disp(sprintf('  lonmin, lonmax, latmin, latmax:'));
    disp(sprintf('  %.2f, %.2f, %.2f, %.2f',lonmin, lonmax, latmin, latmax));
    disp(sprintf('  phmin, phmax, thmin, thmax:'));
    disp(sprintf('  %.4f, %.4f, %.4f, %.4f',phmin, phmax, thmin, thmax));
    
    ngrid = 10*(2^(2*q)) + 2;
    if isub==1
        fA = areaquad(latmin,lonmin,latmax,lonmax);
        disp(sprintf('Patch occupies this fraction of the sphere: %.3e',fA));
        disp(sprintf('This corresponds to a square patch with side length %.3f deg',sqrt(4*pi*fA) / rad));
        ngrid_sub = round(ngrid * fA);
    else
        ngrid_sub = ngrid;
    end
    
    disp(sprintf('getting the gridpoints for q = %i',q));
    disp(sprintf('  number of total possible gridpoints : %i',ngrid));
    disp(sprintf('  maximum number of subset gridpoints : %i',ngrid_sub));
    
    % KEY: get gridpoints within the lon-lat box
    nf = 2^q;
    [gph0,gth0] = mytessa(q,nf,thmin,thmax,phmin,phmax);
    nknt = length(gph0);
    
    disp(sprintf('  actual number of subset gridpoints : %i',nknt));
    
    % compile full set of gridpoints
    gph = [gph ; gph0];
    gth = [gth ; gth0];
    gq = [gq ; q*ones(nknt,1) ];
    nvec(ii) = nknt;
    axmat(ii,:) = [lonmin lonmax latmin latmax];
end

% convert to longitude and latitude
glon = gph/rad;
glat = 90 - gth/rad;

%=========================================================================
% LOCAL FUNCTIONS

function [gph,gth] = mytessa(q,nf,thmin,thmax,phmin,phmax)

% parameters
rad = pi/180;
c72 = cos(rad*72);
base = acos(c72/(1-c72));
dbase = base/nf;

[gph01,gth01] = triangle01(q,nf);
[gph67,gth67] = triangle67(q,nf);
[gph16,gth16] = triangle16(q,nf);

%figure; plot(gph01,gth01,'.');
%figure; plot(gph67,gth67,'.');
%figure; plot(gph16,gth16,'.');

gph = [];
gth = [];

if any([phmin phmax] > pi), ilon360 = 1; else ilon360 = 0; end

disp(sprintf('q = %i, nf = %i, dbase = %.3e deg',q,nf,dbase/rad));

i = 1;
j = 1;
nknt = 0;   % counter for number of knots

% North Pole
arco = 0;
arlo = 0;

% for the North Pole, you just care about the proximity in latitude
if and(arco >= thmin, arco <= thmax)
    disp(sprintf('q = %i: first point is North Pole',q))
    gph = [gph ; arlo];
    gth = [gth ; arco];
    nknt = nknt+1;
end

for i = 2:nf+1
    j1 = 1;
    j2 = i;
    knt = 0;
    for n = 1:5
        if (n > 1) j1 = 2; end
        for j = j1:i
            if and(n==5, j==i) break; end
            knt = knt+1;

            arco = gth01(i,j);
            arlo = gph01(i,j) + (n-1)*rad*72;
            %if (arlo > pi) arlo = arlo-2*pi; end
            if ilon360==1, arlo = wrapTo2Pi(arlo); else arlo = wrapToPi(arlo); end
            if and( and( arlo >= phmin , arlo <= phmax), and(arco >= thmin, arco <= thmax))
                gph = [gph ; arlo];
                gth = [gth ; arco];
                nknt = nknt+1;
            end
        end
    end
end

i1 = nf+1;
for i = 2:nf+1
    knt = 0;
    for n = 1:5
        j1 = 1;
        j2 = nf+1;
        if (n > 1) j1=2; end
        for j = j1:j2
            if and(n==5, j==nf+1) break; end
            knt = knt+1;

            arco = gth67(i,j);     
            arlo = gph67(i,j) + (n-1)*rad*72;
            %if (arlo > pi) arlo = arlo-2*pi; end
            if ilon360==1, arlo = wrapTo2Pi(arlo); else arlo = wrapToPi(arlo); end
            if and( and( arlo >= phmin, arlo <= phmax), and(arco >= thmin, arco <= thmax))
                gph = [gph ; arlo];
                gth = [gth ; arco];
                nknt = nknt+1;
            end
        end
    end
end

for i = 2:nf+1
    knt = 0;
    for n = 1:5
        j1 = 1;
        if (n > 1) j1=2; end
        j2 = nf+2-i;
        for j = j1:j2
            if and(n==5, j==j2) break; end
            knt = knt+1;

            arco = gth16(i,j);
            arlo = gph16(i,j) + (n-1)*rad*72;
            %if (arlo > pi) arlo = arlo-2*pi; end
            if ilon360==1, arlo = wrapTo2Pi(arlo); else arlo = wrapToPi(arlo); end
            if and( and( arlo >= phmin, arlo <= phmax), and(arco >= thmin, arco <= thmax))
                gph = [gph ; arlo];
                gth = [gth ; arco];
                nknt = nknt+1;
            end
        end
    end
end

%------------------------------------------------------------------------

function [tph,tth] = triangle01(q,nf)

% parameters
rad = pi/180;
c72 = cos(rad*72);
base = acos(c72/(1-c72));

tph = zeros(1+2^q,1+2^q);
tth = tph;

dbase = base/nf;
tph(1,1) = 0;
tth(1,1) = 0;

for i = 2:nf+1
    tph(i,1) = 0;
    y = dbase*(i-1);
    tth(i,1) = y;
    x = (cos(y))^2 + c72*(sin(y))^2;
    xx = acos(x);
    xs = xx/(i-1);
    sang = sin(y)*sin(72*rad)/sin(xx);
    ang = asin(sang);
    for j = 2:i
        xm = xs*(j-1);
        cl = cos(y)*cos(xm)+sin(y)*sin(xm)*cos(ang);
        acl = acos(cl);
        tth(i,j) = acl;
        slong = sang*sin(xm)/sin(acl);
        tph(i,j) = asin(slong);
    end
end

%------------------------------------------------------------------------

function [tph,tth] = triangle67(q,nf)

% parameters
rad = pi/180;
c72 = cos(rad*72);
base = acos(c72/(1-c72));

s72 = sin(rad*72);
c144 = cos(rad*144);
s144 = sin(rad*144);
dbase = base/nf;

tph = zeros(1+2^q,1+2^q);
tth = tph;
  
tth(1,1) = base;
tph(1,1) = 0;

 for i = 2:nf+1
     tph(i,1) = 0;
     x = dbase*(i-1);
     ccl = cos(base)*cos(x)+sin(base)*sin(x)*c72;
     alat = acos(ccl);
     tth(1,i) = alat;
     tph(1,i) = asin(s72*sin(x)/sin(alat));
end
  
  for i = 2:nf+1
     x = dbase*(i-1);
     ccl = cos(base)*cos(x)+sin(base)*sin(x)*c144;
     accl = acos(ccl);
     tth(i,1) = accl;
     sgamma = sin(x)*s144/sin(accl);
     gamma = asin(sgamma);
     tph(i,1) = gamma;
     if (i < nf+1)
        top = 72*rad - 2*gamma;
        y = (cos(accl))^2 + cos(top)*(sin(accl))^2;
        ay = acos(y);
        cang = (ccl-ccl*y)/(sin(ay)*sin(accl));
        sang = sin(top)*sin(accl)/sin(ay);
        ang = acos(cang);
        day = ay/(nf+1-i);
        for j = 2:nf+2-i
           yy = day*(j-1);
           cclat = cos(tth(i,1))*cos(yy) + sin(tth(i,1))*sin(yy)*cos(ang);
           aclat = acos(cclat);
           tth(i,j) = aclat;
           sfi = asin(sin(yy)*sang/sin(aclat));
           tph(i,j) = tph(i,1) + sfi;
        end
     end
  end
  
  for i = 1:nf
     ii = nf+2-i;
     for j = 1:nf+1-i
       jj = nf+2-j;
       tth(ii,jj) = pi - tth(i,j);
       tph(ii,jj) = 108*rad - tph(i,j);
     end
  end

%------------------------------------------------------------------------

function [tph1,tth1] = triangle16(q,nf)

% parameters
rad = pi/180;
c72 = cos(rad*72);
base = acos(c72/(1-c72));

tph = zeros(1+2^q,1+2^q);
tth = tph;

xl36 = 36*rad;
dbase = base/nf;
tth(1,1) = 0;
tph(1,1) = 0;

for i = 2:nf+1
    tph(i,1) = 0;
    y = dbase*(i-1);
    tth(i,1) = y;
    x = (cos(y))^2 + c72*(sin(y))^2;
    xx = acos(x);
    xs = xx/(i-1);
    sang = sin(y)*sin(72*rad)/sin(xx);
    ang = asin(sang);
    for j = 2:i
      xm = xs*(j-1);
      cl = cos(y)*cos(xm) + sin(y)*sin(xm)*cos(ang);
      acl = acos(cl);
      tth(i,j) = acl;
      slong = sang*sin(xm)/sin(acl);
      tph(i,j) = asin(slong);
    end
end

for i = 1:nf+1
    ii = nf+2-i;
    for j = 1:i
      tth1(ii,j) = pi - tth(i,j);
      tph1(ii,j) = xl36 + tph(i,j);
    end
end

%=========================================================================
