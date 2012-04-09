%
% function spline_tot = sphereinterp_grid(dlon,dlat,ax0,qparm)
% Carl Tape and Pablo Muse, 04-Jan-2011
%
% This function obtains the centerpoints and scale indices for a set of
% spherical wavelet basis functions, given a set of discrete observations
% on the sphere.
%
% INPUT
%   dlon    longitude of data points
%   dlat    latitude of data points
%   ax0     lon-lat box contained desired region
%   qparm   spherical grid parameters
%        1. qmin
%        2. qmax
% OUTPUT
%   spline_tot = [glon glat q]
%   glon    longitude of gridpoint for basis function
%   glat    latitude of gridpoint for basis function
%   q       grid index (higher q is desner grid)
%
% calls wavelet_thresh.m
% called by sphereinterp.m
%

function spline_tot = sphereinterp_grid(dlon,dlat,ax0,qparm)

disp('------------------------------------------------------------');
disp('entering sphereinterp_grid.m to obtain spherical wavelet basis functions');

earthr = 6371*1e3;      % earth radius (m)
deg = 180/pi;
msize = 4^2;

qmin = qparm{1};
qsec = qparm{2};
qmax = qparm{3};
ntrsh = qparm{4};

lonmin = ax0(1); lonmax = ax0(2);
latmin = ax0(3); latmax = ax0(4);
ndata = length(dlon);

%-------------------------------
% USER PARAMETERS

slabel = 'sphereinterp_grid.m';

% min number of gridpoints for a particular order
%nmin = 1;

% threshold wavelets based on data
%ntrsh = 3;       % KEY: number of evaluations >= qtrsh

% q for the "secular field" -- which combines the estimates from qmin to
% qsec for the multiscale analysis (NOT VERY RELEVANT HERE)
%qsec = round(mean([qmin qmax]));

%-------------------------------   

% compute nominal scalelength of the latlon square of observations
% --> Lscale is the diameter of a circle with area Atot
Aunit = 4*pi * areaquad(ax0(3),ax0(1),ax0(4),ax0(2));
Atot = earthr^2 * Aunit;
%Atot = latlon_area(ax0,earthr);
Lscale = 2*sqrt(Atot/pi);

% % minimum for dsig -- WHAT IS THIS FOR?
% ineg = find( dsig <= 0 );
% for ii = 1:length(ineg)
%     error('WHY DO WE HAVE NEGATIVE SIGMAS?');
%     jj = ineg(ii);
%     fprintf('%6i%10.4f%10.4f%10.4f%10.4f\n',ii,dlon(jj),dlat(jj),d(jj),dsig(jj));
% end

%========================================================
% PROPERTIES OF THE GRIDS AND THE SPHERICAL WAVELETS

q0_vec = [0:12]';

% scalelength (in degrees) for each grid (spline_wang_A.m)
q_scale_deg = [
    63.43494882292201
    31.71747441146100
    15.85873720573050
    7.92936860286525
    3.96468430143263
    1.98234215071631
    0.99117107535816
    0.49558553767908
    0.24779276883954
    0.12389638441977
    0.06194819220988
    0.03097409610494
    0.01548704805247 ];

% angular support of the wavelet of scale 1/2^j in radians,
% j = 1,...,8
% The wavelet scale and the grid resolution are related by j = q-1,
% and q must be greater or equal 2;
%ang_support = [1.78468; 1.04664; 0.548118; 0.277342; ...
%    0.139087; 0.0695956; 0.0348043; 0.017403];
% a tighter criterion for the support: first zero crossing
ang_support = [82.4415 47.3103 24.7075 12.4999 6.26815 3.13642 ...
    1.56851 0.784289 0.392149 0.196075 ]'/deg;

% extrapolate a few more scales using the values above (see plot)
inew = 10:12;
ang_support(inew+1) = 10.^polyval(polyfit([0:9]',log10(ang_support),1),inew);
%figure; plot([0:length(ang_support)-1]',log10(ang_support),'.');
%xlabel(' scale q'); ylabel('log10 [ angular support of scale-q wavelet ]'); grid on;

ang_support_meters = ang_support*earthr;
qtrsh_cos_dist = cos(ang_support);

disp('Support of the spherical wavelets:');
disp('   q      deg           km');
for ix = 1:length(q0_vec)
    disp(sprintf('%4i %10.3f %10.1f',q0_vec(ix),ang_support(ix)*deg,ang_support_meters(ix)/1000));
end
disp('  ');

%---------------------------------------

% the lowest allowable grid order of the basis functions is taken to be
% one whose support is less than twice the length of the length scale
% of the region of observations
ia = find(ang_support_meters < 2*Lscale);
qmin0 = q0_vec(ia(1));
disp(sprintf('minimum allowable grid order is %i',qmin0));
disp(sprintf('   %.2e meters (support of q = %i wavelet) < %.2e meters (2*Lscale)',...
ang_support_meters(qmin0+1),qmin0,2*Lscale));
%qmin = input([' Type min allowable grid order, qmin >= 0 (try ' num2str(qmin0) '): ']);

if qmin < qmin0
    disp(sprintf('WARNING: qmin (%i) < qmin0 (%i)',qmin,qmin0));
    disp(sprintf('consider using qmin = %i',qmin0));
end

% user picks the max allowable grid (finest scale basis functions)
%qmax = input(' Type max allowable grid order, qmax: ');

Dscale_deg = q_scale_deg(qmax+1);

% KEY: get the gridpoints
[glon,glat,gq,nvec,axmat] = getspheregrid(ax0,qmin,qmax);
qvec = unique(gq);
numq = length(qvec);
spline_tot0 = [glon glat gq];
ax1 = axmat(1,:);

% % load the bounds
% ww = 'subgrid_bounds';
% load([dir1 ww '.dat']); temp = eval(ww);
% ax0 = [temp(1,1) temp(1,2) temp(1,3) temp(1,4)];
% ax0 = ax1;
% 
% % load the number of gridpoints in the region for each order grid
% temp = load([dir1 'num_gridpoints.dat']);
% 
% % threshold the set of gridpoints
% iqs = find( and( temp(:,1) >= qmin, temp(:,1) <= qmax) );
% temp = temp(iqs,:);
% ins = find( temp(:,2) >= nmin );
% temp = temp(ins,:);
% qvec = temp(:,1);
% nvec = temp(:,2);
% disp('  '); disp('           q          num'); disp([qvec nvec]);

% obtain index ranges for each q-level
ngrid = sum(nvec);
id    = cumsum(nvec)';
%ifr   = [1 id(1:end-1)+1 1; id ngrid]';
ifr   = [1 id(1:end-1)+1; id]';

% numq = length(qvec);
% spline_tot0 = zeros(ngrid,3);
% figure; hold on;
% disp('  '); disp(' load the gridpoints in this region');
% for iq = 1:numq
%     q = qvec(iq);
%     n = nvec(iq);
%     stit = [' q = ' num2str(q) ' : num = ' num2str(n)]; disp(stit);
% 
%     %ww = ['thph_q' num2str(sprintf('%2.2i',q))];
%     %load([dir1 ww '.dat']); temp = eval(ww);
%     %lon = temp(:,2)*deg; lat = (pi/2-temp(:,1))*deg;
%     plot(lon,lat,'k.');
% 
%     % fill a matrix with gridpoints
%     inds = [ifr(iq,1) : ifr(iq,2)];
%     spline_tot0(inds,:) = [lon lat q*ones(nvec(iq),1) ];
% end
% axis equal, axis(ax0);

% for iq = 1:numq
%     figure;
%     q = qvec(iq);
%     n = nvec(iq);
%     stit = [' q = ' num2str(q) ' : num = ' num2str(n)]; disp(stit);
%     inds = find(gq==q);
%     plot(glon(inds),glat(inds),'k.');
%     title(stit); axis(ax1);
% end

%========================================================
% THRESHOLD INITIAL SET OF BASIS FUNCTIONS

disp('  '); disp(' threshold the gridpoints');
[ikeep, inum] = wavelet_thresh(spline_tot0, qtrsh_cos_dist, ntrsh, dlon, dlat);

spline_tot = spline_tot0(ikeep,:);

%--------------------------------------------------------
% NOTE: The rest is for display purposes only, but some of the indexing
%       vectors may be needed later for plotting.

ngrid = length(spline_tot);
if ngrid==0, error('ngrid = 0: check datapoints and gridpoints'); end
glon = spline_tot(:,1);
glat = spline_tot(:,2);
gq = spline_tot(:,3);

% recompute qvec, nvec, ifr
jj = 0;
%nvec = zeros(qmax+1,1);
nvec = []; qvec = [];
for q = 0:qmax      % ALL possible q orders
    n = length( find( gq == q ) );
    if n >= 1
        jj = jj+1; qvec(jj) = q; nvec(jj) = n;
        qmax0 = q;
    end
end

% if there are no allowable splines in for qmax, then adjust qmax
qmax = qmax0;
qvec = qvec(:); nvec = nvec(:);
numq  = length(qvec);

id    = cumsum(nvec);
ifr   = [ [1 ; id(1:end-1)+1] id];
%inz = find(nvec==0); ifr(inz,:) = 0;
disp('  '); disp('     q   num    id    i1    i2'); disp([qvec nvec id ifr]);

% inds for the highest q-level
inds_qmax = [ifr(end,1) : ifr(end,2)];
%iz = min(inz)-1;
%inds_qmax = [ifr(iz,1) : ifr(iz,2)];

% select the multiscale decomposition based on what you designate to be the
% 'secular field' associated with the plate convergence
% KEY : qvec --> iqvec, ifr --> iqr
% NOTE: the first row (iqvec, iqr) constitutes ALL available splines
%qsec = 5;
%qsec  = input([' Enter max q grid for secular field (' ...
%    num2str(qmin) ' <= qsec <= ' num2str(qmax) '): ']);
iqsec = find(qvec == qsec);
irecs = find(qvec > qsec);
iqvec = [qvec(1) qvec(end) ; qvec(1) qvec(iqsec) ; [qvec(irecs) qvec(irecs)]];
nump = length(iqvec);
iqr = [1 ifr(end,2) ; 1 ifr(iqsec,2) ; ifr(iqsec+1:end,:)];
ipran = [qsec ; qvec(irecs)];

stqran = ['q = ' num2str(qvec(1)) ' to ' num2str(qvec(end)) ' (' num2str(ngrid) ')'];
strsh1 = [num2str(ngrid) ' wavelets / ' num2str(length(spline_tot0)) ' total'];
strsh2 = [' with >= ' num2str(ntrsh) ' stations inside their corresponding spatial supports'];
disp('  '); disp(['GRIDPOINTS '  stqran ':']);
disp([strsh1 strsh2]); disp('  ');

nvec = zeros(nump,1);
for ip = 1:nump
    q1 = iqvec(ip,1); q2 = iqvec(ip,2);
    if1 = iqr(ip,1);  if2 = iqr(ip,2);
    if if1+if2 > 0, nvec(ip) = if2-if1+1; end
    stqtag{ip} = sprintf('q%2.2i_q%2.2i',q1,q2);
    stqs{ip} = [' q = ' num2str(q1) ' to ' num2str(q2)];
    stis{ip} = [' j = ' num2str(if1) ' to ' num2str(if2) ' (' num2str(nvec(ip)) ')'];
    stit = [stqs{ip} ',' stis{ip}]; disp(stit);
end

% indexing for multi-scale strain and multi-scale residual field
id2 = cumsum(nvec(2:end));
iqr2 = [ ones(nump-1,1) id2];
iqvec2 = [qvec(1)*ones(length(ipran),1) ipran];
%iqvec2 = [qvec(1) qvec(iqsec) ; [qvec(1)*ones(length(irecs),1) qvec(irecs)]];

% plot ALL spline gridpoints
figure; hold on;
scatter(glon,glat,msize,gq,'filled');
%scatter(glon,glat,msize,'ko');
colorbar('ytick',qvec);
%plot(glon,glat,'.');
%plot(ax1([1 2 2 1 1]), ax1([3 3 4 4 3]), 'k');
plot(ax0([1 2 2 1 1]), ax0([3 3 4 4 3]), 'k');
axis equal, axis tight;
title({[slabel ' :  ' stqran],[strsh1 strsh2]},'interpreter','none');
xlabel(' Latitude'); ylabel(' Longitude');
orient tall, wysiwyg, fontsize(10)

% plot spline gridpoints by order q
if nump==2
    figure; hold on;
    inds = [iqr(1,1) : iqr(1,2)];
    plot(dlon,dlat,'k+');
    scatter(glon(inds),glat(inds),msize,gq(inds),'filled');
    axis(ax0);
    title({stqs{ip}, stis{ip}});
else
    figure; nc=2; nr=ceil(nump/nc);    
    for ip = 1:nump
        inds = [iqr(ip,1) : iqr(ip,2)];
        subplot(nr,nc,ip); hold on;
        plot(dlon,dlat,'k+');
        if sum(inds) > 0
            %plot(glon(inds),glat(inds),'.');
            scatter(glon(inds),glat(inds),msize,gq(inds),'filled');
        end
        if qvec(end)-qvec(1) > 0
            caxis([qvec(1) qvec(end)]); colorbar('ytick',qvec);
        end
        axis equal
        axis(ax0);
        %axis(ax1);
        title({stqs{ip}, stis{ip}});
    end
end
orient tall, wysiwyg, fontsize(9)

%========================================================
