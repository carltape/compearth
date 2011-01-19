function [xo_mat,yo_mat,i_zone]=utm2ll(xi_mat,yi_mat,i_zone,i_type,datum)
%Convert utm to lon-lat and vice versa
%   [XO,YO,I_ZONE]=utm2ll(XI,YI,I_ZONE,I_TYPE) converts points
%   (XI,YI) in utm zone I_ZONE to points (XO,YO).  XI and YI may be
%   matrices, vectors, or scalars.
%   I_TYPE = 1 -> lon,lat to  utm
%   I_TYPE = 2 -> utm     to  lon,lat
%   
%   [XO,YO,I_ZONE]=utm2ll(XI,YI,I_ZONE,I_TYPE,DATUM) uses
%   reference model DATUM.
%   DATUM = 1 -> NAD27 [default]
%   DATUM = 2 -> WGS84

%****************************************************************
%**   
%**   FILE NAME: utmtoll
%**   
%**   DATE WRITTEN:7/22/93 
%**   
%**   PROGRAMMER:Scott Hensley
%**   
%**   FUNCTIONAL DESCRIPTION: This routine converts between lat
%**   lon and utm coordinates for a datum determined from the input 
%**   a and e2.
%**   
%**   ROUTINES CALLED:none
%**   
%**   NOTES: none
%**   
%**   UPDATE LOG:
%**   MATLAB version, Sarah Minson 06/07/06   
%****************************************************************
        
%     INPUT VARIABLES:
%        real*8 xi,yi                     !input coordinates
%        integer i_zone                   !UTM zone (needed if i_type=2)
%        integer i_type                   !1=lat,lon to utm,2= utm to lat,lon
% old stuff:
%        real*8 r_a                       !ellispoid semi-major axis
%        real*8 r_e2                      !ellipsoid eccentricity squared  
%        real*8 r_v(2)                    !geocentric vector (meters)
%        real*8 r_lat                     !latitude (deg -90 to 90)
%        real*8 r_lon                     !longitude (deg -180 to 180)
%        character*1 a_grid               !UTM North-South grid
   
%       OUTPUT VARIABLES:
%        real*8 xo,yo                     !output coordinates

%       LOCAL VARIABLES:
%        integer i_ft,i_gi
%        real*8 pi,r_dtor
%        real*8 r_ep2,r_k0,r_k
%        real*8 r_fe,r_fn(2)
%        real*8 r_e4,r_e6,r_n,r_t,r_t2,r_c,r_c2,r_ba
%        real*8 r_a2,r_a3,r_a4,r_a5,r_a6 
%        real*8 r_d,r_d2,r_d3,r_d4,r_d5,r_d6
%        real*8 r_lon0,r_lat1,r_m,r_m0,r_mu,r_lat0
%        real*8 r_et,r_e1,r_e12,r_e13,r_e14,r_r,r_vu(2)
%        character*1 a_griddes(20)

%       DATA STATEMENTS:
        r_dtor=1.74532925199d-2;
        i_ft=0;
        a_griddes=['C','D','E','F','G','H','J',...
            'K','L','M','N','P','Q','R','S','T','U',...
            'V','W','X'];
    if nargin < 5; datum=1; end;
    if datum==1
      r_a=6378206.4d0; %equatorial radius, NAD27
	  r_e2=0.00676865799761d0; %ellipsoid eccentricity squared, NAD27
	else
	  r_a=6378137.0; %equatorial radius, WGS84
	  r_pole=6356752.3; %polar radius, WGS84
	  r_e2=(1-r_pole^2/r_a^2); %ellipsoid eccentricity squared 
	end
        r_k0=0.9996d0;    %scale at center 
        r_lat0=0.0d0;
        r_fe=5e5;
        r_fn(1)=0;
        r_fn(2)=1e7;  % N-S hemispheres??

%       PROCESSING STEPS:

        r_ep2 = r_e2/(1.d0 - r_e2);
        r_e4 = r_e2^2;
        r_e6 = r_e2^3;
        r_dtor = pi/180;

    [nrows,ncols] = size(xi_mat);
    xo_mat = zeros(nrows,ncols); yo_mat=zeros(nrows,ncols);
    for ii=1:nrows*ncols
        xi = xi_mat(ii); yi=yi_mat(ii);    
        
        if i_type == 1   %convert lat,lon to UTM
           
           xi=xi*r_dtor;
           yi=yi*r_dtor;
	   if (i_zone==0)
	     i_zone = fix(mod(xi+3.d0*pi,2.d0*pi)/(r_dtor*6.d0)) + 1;
	     i_zone = max(min(i_zone,60),1);
	   end
           r_lon0 = -pi + 6.d0*r_dtor*(i_zone-1) + 3.d0*r_dtor; % central meridian
           
           r_n = r_a/sqrt(1.d0 - r_e2*sin(yi)^2);
           r_t = tan(yi)^2;
           r_t2 = r_t^2;
           r_c = r_ep2*cos(yi)^2;
           r_ba = (xi - r_lon0)*cos(yi);
           r_a2 = r_ba^2;
           r_a3 = r_ba*r_a2; 
           r_a4 = r_ba*r_a3;
           r_a5 = r_ba*r_a4;
           r_a6 = r_ba*r_a5;
           r_m = r_a*((1.0d0-r_e2/4 - 3.0d0*r_e4/64.0d0-5.d0*r_e6/256.d0)*yi - (3.d0*r_e2/8.d0 + 3.d0*r_e4/32.d0 +45.d0*r_e6/1024.d0)*sin(2.d0*yi) +  (15.d0*r_e4/256.d0 +45.d0*r_e6/1024.d0)*sin(4.d0*yi) - (35.d0*r_e6/3072.d0)*sin(6.d0*yi));
           r_m0 = r_a*((1.d0-r_e2/4 - 3.d0*r_e4/64.d0 -5.d0*r_e6/256.d0)*r_lat0 - (3.d0*r_e2/8.d0 + 3.d0*r_e4/32.d0 +45.d0*r_e6/1024.d0)*sin(2.d0*r_lat0) +  (15.d0*r_e4/256.d0 +45.d0*r_e6/1024.d0)*sin(4.d0*r_lat0) - (35.d0*r_e6/3072.d0)*sin(6.d0*r_lat0));
           
           r_v(1) = r_k0*r_n*(r_ba+(1.d0-r_t+r_c)*r_a3/6.d0 +(5.d0-18.d0*r_t+r_t2+72.d0*r_c-58.d0*r_ep2)*r_a5/120.d0);
           r_v(1) = r_v(1) + r_fe;

           r_v(2) = r_k0*(r_m - r_m0 + r_n*tan(yi)*( r_a2/2.d0 + (5.d0-r_t+9.d0*r_c+4.d0*r_c^2)*(r_a4/24.d0) + (61.d0-58.d0*r_t+r_t2+600.d0*r_c-330.d0*r_ep2)*(r_a6/720.d0) ));
%           if yi >= 0
%              r_v(2) = r_v(2) + r_fn(1)
%           else
%              r_v(2) = r_v(2) + r_fn(2)
%           end

           r_k = r_k0*(1.d0+(1.d0+r_ep2*cos(yi)^2)*(r_v(1)-r_fe)^2/(2.d0*(r_k0^2)*r_n^2));

           i_gi = fix((yi/r_dtor+80.d0)/8.d0) + 1;
           i_gi = max(min(i_gi,20),1);
           a_grid = a_griddes(i_gi);

           xo=r_v(1);
           yo=r_v(2);
           
        elseif i_type == 2         %convert UTM to lat,lon 

           r_vu(1) = xi - r_fe;
           r_vu(2) = yi;
%           if a_grid < 'M'
%              r_vu(2) = yi - r_fn(2)
%           end
           r_lon0 = -pi + 6.d0*r_dtor*(i_zone-1) + 3.d0*r_dtor;
           
           r_et = sqrt(1.d0-r_e2);
           r_e1 = (1.d0-r_et)/(1.d0+r_et);
           r_e12 = r_e1^2;
           r_e13 = r_e1*r_e12;
           r_e14 = r_e1*r_e13;
           r_m = r_vu(2)/r_k0;
           r_mu = r_m/(r_a*(1.d0-r_e2/4.d0-3.d0*r_e4/64.d0-5.d0*r_e6/256.d0));
           r_lat1 = r_mu + (3.d0*r_e1/2.d0-27.d0*r_e13/32.d0)*sin(2.d0*r_mu)+(21.d0*r_e12/16.d0-55.d0*r_e14/32.d0)*sin(4.d0*r_mu) +(151.d0*r_e13/96.d0)*sin(6.d0*r_mu) +(1097.d0*r_e14/512.d0)*sin(8.d0*r_mu);

           r_n = r_a/sqrt(1.d0 - r_e2*sin(r_lat1)^2);
           r_r = (r_a*(1.d0-r_e2))/sqrt(1.d0 - r_e2*sin(r_lat1)^2)^3;
           r_t = tan(r_lat1)^2;
           r_t2 = r_t^2;
           r_c = r_ep2*cos(r_lat1)^2;
           r_c2 = r_c^2;
           r_d = r_vu(1)/(r_n*r_k0);
           r_d2 = r_d^2;
           r_d3 = r_d2*r_d;
           r_d4 = r_d3*r_d;
           r_d5 = r_d4*r_d;
           r_d6 = r_d5*r_d;
 
           yo = r_lat1 - (r_n*tan(r_lat1)/r_r)*(r_d2/2.d0+(5.d0+3.d0*r_t+10.d0*r_c-4.d0*r_c2-9.d0*r_ep2)*r_d4/24.d0 +(61.d0+90*r_t+298.d0*r_c+45.d0*r_t2-252.d0*r_ep2-3.d0*r_c2)*(r_d6/720.d0));
           xo = r_lon0 + (r_d - (1.d0+2.d0*r_t+r_c)*r_d3/6.d0 +(5.d0-2.d0*r_c+28.d0*r_t-3.d0*r_c2+8.d0*r_ep2+24.d0*r_t2)*(r_d5/120.d0))/cos(r_lat1);
           xo=xo/r_dtor;
           yo=yo/r_dtor;

        end
        
        xo_mat(ii)=xo; yo_mat(ii)=yo;
    end
        
        
