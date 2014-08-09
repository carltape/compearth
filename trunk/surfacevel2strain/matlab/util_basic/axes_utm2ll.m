function ax_out = axes_utm2ll(ax_in,s_zone,i_type,ellipsoid) 
%AXES_UTM2LL convert between a bounding region in lat-lon and in UTM

xmin0 = ax_in(1);
xmax0 = ax_in(2);
ymin0 = ax_in(3);
ymax0 = ax_in(4);

% % if no ellipsoid is given then use Matlab default for that zone
% if nargin == 3
%     disp('no ellipsoid is provided as input');
%     [ellipsoid,estr] = utmgeoid(s_zone)
% end

if nargin==3
    [xmin,ymin] = utm2ll(xmin0,ymin0,s_zone,i_type);
    [xmax,ymax] = utm2ll(xmax0,ymax0,s_zone,i_type);
else
    [xmin,ymin] = utm2ll(xmin0,ymin0,s_zone,i_type,ellipsoid);
    [xmax,ymax] = utm2ll(xmax0,ymax0,s_zone,i_type,ellipsoid);
end

ax_out = [xmin xmax ymin ymax];
     
%==========================================================================

% EXAMPLE
if 0==1
    ax_ll = [-121.6 -114.7 32.2 36.8];
    ax_utm = axes_utm2ll(ax_ll,'11S',0);
    ax_llcheck = axes_utm2ll(ax_utm,'11S',1);
    ax_ll, ax_utm, ax_llcheck
end

%==========================================================================
