%
% function [dh,dv] = fault_2D(x,slip,W,dip,depth,fault_type)
%
% analytical solution for infinitely long fault in a homogeneous 1/2 space
%
% dh            = horizontal displacement
% dv            = vertical displacement (positive down)
% x             = horizontal position
% slip          = slip
% dip           = dip in degrees
% W             = down dip fault length
% depth         = depth to top of fault
% fault_type    = 1 for strike slip, 2 for dip slip
%
% dh is displacement in horizontal direction
% perpendicular to strike for dip slip and 
% parallel to strike for strike slip
% 
% NOTICE THE AXES CONVENTION:
%   x is east
%   y is south
%   z is down
% fault is centered on x=0 and strikes in the +- y direction
% dip is clockwise from positive x axis
%
% NOTICE THE UNITS:
% W, depth must be in same units as x(e.g., km)
% output (dh, dv) in same units as slip (e.g., cm)
% 
% fault_type = 1 for strike slip, 2 for dip slip
% for dip slip fault, slip > 0 (0 < dip < 90) is THRUST fault
%
% depth is measured at top of fault
% 
% does not appear to be robust for handling x as row or column matrix
%
% currently only handles dip=90 for ss faults
%
% calls xxx
% called by xxx
%

function [dh,dv] = fault_2D(x,slip,W,dip,depth,fault_type)

dip = dip*pi/180;

if(depth<0)
    disp('depth must be a positive quantity');
    return
end
 
if(fault_type == 1)
	d1 = W+depth;
	d2 = depth;
	dh = -slip*(atan(x/d1)-atan(x/d2))/pi;
	dv = zeros(size(dh));
	return
end

%factor of pi in second term of each U?
if(fault_type == 2)
	a = depth/sin(dip);
	D1 = x.^2 + (W+a)^2 - 2*x*(W+a)*cos(dip);
	D2 = x.^2 + a^2 - 2*x*a*cos(dip);
    
	%horizontal
	dh = (slip/pi) * ...
     ( (W+a)*sin(dip)*((W+a)-x*cos(dip))./D1 + ...
     cos(dip)*atan( (x-(W+a)*cos(dip))/((W+a)*sin(dip)) ) - ...
     a*sin(dip)*(a-x*cos(dip))./D2 - ...
     cos(dip)*atan( (x-a*cos(dip))/(a*sin(dip)) ) );
        
	%vertical
	dv = (slip*sin(dip)/pi) * ...
     ( (x*(W+a)*sin(dip)./D1) + ...
       atan( (x-(W+a)*cos(dip))/((W+a)*sin(dip)) ) - ...
       (x*a*sin(dip)./D2) - ...
       atan( (x-a*cos(dip))/(a*sin(dip)) ) );
	
	return
end
  
