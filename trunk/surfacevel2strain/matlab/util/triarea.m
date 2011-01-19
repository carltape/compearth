%
% function tri = triarea(v1, v2, v3, isphere)
% Carl Tape, 04-Feb-2003
% 
% INPUT
%   v1,v2,v3   input UNIT vectors, ORDERED counter-clockwise
%   isphere    = 1: arc distance and spherical areas
%              = 0: chord distances and planar areas
% OUTPUT
%   tri        area of triangle described by v1,v2,v3
%
% EXAMPLE: triarea([0 0 1],[1 0 0],[0 1 0],1) - 4*pi/8
%
% SEE MATLAB FUNCTION areaint IN EXAMPLE BELOW.
%
% calls xxx
% called by latlon_area.m, tggriddims.m
%

function tri = triarea(v1, v2, v3, isphere)

% probably we want to ensure that v1,v2,v3 are ordered unit vectors

if isphere == 1
    % a,b,c are the arc-sides of the triangles
	a = acos(dot(v3, v2));
	b = acos(dot(v1, v2));
	c = acos(dot(v1, v3));
	s = (a + b + c)/2;
    
	tri = 4*atan(sqrt(tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)));
    
    %disp(['     a          b        c        s       tri]);
    %disp([a b c s tri]);
    
else
    x = v3 - v2;
    y = v1 - v2;
    tri = 0.5 * norm(cross(x,y));
end

%------------------------------------------------------------------

if 0==1
    triarea([0 0 1],[1 0 0],[0 1 0],1) - 4*pi/8
    
    % Matlab function areaint --
    % Accuracy of the integration method is inversely proportional
    % to the distance between lat/long points.
    nvec = 2.^[1:10]';
    numn = length(nvec);
    err_matlab = zeros(numn,1);
    for ii = 1:numn
        n = nvec(ii);
        
        % lat-lon points describing an octahedral piece of the sphere
        lat = [linspace(90,0,n) linspace(0,0,n) linspace(0,90,n)];
        lon = [linspace(0,0,n) linspace(0,90,n) linspace(90,90,n)];
        err_matlab(ii) = 1/8 - areaint(lat,lon);
    end
    figure; loglog(nvec,abs(err_matlab),'.'); grid on;
    xlabel(' number of points on a side of the triangular patch');
    ylabel(' error in the computed area of the patch');
end

%=====================================================================
