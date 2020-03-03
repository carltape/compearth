function [az_deg,dt,Z] = dest2split(d1,d2)
% convert estimated scalar fields back to splitting parameters

Z = d1 + 1i*d2;

% recover dt and splitting angle
R = abs(Z);
Theta_rad = angle(Z);   % radians
theta_rad = Theta_rad / 2;
r = sqrt(R);
dt = r;
th_deg = theta_rad*180/pi;
az_deg = ph2az(th_deg);
% to put the azimuth as -90 to 90
az_deg(find(az_deg > 90)) = az_deg(find(az_deg > 90)) - 180;
