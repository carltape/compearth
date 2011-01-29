%
% function A = latlon_area(ax0,srad)
% Carl Tape, 25-July-2008
% 
% This function computes the area of a square patch on the sphere described
% by two bounding latitude lines and two bounding longitude lines.
%
% Matlab function areaquad.m makes latlon_area.m OBSOLETE.
%
% See examples below.
% 
% calls triarea.m
% called by xxx
%

function A = latlon_area(ax0,srad)

if 0==1     % use Matlab built-in function
    Aunit = 4*pi * areaquad(ax0(3),ax0(1),ax0(4),ax0(2));
    A = srad^2 * Aunit;
    
else
    lonmin = ax0(1); lonmax = ax0(2);
    latmin = ax0(3); latmax = ax0(4);

    % order points counter-clockwise
    Plat = [latmin latmin latmax latmax]';
    Plon = [lonmin lonmax lonmax lonmin]';
    Pxyz = latlon2xyz(Plat,Plon,1);

    % divide square into a NW triangular patch and a SE triangular patch
    % --> these are computed for the UNIT SPHERE
    Anw = triarea(Pxyz(:,1), Pxyz(:,3), Pxyz(:,4), 1);
    Ase = triarea(Pxyz(:,1), Pxyz(:,2), Pxyz(:,3), 1);

    A = srad^2 * (Anw + Ase);
end

%--------------------------------------
% EXAMPLES

if 0==1
    format long
    
    % quarter-wedge of the sphere
    srad = 1;
    Atot = 4*pi*srad^2;
    ax0 = [0 90 -89.99 89.99];
    A = latlon_area(ax0,srad) / Atot
    
    % socal tomography region
    srad = 6371*1e3;
    ax0 = [-121.6 -114.7 32.2 36.8];
    A = latlon_area(ax0,srad)
end

%=====================================================================
