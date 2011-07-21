%
% run_getspheregrid.m
% Carl Tape, 01-Jan-2011
%
% Matlab version of Fortran code getsubgrids.f90 associated with the
% software package SURFACEVEL2STRAIN (Tape et al., 2009, GJI).
%
% calls getspheregrid.m
% called by xxx
%

clc
clear
close all
format short, format compact

deg = 180/pi;
axf = [-180 180 -90 90];

% input
%ax0 = [-180 180 -90 90]; qmin = 0; qmax = 3;
%ax0 = [-180 180 -90 90]; qmin = 4; qmax = qmin;
%ax0 = [145 200 40 70]; qmin = 4; qmax = qmin;   % alaska
%ax0 = [185 240 40 70]; qmin = 4; qmax = qmin;   % alaska
%ax0 = [-122 -113 30 38]; qmin = 7; qmax = qmin;   % socal
ax0 = [-122 -113 30 38]; qmin = 0; qmax = 6;   % socal
%ax0 = [-128 -100 30 51]; qmin = 0; qmax = 7;   % pacific

if any(ax0(1:2) > 180), ilon360=1; else ilon360=0; end
if ilon360==1
    axf = [0 360 -90 90]; lontick = [0:60:360];
else
    axf = [-180 180 -90 90]; lontick = [-180:60:180];
end
lattick = [-90:30:90];

[glon,glat,gq,nvec,axmat] = getspheregrid(ax0,qmin,qmax);
qvec = unique(gq);
nq = length(qvec);

disp('number of gridpoints for each q');
disp([ [qmin:qmax]' nvec]);
disp('lon-lat bounds for each q:');
disp(axmat);

%dlon = glon*deg;
%dlat = 90-glat*deg;

figure;
nc=2; nr=ceil(nq/nc);
if nq==1, nr=1; nc=1; end
for ii = 1:nq
    q = qvec(ii);
    inds = find(gq==q);
    subplot(nr,nc,ii); hold on;
    plot(glon(inds), glat(inds), '.');
    title(sprintf('q = %i, n = %i/%i',q,length(inds),length(glon)));
    plot(ax0([1 2 2 1 1]),ax0([3 3 4 4 3]),'r');
    %axis(ax0);
    set(gca,'xtick',lontick,'ytick',lattick); axis(axf); grid on;
end

%=========================================================================
