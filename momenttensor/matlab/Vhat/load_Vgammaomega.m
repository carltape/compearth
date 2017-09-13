function [VVp,Wp,Gp,VV,W,G] = load_Vgammaomega(USE_POSITIVE_GAMMA)
%LOAD_VGAMMAOMEGA load pre-saved library of V_gamma(omega) calculations
%
% W. Tape and C. Tape, 2017, GJI
% Volume in moment tensor space in terms of distance
%

% USE_MAT_FILES
% true:  use pre-saved Matlab files that were saved in run_Vhat_lib.m by calling Vgammaomega.m;
%        these include calculation times for each Vhat_gamma(omega).
%        variable names: omega V tocv Vp Vn tocvp tocvn
% false: use pre-saved ASCII files that were written from the Matlab files
%        by calling Vomega.m
USE_MAT_FILES = false;

deg = 180/pi;
if nargin==0
    USE_POSITIVE_GAMMA = true;
end
if USE_POSITIVE_GAMMA
    gtag = 'gammap';
else
    gtag = 'gamman';
end

path_Vhat;
if USE_MAT_FILES
    % base directory
    bdir = strcat(Vhatdir_write,'/data/');
    atag = 'atol10'; rtag = 'rtol18';
    ddir = sprintf('%s%s_iterated_%s_%s/',bdir,gtag,atag,rtag);
else
    bdir = Vhatdir_presaved;
    ddir = sprintf('%sVhat_%s/',bdir,gtag);
end

% omega values
disp('using default omega values from library');
fname = sprintf('%somegavec.dat',ddir);
if ~exist(fname,'file'), error('%s does not exist',fname); end
omegalib_deg = load(fname);
omegalib_rad = omegalib_deg/deg;
nomega = length(omegalib_rad);
nomegap = nomega-1;

% gamma values
fname = sprintf('%sgammavec.dat',ddir);
if ~exist(fname,'file'), error('%s does not exist',fname); end
temp = load(fname);
gamma_deg = temp(:,1);
gamma_rad = gamma_deg/deg;
ngamma = length(gamma_rad);

% Tape and Tape, 2016, Eq 46b // Tape and Tape, 2017, Eq 38
% we only need to consider positive gamma and positive delta
%gamma_rad = abs(gamma_rad);
%delta_rad = abs(delta_rad);

VV = NaN(ngamma,nomega);
VVp = NaN(ngamma,nomega-1);

for ii=1:ngamma
    glib_deg = gamma_deg(ii);
    gtag = sprintf('100gamma_%4.4i',round(100*glib_deg));
    if USE_MAT_FILES
        % run_Vhat_lib.m: omega V tocv Vp Vn tocvp tocvn
        fname = sprintf('%sVhat_%s.mat',ddir,gtag);
        %if ~exist(fname,'file'), error('%s does not exist',fname); end
        if ~exist(fname,'file')
            warning('%s does not exist',fname);
            continue
        end
        load(fname);
        Vlib = real(V(:));
        omegalib_rad = omega(:);

        % create points for evaluating derivative
        % (shift omega points by dw/2, then cut the end point)
        dw = omegalib_rad(2) - omegalib_rad(1);
        %omegalibp_rad = dw/2 + omegalib_rad;
        %omegalibp_rad(end) = [];
        % numerical derivative (note: omega points are uniformly spaced)
        Vlibp  = diff(Vlib) / dw;
    else
        % Vhat
        fname = sprintf('%sVhat_atol10_%s_100delta_0000.dat',ddir,gtag);
        xx = load(fname);
        Vlib = xx(:,2);
        % Vhat-prime
        fname = sprintf('%sVhatp_atol10_%s_100delta_0000.dat',ddir,gtag);
        xxp = load(fname);
        Vlibp = xxp(:,2);
    end
    
    VV(ii,:)  = Vlib(:)';
    VVp(ii,:) = Vlibp(:)';
end

if USE_MAT_FILES
    omegalibp_rad = dw/2 + omegalib_rad;
    omegalibp_rad(end) = [];
else
    omegalib_rad = xx(:,1)/deg;
    omegalibp_rad = xxp(:,1)/deg;
end

%W = repmat( omegalib_rad(:), 1, ngamma );
%G = repmat( gamma_rad(:)', nomega, 1 );
%Wp = repmat( omegalibp_rad(:), 1, ngamma );
%Gp = repmat( gamma_rad(:)', nomega-1, 1 );

W = repmat( omegalib_rad(:)', ngamma, 1 );
G = repmat( gamma_rad(:), 1, nomega );
Wp = repmat( omegalibp_rad(:)', ngamma, 1 );
Gp = repmat( gamma_rad(:), 1, nomega-1 );

%==========================================================================
% EXAMPLES

if 0==1
    clims = [0 1]; deg = 180/pi;
    tic; [Vp,Wp,Gp,V,W,G] = load_Vgammaomega; toc
    %figure; pcolor(W*deg,G*deg,V); shading flat; axis equal; axis([0 180 0 30]); caxis(clims);
    %figure; surf(W*deg,G*deg,V); shading flat;
    figure; pcolor(Wp*deg,Gp*deg,Vp); shading flat; axis equal; axis([0 180 0 30]); caxis(clims);
    figure; surf(Wp*deg,Gp*deg,Vp); shading flat; caxis(clims); axis([0 180 0 30 clims]);
    
    % vectors that parameterize the arrays
    omegalib_rad = W(1,:)';
    omegalibp_rad = Wp(1,:)';
    gamma_rad = G(:,1);
    
    % interpolate to get Vhat_gamma(omega)
    otar = [10:5:110]'/deg;
    gtar = 10/deg;
    interp2(W,G,V,otar,gtar*ones(size(otar)))
    
    % gamma < 0
    USE_POSITIVE_GAMMA = false;
    tic; [Vp,Wp,Gp,V,W,G] = load_Vgammaomega(USE_POSITIVE_GAMMA); toc
    %figure; pcolor(W*deg,G*deg,V); shading flat; axis equal; axis([0 180 0 30]); caxis(clims);
    %figure; surf(W*deg,G*deg,V); shading flat;
    figure; pcolor(Wp*deg,Gp*deg,Vp); shading flat; axis equal; axis([0 180 -30 0]); caxis(clims);
    figure; surf(Wp*deg,Gp*deg,Vp); shading flat; caxis(clims); axis([0 180 -30 0 clims]);
end

%==========================================================================
