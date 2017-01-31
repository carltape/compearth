function [Vhatp,Vhat,otar,omegacrit] = Vomega(gamma_rad,delta_rad,omega_rad)
%function [Vhatp,Vhat,omegalibp_rad,omegalib_rad,omegacrit] = Vomega(gamma_rad,delta_rad,omega_rad)
%VOMEGA fractional volume as a function of omega for fixed lune point
%
% INPUT
%   gamma_rad       lune longitude, in radians
%   delta_rad       lune latitude, in radians
%   omega_rad       OPTIONAL: angle between moment tensors, in radians
%
% OUTPUT
%   Vhatp           set of Vhat'(omega) for each (gamma,delta) [precomputed]
%   Vhat            set of Vhat(omega) for each (gamma,delta) [precomputed]
%   otar            omega angles for precomputed Vhat' and Vhat' curves
%   omegacrit       critical angles for precomputed Vhat curves
%
% SEE EXAMPLES AT BOTTOM.
%
% calls load_Vgammaomega.m, load_omegacrit.m
%
% Carl Tape, 2017-01-19
%

deg = 180/pi;

%==========================================================================
% USER INPUT

% USE_POSITIVE_GAMMA
% true:  use Vhat(omega) files for gamma >= 0 always
% false: use Vhat(omega) files for gamma < 0, if an input gamma < 0 is specified
USE_POSITIVE_GAMMA = true;

EXTRAPVAL = 0;

USE_OMEGA_MIDPOINTS = false;

% WRITE_ASCII_FILES
% true:  write Vhat(omega) and Vhat'(omega) to ascii files (must have USE_MAT_FILES = true)
% false: do not write any ascii files
% To write out a set of text files, try this:
%    gamma = [0:0.5:30]*pi/180; for ii=1:length(gamma), Vomega(gamma(ii),0); end
%    gamma = -[0:0.5:30]*pi/180; gamma(end) = []; for ii=1:length(gamma), Vomega(gamma(ii),0); end
%WRITE_ASCII_FILES = false;

bdir = '/home/carltape/PROJECTS/cmt/Vhat/compearth/';
ddir = sprintf('%sVhat_gammap/',bdir);
%cdir = sprintf('%somega_crit/',bdir);

%==========================================================================

if nargin==2
    disp('using default omega values from library');
    fname = sprintf('%somegavec.dat',ddir);
    if ~exist(fname,'file'), error('%s does not exist',fname); end
    omegalib_radians = load(fname);
    nomega = length(omegalib_radians);
    nomegap = nomega-1;
    if USE_OMEGA_MIDPOINTS
        nomega = nomegap;
    end
else
    disp('input omega values detected');
    nomega = length(omega_rad);
    nomegap = nomega;
    if nargout > 2, error('only 2 output arguments allowed'); end
end

if numel(gamma_rad)~=numel(delta_rad)
    error('gamma and delta must have same number of entries');
end
nlune = length(gamma_rad);

% load library of Vhat (and Vhat-prime) values
[Vp,Wp,Gp,V,W,G] = load_Vgammaomega(USE_POSITIVE_GAMMA);
omegalib_rad = W(1,:);
omegalibp_rad = Wp(1,:);
%gammalib_rad = G(:,1);

[omegacritdev,gamma_deg_crit] = load_omegacrit;

Vhatp = NaN(nlune,nomega);
Vhat  = NaN(nlune,nomega);
omegacrit = cell(nlune,1);

% loop over lune points
for kk=1:nlune
    % lune point
    gamma = gamma_rad(kk);
    delta = delta_rad(kk);
    disp(sprintf('gamma = %.2f deg, delta = %.2f deg',gamma*deg,delta*deg));

    %-------------
    % critical values
    % note: for now we just use the closest gamma, though we ought to be
    %       open plot_interpolating on gamma as well (like we do for Vhat below)
    
    % get closest precomputed Vhat_gamma(omega) (see run_Vhat_lib.m)
    [~,imin] = min(abs(abs(gamma) - gamma_deg_crit/deg));
    disp(sprintf('closest pre-computed gamma (for critical values) is %.2f deg (diff %.2f deg)',...
       gamma_deg_crit(imin), gamma*deg - gamma_deg_crit(imin) ));
    
    % Tape and Tape (2016), Eqs 44, 45, 49
    otemp = omegacritdev{imin};
    ocritdev_rad = otemp(:,1)/deg;
    ocrit_rad = omegadev2omega(delta_rad,ocritdev_rad);
    otemp(:,1) = ocrit_rad*deg;
    otemp(:,2) = otemp(:,2) .* cos(ocrit_rad/2) ./ ...
            sqrt( sin(pi/2-delta_rad(kk)).^2 - sin(ocrit_rad/2).^2 );
    otemp(1,2:3)   = [0 0];
    otemp(end,2:3) = [0 1];
    omegacrit{kk} = otemp;
    
    %-------------
    % Vhat(omega) and Vhat'(omega)
    
    % note that if the input omega values have omega > omega4,
    % the calculation will lead to complex values or Inf
    if nargin==3
        otar = omega_rad(:)';
    else
        if USE_OMEGA_MIDPOINTS
            otar = omegalibp_rad;  % 0.05,0.15,...,179.90,179.95
        else
            otar = omegalib_rad;   % 0,0.1,0.2,...,179.9,180.0
        end
    end
    
    % WARNING: THESE OMEGAS ARE NOT UNIFORMLY SPACED
    omegadev_rad = real( omega2omegadev(delta_rad(kk),otar) );
    %figure; plot(omegadev_rad, '.'); error
    
    % interpolate to get Vhat_gamma(omega) and Vhat'(omega)
    % Tape and Tape (2016), Eq 45
    Vx = interp2(W,G,V,omegadev_rad,gamma*ones(size(omegadev_rad)),'linear');
    Vdevpx = interp2(Wp,Gp,Vp,omegadev_rad,gamma*ones(size(omegadev_rad)),'linear',EXTRAPVAL);
    
    % Tape and Tape (2016), Eq 49
    wvec = cos(otar/2) ./ sqrt( sin(pi/2-delta_rad(kk)).^2 - sin(otar/2).^2 );
    Vpx = Vdevpx .* wvec;
    %if ~isreal(wvec), error('wvec must be real'); end
    
    % we could safeguard the omega values to avoid using real()
    Vhat(kk,:)  = real(Vx);
    Vhatp(kk,:) = real(Vpx);
    
%     if USE_MAT_FILES
%         % save text files of Vhat and Vhatp
%         if and(nargin==2,WRITE_ASCII_FILES)
%             ofile = sprintf('%sVhat_%s_%s_%s.dat',ddir,atag,gtag,dtag)
%             fid = fopen(ofile,'w');
%             for ii=1:length(omegalib_rad)
%                 fprintf(fid,'%20.12f %20.12f\n',omegalib_rad(ii)*deg,Vx(ii));
%             end
%             fclose(fid);
%             ofile = sprintf('%sVhatp_%s_%s_%s.dat',ddir,atag,gtag,dtag)
%             fid = fopen(ofile,'w');
%             for ii=1:length(omegalibp_rad)
%                 fprintf(fid,'%20.12f %20.12f\n',omegalibp_rad(ii)*deg,Vpx(ii));
%             end
%             fclose(fid);
%             %figure; plot(omegalibp_rad*deg,Vpx,'.-');
%         end
%         %figure; plot(omegalibp_rad*deg,Vpx,'.-');
%     end
end

if nlune==1
   omegacrit = omegacrit{1};
end

%--------------------------------------------------------------------------

function omegadev_rad = omega2omegadev(delta_rad,omega_rad)

beta_rad = pi/2 - delta_rad;
% Tape and Tape (2016), Eq 44
omegadev_rad = 2 * asin( sin(omega_rad/2) / sin(beta_rad) );

%--------------------------------------------------------------------------

function omega_rad = omegadev2omega(delta_rad,omegadev_rad)

beta_rad = pi/2 - delta_rad;
% Tape and Tape (2016), Eq 44
omega_rad = 2 * asin( sin(beta_rad) .* sin(omegadev_rad/2) );

%==========================================================================
% EXAMPLES

if 0==1
    %% the mesa (DC)
    figure; plot_Vomega(0,0); axis equal, axis([0 pi 0 1]);
    % deviatoric source types
    gvec = [0:4:28]; delta_deg = 0; figure;
    for xx = 1:length(gvec), plot_Vomega(gvec(xx),delta_deg); end
    axis equal, axis([0 pi 0 2.1]); title('');
    % gamma = 10, vary delta
    gamma_deg = 10; dvec = [0:15:75]; figure;
    for xx = 1:length(dvec), plot_Vomega(gamma_deg,dvec(xx)); end
    axis equal, axis([0 pi 0 5.2]); title('');
    % source types that are near the lune boundary
    gamma_deg = 29; dvec = [0:15:75]; figure;
    for xx = 1:length(dvec), plot_Vomega(gamma_deg,dvec(xx)); end;
    axis equal, axis([0 pi 0 pi]); title('');
    
    %% test interpolation with gamma
    gvec = [26 26.25 26.5]; figure;
    for xx = 1:length(gvec), plot_Vomega(gvec(xx),0); end;
    axis equal, axis([0 pi 0 pi]);
    
    %% densify gamma, use default omega
    dg = 0.1; deg = 180/pi; view3 = [-15,6];
    delta0 = 0; clims = [0 1];
    %delta0 = 40; clims = [0 2];
    gamma_rad = [0:dg:30-dg]/deg;
    delta_rad = delta0/deg*ones(size(gamma_rad));
    [Vhatp,Vhat,otar] = Vomega(gamma_rad,delta_rad);
    [W,G] = meshgrid(otar,gamma_rad);
    figure; surf(W*deg,G*deg,Vhatp); shading flat; caxis(clims); axis([0 180 0 30 clims]); view(view3);
    xlabel('\omega, deg'); ylabel('\gamma, deg'); zlabel('Vhat''(\omega)');
    title(sprintf('\\delta = %.2f deg',delta0));
    % compare with uninterpolated curves
    [Vp,Wp,Gp] = load_Vgammaomega;
    figure; surf(Wp*deg,Gp*deg,Vp); shading flat; caxis(clims); axis([0 180 0 30 clims]); view(view3);
    % now densify both gamma and omega
    otar = [0:0.05:180]/deg;
    tic, [Vhatp,Vhat] = Vomega(gamma_rad,delta_rad,otar); toc
    [W,G] = meshgrid(otar,gamma_rad);
    figure; surf(W*deg,G*deg,Vhatp); shading flat; caxis(clims); axis([0 180 0 30 clims]); view(view3);
end

%==========================================================================
