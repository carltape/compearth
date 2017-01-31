function [omegacrit,gamma_deg] = load_omegacrit
%LOAD_OMEGACRIT load pre-computed critical values omegac, and Vhat(omegac)
% 
% Load pre-computed values of Vhat(omegac) where omegac is a critical
% value for the deviatoric points on the lune.
%
%

deg = 180/pi;
bdir = '/home/carltape/PROJECTS/cmt/Vhat/compearth/';
cdir = sprintf('%somega_crit/',bdir);
ddir = sprintf('%sVhat_gammap/',bdir);

% gamma values
fname = sprintf('%sgammavec.dat',ddir);
if ~exist(fname,'file'), error('%s does not exist',fname); end
temp = load(fname);
gamma_deg = temp(:,1);
ngamma = length(gamma_deg);

% Tape and Tape, 2016, Eq 46b
% we only need to consider positive gamma and positive delta
%gamma_rad = abs(gamma_rad);
%delta_rad = abs(delta_rad);

omegacrit = cell(ngamma,1);

for ii=1:ngamma
    glib_deg = abs(gamma_deg(ii));  % ensure positive gamma
    ctag = sprintf('100gamma_%4.4i',round(100*glib_deg));
    
    % get precomputed critical angles for gamma = 0 (see run_Vhat_lib.m)
    % NOTE THAT THESE ARE NOT THE CRITICAL ANGLES FOR A NONZERO DELTA
    fname = sprintf('%somegacrit_%s.dat',cdir,ctag);
    %if ~exist(fname,'file'), error('%s does not exist',fname); end
    if ~exist(fname,'file')
        warning('%s does not exist',fname);
        continue
    end
    otemp = load(fname);
    omegacrit{ii} = otemp;
end

%==========================================================================
% EXAMPLES

if 0==1
    % critical values in gamma-omega space
    [omegacrit,gamma_deg] = load_omegacrit;
    ngamma = length(omegacrit);
    figure; hold on; axis equal, axis([-1 181 -1 31]);
    for ii=1:ngamma
        otemp = omegacrit{ii};
        plot(otemp(:,1),gamma_deg(ii),'k.');
    end
    
    % V_gamma'(omega) for all critical values
    figure; hold on; axis([-1 181 0 3]);
    for ii=1:ngamma
        otemp = omegacrit{ii};
        plot(otemp(:,1),otemp(:,2),'k.');
    end
    xlabel('\omega_c, deg'); ylabel('Vhat''(\omega_c)');
    fontsize(16)
end

%==========================================================================
