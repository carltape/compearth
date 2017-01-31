function [Vhatp,Vhat,omegacrit_deg] = plot_Vomega(gamma_deg,delta_deg)
%PLOT_VOMEGA plot V'(omega) for a fixed lune point (source type)
%
% INPUT
%   gamma_deg       lune longitude, degrees
%   delta_deg       lune latitude, degrees
%
% OUTPUT
%   Vhatp           Vhat'(omega)
%   Vhat            Vhat(omega)
%   omegacrit_deg   critical angles
%
   
if length(gamma_deg)~=1, error('only 1 gamma point allowed'); end
if length(delta_deg)~=1, error('only 1 delta point allowed'); end

deg = 180/pi;

% GENERAL OPTIONS
bincludecrit = false;   % include the most accurate Vhatp and Vhat values
                        % for the critical points

% PLOTTING OPTIONS
bplotrad = true;       % plot x-axis in radians, which is useful for representing
                        %    probability density functions (equal area)
bplotcritpoints = true; % plot critical angles w1, w2, w3, w4
bplotcritlines = false; % 
bplotV = false;         % plot Vhat(omega) in addition to Vhat'(omega)
%Vpstyle = 'b.-';
Vpstyle = 'b-';
Vstyle = 'r-';
msize = 6;              % marker size for critical angles
lwid = 1.5;             % line width for curves
xlims = [-1 181]/deg;

% this will return the default precomputed Vhat'(omega) curves
% NOTE: This will allow for interpolation with gamma but NOT interpolation
%       with omega (for that, just use Vomega.m separately).
[Vhatp,Vhat,omega_rad,omegacrit] = Vomega(gamma_deg/deg,delta_deg/deg);

%Vhatp = real(Vhatp);
%Vhat = real(Vhat);

if any([bplotcritpoints bplotcritlines bincludecrit])
    %otemp = omegacrit{1};
    otemp = omegacrit;
    omegacrit_deg = otemp(:,1);
    Vpcrit = otemp(:,2);
    Vcrit = otemp(:,3);
    %ncrit = length(omegacrit_deg);
    
    % add in the critical angles
    if bincludecrit
        % Vhat'
        vtemp1 = [omegap_rad ; omegacrit_deg/deg];
        vtemp2 = [Vhatp ; Vpcrit];
        vtemp = sortrows([vtemp1 vtemp2],1);
        omegap_rad = vtemp(:,1);
        Vhatp = vtemp(:,2);
        % Vhat
        vtemp1 = [omega_rad ; omegacrit_deg/deg];
        vtemp2 = [Vhat ; Vcrit];
        vtemp = sortrows([vtemp1 vtemp2],1);
        omega_rad = vtemp(:,1);
        Vhat = vtemp(:,2);
    end
end

% plot Vhat'(omega)
plot(omega_rad,Vhatp,Vpstyle,'linewidth',lwid);
hold on;

if bplotV
    plot(omega_rad,Vhat,Vstyle,'linewidth',lwid);
    if bplotcritlines
        [xmat,ymat] = horzlines(Vcrit,xlims(1),xlims(2));
        plot(xmat,ymat,'k--');
    end
    if bplotcritpoints
        plot(omegacrit_deg/deg,Vcrit,'ko','markerfacecolor','r','markersize',msize);
    end
end

if bplotcritlines
    axx = axis; yran = axx(3:4);
    %yran = [0 2.5];
    [xmat,ymat] = vertlines([omegacrit_deg/deg],yran(1),yran(2));
    plot(xmat,ymat,'k--');
end
if bplotcritpoints
    plot(omegacrit_deg/deg,Vpcrit,'ko','markerfacecolor','c','markersize',msize);
    %plot(omegacrit_deg/deg,Vpcrit,'ro','markersize',msize);
end

%axis equal, axis(ax0);
%xlim([-1 max(omegacrit_deg)+1]/deg);
if bplotrad
    xlabel('\omega, radians');
    %xticks = [0:30:180]/deg;
    %set(gca,'xtick',xticks);
    %axis equal, axis tight         % PRESERVE AREA
else
    xlim(xlims);
    xlabel('\omega, degrees');
    xticks = [0:30:180];
    set(gca,'xtick',xticks/deg,'xticklabel',numvec2cell(xticks));
end

title(sprintf('\\gamma = %.1f deg,  \\delta = %.1f deg',gamma_deg,delta_deg));
ylabel('Vhat''(\omega)');

%==========================================================================
