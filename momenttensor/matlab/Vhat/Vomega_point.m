function [Vhatp,Vhat] = Vomega_point(gamma_deg,delta_deg,omega_deg,atol)
%VOMEGA fractional volume as a function of omega for fixed lune point
%
% INPUT
%   gamma_deg       lune longitude, in degrees
%   delta_deg       lune latitude, in degrees
%   omega_deg       angle between moment tensors, in degrees
%
% OUTPUT
%   Vhatp           Vhat'(omega)
%   Vhat            Vhat(omega)
%
% See examples in run_Vomega.m
%
% calls Vgammaomega.m
%

if length(gamma_deg)~=1, error('only 1 gamma point allowed'); end
if length(delta_deg)~=1, error('only 1 delta point allowed'); end
if nargin==3, atol = 1e-12; end
n = length(omega_deg);

deg = 180/pi;
%DOMEGA = atol*1e2;  % KEY COMMAND (depends on tolerances in Vgammaomega.m)
DOMEGA = 1e-6;

% radians
gamma = gamma_deg/deg;
delta = delta_deg/deg;
omega = omega_deg/deg;

Vhat  = NaN(n,1);
Vhatp = NaN(n,1);

disp(sprintf('Vomega_point.m: atol = %.1e, domega = %.2e radians',atol,DOMEGA));
for ii=1:n
    otemp = omega(ii);
    disp(sprintf('%i/%i: omega = %.2f deg',ii,n,otemp*deg));
    Vhat(ii) = Vx(gamma,delta,otemp,atol);
    
    % centered on omega (2 additional evaluations)
    % note: this is the most appropriate formula to use for the derivative,
    % though it does require two additional function evaluations
    Vhatp(ii) = ( Vx(gamma,delta,otemp+DOMEGA/2,atol) - ...
                  Vx(gamma,delta,otemp-DOMEGA/2,atol) ) / DOMEGA;
              
    % to the right (1 additional evaluation)
    %Vhatp(ii) = ( Vx(gamma,delta,otemp+DOMEGA,atol) - ...
    %              Vhat(ii) ) / DOMEGA;
    
    % to the left (1 additional evaluation)
    %Vhatp(ii) = ( Vhat(ii) - ...
    %             Vx(gamma,delta,otemp-DOMEGA,atol) ) / DOMEGA;
end

%==========================================================================

function V = Vx(gamma,delta,omega,atol)

beta = pi/2 - delta;
% Tape and Tape (2016), Eq 44
omegadev = 2 * asin( sin(omega/2) / sin(beta) );
% Tape and Tape (2016), Eq 45
V = Vgammaomega(gamma,omegadev,atol);

%==========================================================================
