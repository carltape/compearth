function Mw = m02mw(imag,M0)
%M02MW converts from scalar seismic moment to moment magnitude
%
% INPUT
%   imag    =1 for Kanamori (1977), =2 for GCMT formula
%   M0      scalar seismic moment in N-m (not dyne-cm)
%
% See Latex notes cmt.pdf.
%
% EXAMPLE (Denali earthquake):
%    M0 = 7.48*1e20; Mw = m02mw(1,M0)
%
% Carl Tape, 2007-10-31
%

% convert moment tensor from N-m to dyne-cm, since these formulas assume dyne-cm
M0 = 1e7 * M0;

if imag==1
    % Kanamori 1978
    %k = -10.7;
    
    % Kanamori 1977
    %k = -(2/3)*(11.8 - log10(1/2e4));   % for M0 in units of dyne-cm
    %k = -(2/3)*(11.8 - log10(5e-5));    % dyne-cm
    %k = -(2/3)*(16.8 - log10(5));       % dyne-cm
    %k = -11.2 + (2/3)*log10(5);         % dyne-cm
    %k = -10.7340
    
    % Kanamori 1977 or 1978
    %Mw = (2/3) * log10(M0) + k;
    
    % Kanamori 1977 (modified form, but exact)
    A = 2 / (3*log(10));
    K = 0.2 * 10^16.8;
    Mw = A*log(M0/K);
    
elseif imag == 2
    % Harvard CMT
    % k = 2/3 * 16.1
    Mw = (2/3) * (log10(M0) - 16.1);
    
else
    error('imag must be 1 or 2.');
end

%----------------------------------

if 0==1
    M0 = 2.7 * 1e23;   % Cipar and Kanamori (1974), N-m
    Mw = m02mw(2,M0)
    Mw = m02mw(1,M0)
    
    % check direct calculation for M0 in N-m
    k = -(2/3)*(11.8 - log10(1/2e4) - 7)
    Mw = (2/3)*log10(M0) + k
end

%=====================================================
