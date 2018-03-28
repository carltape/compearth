function M0 = mw2m0(imag,Mw)
%MW2M0 convert Mw to M0 using published formulas
%
% This function inputs a column vector of moment magnitudes (Mw)
% and outputs a column vector of moments (M0) in units of N-m.
%
% Carl Tape, 2003-09-26
%

if imag==1
    % Kanamori 1977
    M0 = 10.^( (3/2)*Mw + (11.8 - log10(5e-5)) );
    
elseif imag == 2
    % Harvard CMT
    M0 = 10.^( (3/2)*Mw + 16.1 );
    
else
    error('imag must be 1 or 2.');
end

% formulas are designed for units of dyne-cm, not N-m
M0 = M0 / 1e7;
