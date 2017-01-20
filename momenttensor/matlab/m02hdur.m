function hdur = m02hdur(M0)
%M02HDUR input seismic moment in N-m and output half-duration in seconds
%
% Equation: empirical formula stated on p. 178 of Dahlen and Tromp (1999), p. 178

M0 = M0(:);
hdur = 2.4e-6 * M0.^(1/3);

%--------------------------------------------------------------------------

if 0==1
    Mw = linspace(2,5,100);
    M0 = mw2m0(1,Mw);
    hdur = m02hdur(M0);
    
    figure; nr=2; nc=1;
    
    subplot(nr,nc,1); plot(Mw,hdur,'b'); grid on;
    xlabel('Moment magnitude, Mw');
    ylabel('Empirical half-duration, hdur, s');
    
    % SPECFEM2D: f0 = 1/hdur
    subplot(nr,nc,2); plot(Mw,1./hdur,'b'); grid on;
    xlabel('Moment magnitude, Mw');
    ylabel('Empirical central frequency, f0, Hz');
end

%==========================================================================