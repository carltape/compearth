function [nu,alpha,N1,N2,lam] = CMT2faultpar(M,idisplay)
%CMT2FAULTPAR convert (full) moment tensor to fault parameters
%
% Fault normal vectors N1 and N2 will have the same basis as the input
% tensor M.
%
% INPUT
%   M       6 x n set of moment tensors
%   idisplay
%
% OUTPUT
%   nu      Poisson parameter
%   alpha   angle between N1 and N2 (fault planes)
%   N1      unit normal vector, fault plane 1 (same basis as M)
%   N2      unit normal vector, fault plane 2 (same basis as M)
%   lam     eigenvalues of M
%
% Note: N1 and N2 can point in any direction. You can flip the sign of BOTH
%       N1 and N2 without changing the system.
%
% See Tape and Tape (2013), "The classical model for moment tensors",
%    including the sample calculation in Appendix A (see TT2013AppA.m).
%
% Carl Tape, 12/21/2013
%

if nargin==1, idisplay=0; end

% make sure M is 6 x n
[M,n] = Mdim(M);

% get nu and alpha from the eigenvalues
isort = 1;
[lam,U] = CMTdecom(M,isort);
[nu,alpha] = lam2nualpha(lam);

% Eq 36
% N1 = U Y_{-alpha/2} e_1
% N2 = U Y_{+alpha/2} e_1
N1 = NaN(3,n);
N2 = NaN(3,n);
for ii=1:n
   Ux = U(:,:,ii);
   Y1 = rotmat(-alpha(ii)/2,2);
   Y2 = rotmat( alpha(ii)/2,2);
   N1(:,ii) = Ux * Y1 * [1 0 0]';
   N2(:,ii) = Ux * Y2 * [1 0 0]';
%    if ~isreal(N1(:,ii))
%        ii, lam(:,ii), alpha(ii), Ux, Y1, N1(:,ii)
%        error('N1 is complex');
%    end
%    if ~isreal(N1(:,ii))
%        ii, lam(:,ii), alpha(ii), Ux, Y2, N2(:,ii)
%        error('N2 is complex');
%    end
end

if idisplay==1
   for ii=1:n
        disp('-------------------------------------');
        %disp(sprintf('CMT (up-south-east): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',M(:,ii)'));
        disp(sprintf('M : %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',M(:,ii)'));
        disp(sprintf('  : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',M(:,ii)'/max(abs(M(:,ii)))));
        disp('eigenvectors (U = [u1 u2 u3]):');
        
        Ux = U(:,:,ii)
        disp('eigenvalues:');
        disp(sprintf('%11.3e%11.3e%11.3e',lam(:,ii)));
        UUt = Ux*Ux'
        detU = det(Ux)

        disp('fault vectors (and magnitudes):');
        disp(sprintf('   N1: %8.4f%8.4f%8.4f : %8.4e',N1(:,ii),sqrt(N1(1,ii)^2+N1(2,ii)^2+N1(3,ii)^2)));
        disp(sprintf('   N2: %8.4f%8.4f%8.4f : %8.4e',N2(:,ii),sqrt(N2(1,ii)^2+N2(2,ii)^2+N2(3,ii)^2)));
        disp(sprintf('alpha = %.2f [check acos(n1.n2) = %.2f]',alpha(ii),180/pi*acos(dot(N1(:,ii),N2(:,ii)))));
        disp(sprintf('nu = %.2f',nu(ii)));
   end
end

%==========================================================================

if 0==1
    %% Tape and Tape (2013), Appendix A (and Tables A1 and A2)
    M = [  3.108304932835845
           3.044425632830430
           3.382269434333724
          -4.855033301709626
          -1.949280336439431
           1.110527600460120  ];
    [nu,alpha,N1,N2,lam] = CMT2faultpar(M,1);
    
    %% sample matrices to use to get a set of fault vectors
    n = 1000;  % 1 or 1000
    M = -1 + 2*rand(6,n);
    [nu,alpha,N1,N2,lam] = CMT2faultpar(M,1);
    % note: distribution of alpha is NOT flat
    figure; plot_histo(alpha,[0:5:180]);
    
    %% let us assume that M (and the fault vectors) are in the south-east-up basis
    % now get spherical coordinates
    % ph will be w.r.t. south in the positive sense,
    % so ph=0 points south (az=180) and ph=90 points east (az=90)
    % four possible fault normals
    [th1,ph1]   = xyz2tp(N1);
    [th1m,ph1m] = xyz2tp(-N1);
    [th2,ph2]   = xyz2tp(N2);
    [th2m,ph2m] = xyz2tp(-N2);
    thall = 180/pi*[th1 ; th1m ; th2 ; th2m];
    phall = 180/pi*[ph1 ; ph1m ; ph2 ; ph2m];
    % convert phi to azimuth, measured counter-clockwise from north
    azall = wrapTo360(-phall + 180);
    % convert to plunge
    plall = thall - 90;
    % check that this makes sense
    figure; plot(phall,azall,'.'); grid on; xlabel('phi'); ylabel('azimuth');
    % get xy points for plotting
    % separate vectors into two sets: those pointing up vs down
    iup = find(thall <= 90);
    idn = find(thall > 90);
    % these two plots should be antipodes of each other
    [Ux,Uy] = pa2xy(plall(iup),azall(iup)); xlabel('upper hemisphere');
    [Dx,Dy] = pa2xy(plall(idn),azall(idn)); xlabel('lower hemisphere');
end

%==========================================================================
