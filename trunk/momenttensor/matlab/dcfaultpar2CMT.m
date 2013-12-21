function [MDC,k1,d1,n1,p1,p2,p3] = dcfaultpar2CMT(kap1,theta1,sig1,idisplay)
%DCFAULTPAR2CMT converts fault parameters (strike, dip, rake) to double couple moment tensor
%
% This is only meaningful for double couple moment tensors, for which the
% slip vector is in the fault plane (and therefore the rake angle is meaningful).
%
% INPUT:
%   kap1        strike (0 to 360)
%   theta1      dip (0 to 90)
%   sig1        rake (-180 to 180)
%
% OUTPUT:
%   MDC         6 x n "best" double couple moment tensor in CMT convention
%                   NOTE: M0 = 1 (and matrix norm is sqrt(2))
%   k1          strike vector (3 x n) in south-east-up convention
%   d1          slip vector (3 x n) in south-east-up convention
%   n1          normal vector (3 x n) in south-east-up convention
%   p1,p2,p3    eigenvectors (each 3 x n) in south-east-up convention
%   
% See inverse program CMT2dcfaultpar.m
%
% calls dcfaultvec2faultpar.m, CMT2dcfaultvec.m
%
% Carl Tape, 31-Mar-2011
%

%deg = 180/pi;
if nargin==3, idisplay=0; end
n = length(kap1);

% convert fault parameters to fault vectors
F1 = [kap1(:) theta1(:) sig1(:)];
V1 = dcfaultvec2faultpar(F1,0,idisplay);

% convert fault vectors to moment tensor
[MDC,U] = CMT2dcfaultvec(V1,0,idisplay);

% fault vectors in south-east-up convention
k1 = V1(1:3,:);
d1 = V1(4:6,:);
n1 = V1(7:9,:);

% eigenvectors in south-east-up convention
p1 = U(1:3,:);
p2 = U(4:6,:);
p3 = U(7:9,:);

if idisplay==1
    for ii = 1:n
        disp('-------------------------------------');
        %disp(sprintf('CMT (up-south-east): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',M(:,ii)'));
        %disp(sprintf('Mdev (south-east-up): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',Mdev(:,ii)'));
        %disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',Mdev(:,ii)'/max(abs(Mdev(:,ii)))));
        disp(sprintf('MDC  (south-east-up): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',MDC(:,ii)'));
        disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',MDC(:,ii)'/max(abs(MDC(:,ii)))));
        
        disp('eigenvectors (U = [v1 v2 v3]):');
        U = [p1(:,ii) p2(:,ii) p3(:,ii)]
        %disp('eigenvalues:');
        %disp(sprintf('%11.3e%11.3e%11.3e',eigvals(:,ii)));
        UUt = U*U'
        detU = det(U)
        disp('index, strike, dip, rake:');
        disp(sprintf('%6i(1)%6.1f%6.1f%6.1f',ii,kap1(ii),theta1(ii),sig1(ii)));
        %disp(sprintf('%6i(2)%6.1f%6.1f%6.1f',ii,kap2(ii),del2(ii),lam2(ii)));
        disp('fault vectors (and magnitudes):');
        disp(sprintf('   n1: %8.4f%8.4f%8.4f : %8.4e',n1(:,ii),sqrt(n1(1,ii)^2+n1(2,ii)^2+n1(3,ii)^2)));
        disp(sprintf('   k1: %8.4f%8.4f%8.4f : %8.4e',k1(:,ii),sqrt(k1(1,ii)^2+k1(2,ii)^2+k1(3,ii)^2) ));
        disp(sprintf('   d1: %8.4f%8.4f%8.4f : %8.4e',d1(:,ii),sqrt(d1(1,ii)^2+d1(2,ii)^2+d1(3,ii)^2)));
        %disp(sprintf('   n2: %8.4f%8.4f%8.4f : %8.4e',n2(:,ii),sqrt(n2(1,ii)^2+n2(2,ii)^2+n2(3,ii)^2)));
        %disp(sprintf('   k2: %8.4f%8.4f%8.4f : %8.4e',k2(:,ii),sqrt(k2(1,ii)^2+k2(2,ii)^2+k2(3,ii)^2)));
        %disp(sprintf('   d2: %8.4f%8.4f%8.4f : %8.4e',d2(:,ii),sqrt(d2(1,ii)^2+d2(2,ii)^2+d2(3,ii)^2)));
    end
end

%=======================================================================
