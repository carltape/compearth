function [MDC,kap1,theta1,sig1,kap2,theta2,sig2,k1,d1,n1,k2,d2,n2,U,lams] = CMT2faultpar(M,idisplay)
%
% This converts a seismic moment tensor into two sets of fault parameters
% for its double couple.
%
% For moment tensors with strong CLVD component, the double-couple
% representaion of fault parameters has little physical meaning.
%
% INPUT
%   M       6 x n moment tensors in up-south-east (CMT) convention
%           M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
% OUTPUT
%   MDC             6 x n closest double couple moment tensor
%   kap1,kap2       strike (0 to 360)
%   theta1,theta2   dip (0 to 90)
%   sig1,sig2       rake (-180 to 180)
%   k1,k2           strike vector (3 x n) in south-east-up convention
%   d1,d2           slip vector (3 x n) in south-east-up convention
%   n1,n2           normal vector (3 x n) in south-east-up convention
%
% See inverse program faultpar2CMT.m.
%
% calls CMT2faultvec.m, faultvec2faultpar.m, swap.m
% called by test_CMT2faultpar.m
%

deg = 180/pi;
if nargin==1, idisplay=0; end

% moment tensor to fault vectors
% note: we could return the eigenbasis U here also (U = [p1 p2 p3])
[V1,V2,MDC,U,lams] = CMT2faultvec(M,1,idisplay);

% fault vectors to fault parameters
F1 = faultvec2faultpar(V1,1,idisplay);  % plane 1
F2 = faultvec2faultpar(V2,1,idisplay);  % plane 2

% fault vectors
k1 = V1(1:3,:); d1 = V1(4:6,:); n1 = V1(7:9,:);
k2 = V2(1:3,:); d2 = V2(4:6,:); n2 = V2(7:9,:);

% fault parameters (strike, dip, rake)
kap1 = F1(:,1); theta1 = F1(:,2); sig1 = F1(:,3);
kap2 = F2(:,1); theta2 = F2(:,2); sig2 = F2(:,3);

%-----------------------
% to match CMT output, the first plane is taken to be the SHALLOW dip
% note: is there a better way to swap elements of vectors in matlab?

% % find dips of plane 1 that are greater than dips of plane 2
% iswap = find(theta1 > theta2);
% 
% D1 = [kap1 theta1 sig1 k1' d1' n1'];
% D2 = [kap2 theta2 sig2 k2' d2' n2'];
% [D1,D2] = swap(D1,D2,iswap);
% 
% kap1 = D1(:,1); theta1 = D1(:,2); sig1 = D1(:,3);
% k1 = D1(:,4:6)';
% d1 = D1(:,7:9)';
% n1 = D1(:,10:12)';
% 
% kap2 = D2(:,1); theta2 = D2(:,2); sig2 = D2(:,3);
% k2 = D2(:,4:6)';
% d2 = D2(:,7:9)';
% n2 = D2(:,10:12)';

if idisplay==1
    n = length(kap1);
    for ii = 1:n
        disp('-------------------------------------');
        disp(sprintf('CMT  (up-south-east): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',M(:,ii)'));
        %disp(sprintf('Mdev (south-east-up): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',Mdev(:,ii)'));
        %disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',Mdev(:,ii)'/max(abs(Mdev(:,ii)))));
        disp(sprintf('MDC  (up-south-east): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',MDC(:,ii)'));
        %disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',MDC(:,ii)'/max(abs(MDC(:,ii)))));
        
        %disp('eigenvectors (U = [v1 v2 v3]):');
        %U = U(:,:,ii)
        %disp('eigenvalues:');
        %disp(sprintf('%11.3e%11.3e%11.3e',eigvals(:,ii)));
        %UUt = U*U'
        %detU = det(U)
        disp('index, strike, dip, rake:');
        disp(sprintf('%6i(1)%6.1f%6.1f%6.1f',ii,kap1(ii),theta1(ii),sig1(ii)));
        disp(sprintf('%6i(2)%6.1f%6.1f%6.1f',ii,kap2(ii),theta2(ii),sig2(ii)));
        disp('fault vectors (and magnitudes):');
        disp(sprintf('   n1: %8.4f%8.4f%8.4f : %8.4e',n1(:,ii),sqrt(n1(1,ii)^2+n1(2,ii)^2+n1(3,ii)^2)));
        disp(sprintf('   k1: %8.4f%8.4f%8.4f : %8.4e',k1(:,ii),sqrt(k1(1,ii)^2+k1(2,ii)^2+k1(3,ii)^2) ));
        disp(sprintf('   d1: %8.4f%8.4f%8.4f : %8.4e',d1(:,ii),sqrt(d1(1,ii)^2+d1(2,ii)^2+d1(3,ii)^2)));
        disp(sprintf('   n2: %8.4f%8.4f%8.4f : %8.4e',n2(:,ii),sqrt(n2(1,ii)^2+n2(2,ii)^2+n2(3,ii)^2)));
        disp(sprintf('   k2: %8.4f%8.4f%8.4f : %8.4e',k2(:,ii),sqrt(k2(1,ii)^2+k2(2,ii)^2+k2(3,ii)^2)));
        disp(sprintf('   d2: %8.4f%8.4f%8.4f : %8.4e',d2(:,ii),sqrt(d2(1,ii)^2+d2(2,ii)^2+d2(3,ii)^2)));
    end
end

%==========================================================================
