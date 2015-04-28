function omega = CMT2omega(M1,M2)
%CMT2OMEGA compute the omega angles between two moment tensors
%
% INPUT
%   M1,M2     6 x n moment tensors: M = [M11 M22 M33 M12 M13 M23]
%
% OUTPUT
%   omega     9-angle between moment tensors
%
% If M2 is missing, then a reference M is assumed for M1.
% If M1 or M2 is size 6x1, then all moment tensors in the other set will eb
% measured w.r.t. the 6x1 matrix.
%
% omega is a measure between two moment tensors and combines differences in
% eigenvalues (source types) and orientation, but not magnitude.
% omegadc, the angle introduced in TapeTape2012 ("Angle between principal
% axis triples"), measures a difference in orientation only.
% omegadc is the 9-angle (or omega angle) between the closest double couples (see CMT2omegadc_xi0.m).
% 
% Carl Tape 2/20/2015

bdisplay = false;
bfigure = true;

if nargin==1        % M1 only specified
    [M,n] = Mdim(M1);
    Mref = 1/sqrt(2)*[1 0 -1 0 0 0]';
    M2   = repmat(Mref,1,n);
    warning('measuring omega from the DC to all other MTs in the set');
    
else                % M1 and M2 specified
    [M1,n1] = Mdim(M1);
    [M2,n2] = Mdim(M2);
    if or(n1==1,n2==1)
        if n1~=n2
            warning('measuring omega from one MT to all other MTs in the set');
        end
        if n1==1    % M1 dimension 1
           n = n2;
           M1
           M1 = repmat(M1,1,n);
        else        % M2 dimension 1
           n = n1;
           M2
           M2 = repmat(M2,1,n);
        end
        
    elseif n1~=n2   % M1 and M2 have different dimension (but neither is dimension 1)
        n1, n2
        error('M1 or M2 must have only one MT if they are not the same size');
        
    else            % M1 and M2 have same dimension
        n = n1;
    end
end
  
% OMEGA
% convert to matrix
M1mat = Mvec2Mmat(M1,1);
M2mat = Mvec2Mmat(M2,1);
cosom = zeros(n,1);
for ii=1:n
   M1x = M1mat(:,:,ii);
   M2x = M2mat(:,:,ii);
   cosom(ii) = dot(M1x(:),M2x(:)) / (norm(M1x(:))*norm(M2x(:)));
end
% correction for possible numerical errors
% this correction is needed for comparing U's that are very close to each other
ipos = cosom > 1;
ineg = cosom < -1;
disp(sprintf('%i/%i dot products > 1',sum(ipos),n));
disp(sprintf('%i/%i dot products < 1',sum(ineg),n));
cosom(ipos) = 1;
cosom(ineg) = -1;
omega = acos(cosom) * 180/pi;

if bdisplay==1
    disp('CMT2omega.m all omega and xi0:');
   for ii=1:n
      disp(sprintf('%12i/%6i: omega = %8.4f',ii,n,omega(ii))); 
   end
end

if and(bfigure, n>1)
   OMAX = 0.07;  % will depend on bin size
   figure; hold on; plot_histo(omega,[0:5:180]);
   plot([90 90],[0 OMAX],'r','linewidth',2);
   set(gca,'xtick',0:10:180); %axis([0 120 0 OMAX]);
   xlabel('omega angle describing difference in orientation');
   title(sprintf('OMEGA: min = %.2f, max = %.2f',min(omega),max(omega)));
end

%==========================================================================
