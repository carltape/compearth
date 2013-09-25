function xi = U2xi(U)
%U2xi compute rotation angle ([0,180]) between two bases
%
% INPUT
%   U      3 x 3 x n set of rotation matrices
%
% OUTPUT   
%   xi     n x 1 set of rotation angles
%
% Note that this is the xi angle, NOT the xi0 angle, which is the
% minimum rotation angle. See U2xi0.m to get both.
% 
% Carl Tape, 8/12/2012
%

[~,~,n] = size(U);

% trace
trU = U(1,1,:) + U(2,2,:) + U(3,3,:);
trU = trU(:);

% TapeTape2013, Eq 2a
cosxi = (-1 + trU)/2;
xiA = acos(cosxi)*180/pi;

% correction for possible numerical errors
ipos = cosxi > 1;
ineg = cosxi < -1;
disp(sprintf('%i/%i dot products > 1',sum(ipos),n));
disp(sprintf('%i/%i dot products < 1',sum(ineg),n));
% this correction is needed for comparing U's that are very close to each other
cosxi(cosxi > 1) = 1;
cosxi(cosxi < -1) = -1;
xi = acos(cosxi)*180/pi;

if ~isreal(xiA)
    cosxiA = (-1 + trU)/2;
   disp('IMAGINARY ENTRIES (U2xi.m):');
   for ii=1:n
      if ~imag(xiA(ii))==0
          disp(sprintf('%i/%i:',ii,n));
          disp('U:'); U(:,:,ii)
          disp('tr(U), cos(xi), xi:'); trU(ii),cosxiA(ii),xiA(ii)
          disp('-----> cos(xi), xi:'); cosxi(ii),xi(ii)
      end
   end
end

%==========================================================================
