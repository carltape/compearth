function xi = U2xi(U1,U2)
%U2xi compute rotation angle ([0,180]) between two bases
%
% INPUT
%   U1     3 x 3 x n set of rotation matrices
%   U2     3 x 3 x n set of rotation matrices
%
% OUTPUT   
%   xi     n x 1 set of rotation angles
%
% Note that these are the xi angles, NOT the xi_0 angles, which are the
% minimum rotation angles.
%
% See also U2q.m, which also provides xi.
% 
% Carl Tape, 8/12/2012
%

% check dimensions
[~,~,n1] = size(U1);
[~,~,n2] = size(U2);
if n1~=n2, error('U1 and U2 must be equal in dimension'); end
n = n1;

% compute U = U1' * U2
U = 0*U1;
for ii=1:n
    U(:,:,ii) = U1(:,:,ii)' * U2(:,:,ii);
end

% trace
trU = U(1,1,:) + U(2,2,:) + U(3,3,:);
trU = trU(:);

% TapeTape2013, Eq 2a
cosxi0 = (-1 + trU)/2;
xi0 = acos(cosxi0)*180/pi;

% correction for numerical errors
cosxi = cosxi0;
ipos = cosxi > 1;
ineg = cosxi < -1;
disp(sprintf('%i/%i dot products > 1',sum(ipos),n));
disp(sprintf('%i/%i dot products < 1',sum(ineg),n));
% this correction is needed for comparing U's that are very close to each other
cosxi(cosxi > 1) = 1;
cosxi(cosxi < -1) = -1;

xi = acos(cosxi)*180/pi;

if ~isreal(xi0)
   disp('IMAGINARY ENTRIES (U2xi.m):');
   for ii=1:n
      if ~imag(xi0(ii))==0
          disp(sprintf('%i/%i:',ii,n));
          disp('U1:'); U1(:,:,ii)
          disp('U2:'); U2(:,:,ii)
          disp('U = U1T*U2:');
          disp('tr(U), cos(xi), xi:'); trU(ii),cosxi0(ii),xi0(ii)
          disp('-----> cos(xi), xi:'); cosxi(ii),xi(ii)
      end
   end
end

%==========================================================================
