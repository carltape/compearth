function Uout = Udetcheck(Uin)
%UDETCHECK ensure that U is a rotation matrix with det(U) = 1
% 
% Carl Tape, 4/11/2011

[a,b,n] = size(Uin);
if any([a b]~=3), error('U must be 3 x 3 x n'); end

Uout = Uin;
for ii=1:n
   U0 = Uout(:,:,ii);
   if det(U0) < 0
       %disp(sprintf('det(U) < 0: flipping sign of 2nd column (%i/%i)',ii,n));
       Uout(:,:,ii) = [U0(:,1) -U0(:,2) U0(:,3)];
       %detU = det(Uout(:,:,ii)) 
       %UUt = Uout(:,:,ii) * Uout(:,:,ii)'
   end
end

% Uout = -Uin;
% for ii=1:n
%    U0 = Uout(:,:,ii);
%    if det(U0) < 0
%        Uout(:,:,ii) = [U0(:,1) U0(:,2) -U0(:,3)];
%    end
% end

% check
for ii=1:n
   U0 = Uout(:,:,ii);
   if det(U0) < 0
       detU = det(Uout(:,:,ii)) 
       UUt = Uout(:,:,ii) * Uout(:,:,ii)'
       error('det(U) < 0');
   end
end

% % eigenbasis convention (this seems to lead to the best agreement with
% % the CMT catalog output)
% U = -U;
% for ii=1:n
%    U0 = U(:,:,ii);
%    if det(U0) < 0
%        U(:,:,ii) = [U0(:,1) U0(:,2) -U0(:,3)];
%        %detU = det(U(:,:,ii)) 
%        %UUt = U(:,:,ii) * U(:,:,ii)'
%    end
% end

%==========================================================================