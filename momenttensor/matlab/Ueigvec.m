function [Uout,igreen,ieps] = Ueigvec(Uin,EPSVAL)
%UEIGVEC converts input basis into a special convention
%
% A positive determinant ensures that the eigenbasis is equivalent to a
% rotation.
% 
% calls iUgreen.m, Udetcheck.m, rotmat.m, rotmat2rotvec.m
% called by CMT2faultvec.m
%
% Carl Tape, 11-March-2011
%


[a,b,n] = size(Uin);
if any([a b]~=3), error('U must be 3 x 3 x n'); end

% tolerance for determining whether U is a rotation matrix
if nargin==1, EPSVAL = 1e-6; end

% ensure that all matrices are ROTATION matrices
%Uin, disp('determinant check');
Uin = Udetcheck(Uin);

% initialize
igreen = NaN*ones(n,1);
Uout = Uin;

ieps = 0;   % counter for identity matrix assignments (EPSVAL)
for ii = 1:n
    U0 = Uin(:,:,ii);
    
    % loop over four possible eigenbases
    for kk=1:4
        switch kk
            case 1, U = U0;
            case 2, U = U0 * rotmat(180,1);
            case 3, U = U0 * rotmat(180,2);
            case 4, U = U0 * rotmat(180,3);
        end
        %disp('----------------------------'); disp(sprintf('case %i/4',kk)); U
        
        % get rotation vector
        % note: rotangle not needed here
        [rotaxis,umag] = rotmat2rotvec(U,EPSVAL);
        if umag==0, ieps=ieps+1; end
        
%         % for a rotation matrix there will be one eigenvector with
%         % eigenvalue =1, and two eigenvectors with complex eigenvalues
%         [V,D] = eig(U);
%         ipick = find( abs(diag(D) - 1) <= EPSVAL);
%         if length(ipick) == 1
%             v = V(:,ipick);     % rotation axis candidate
%         else
%             error('Ueigvec.m: must have exactly one eigenvector with eigenvalue =1');
%         end
%         
%         % wperp for candidate rotation axis
%         w = Wperp(v);
%         
%         % candidate rotation angle
%         rotangle_candidate = fangle_signed(w, U*w, v);
%         
%         % rotation axis
%         rotaxis  = sin( rotangle_candidate/2 / deg) * v;
%         
%         % rotation angle (not needed here)
%         rotangle = fangle(w, U*w);
        
        % check if vector is "in the green region" (latex notes)
        igreen0 = iUgreen(rotaxis);
        if ~isnan(igreen0)
            igreen(ii) = igreen0;
            Uout(:,:,ii) = U;
        end
    end
end

disp(sprintf('%i/%i input matrices are identity matrices',ieps,n));

%============================

% [a,b,n] = size(Uin);
% if any([a b]~=3), error('U must be 3 x 3 x n'); end
% 
% Uout = Uin;     % initialize
% 
% for ii = 1:n
%     % four possible eigenbases
%     U1 = Uin(:,:,ii);
%     U2 = U1 * rotmat(180,1);
%     U3 = U1 * rotmat(180,2);
%     U4 = U1 * rotmat(180,3);
% 
%     % compute euler vector
%     e1 = rotmat2euler(U1);
%     e2 = rotmat2euler(U2);
%     e3 = rotmat2euler(U3);
%     e4 = rotmat2euler(U4);
%     erad1 = e1(3)*pi/180;
%     erad2 = e2(3)*pi/180;
%     erad3 = e3(3)*pi/180;
%     erad4 = e4(3)*pi/180;
%     exyz1 = euler_convert([e1(1:2) ; sin(erad1/2) ],0);
%     exyz2 = euler_convert([e2(1:2) ; sin(erad2/2) ],0);
%     exyz3 = euler_convert([e3(1:2) ; sin(erad3/2) ],0);
%     exyz4 = euler_convert([e4(1:2) ; sin(erad4/2) ],0);
%     
%     % pick the rotation vector with wx > 0 and wz > 0
%     eall = [exyz1 exyz2 exyz3 exyz4];
%     i0 = find(and(eall(1,:) > 0, eall(3,:) > 0));
%     if ~isempty(i0)
%         if length(i0) > 1
%             %ii, i0, eall
%             %error('more than one U matrix satisfying criteria');
%             disp(sprintf('WARNING (Ueigvec:m): %i/4 U matrices satisfying criteria for index %i',length(i0),ii));
%             i0 = i0(1);     % just take the first one
%         else
%             switch i0
%                 case 1, Uout(:,:,ii) = U1;
%                 case 2, Uout(:,:,ii) = U2;
%                 case 3, Uout(:,:,ii) = U3;
%                 case 4, Uout(:,:,ii) = U4;
%             end
%         end
%     end
% end

%============================

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

%=====================================================