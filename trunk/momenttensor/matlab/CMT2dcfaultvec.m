function [Xout1,Xout2,Xout3,Xout4,Xout5] = CMT2dcfaultvec(Xin,itype,idisplay)
%CMT2DCFAULTVEC converts between fault vectors and double couple moment tensors
%
% A general moment tensor (non-double couple, isotropic) can be used as input.
% For moment tensors with strong CLVD component, the double-couple
% representation of fault planes has little physical meaning.
% See CMT2faultvec.m for the general case.
%
% itype = 1     convert moment tensors to two sets of fault vectors
%  INPUT:
%    Xin        6 x n moment tensors in CMT convention
%               M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%  OUTPUT:
%    Xout1      9 x n fault vectors for plane 1 [k1 ; d1 ; n1];
%    Xout2      9 x n fault vectors for plane 2 [k2 ; d2 ; n2];
%    Xout3      6 x n double couple moment tensors in CMT convention
%
% itype = 0     convert ONE SET of fault vectors to moment tensors
%  INPUT:
%    Xin        9 x n fault vectors [k1 ; d1 ; n1];
%  OUTPUT:
%    Xout1      6 x n double couple moment tensors in CMT convention
%    Xout2      9 x n eigenvectors [p1 ; p2 ; p3]
%
% calls
%   CMTdecom_iso.m
%   CMTdecom_dcclvdm
%   CMTdecom.m
%   transform_MT.m
%   Ueigvec.m
%   unit.m
%   swap.m
%   fvec2fmat.m, fmat2fvec.m
% called by CMT2dcfaultpar.m, dcfaultpar2CMT.m
%
% Carl Tape, 31-Mar-2011
%

%deg = 180/pi;
if nargin==2, idisplay=0; end

% convert M from up-south-east to south-east-up
% OBSOLETE: use convert_MT.m instead
%Bconvert = [0 0 1 ; 1 0 0 ; 0 1 0];

% DC in eigenbasis (see isort)
IDC = [1 0 -1 0 0 0]';

% OPTION 1: moment tensor to fault vectors
if itype==1
    M = Xin;

    % make sure M is 6 x n
    M = Mdim(M);
    
    % cartesian basis -- same as DT, p. 832
    ivec = [1 0 0]';    % south: x = r sin(th) cos(ph)
    jvec = [0 1 0]';    % east:  y = r sin(th) sin(ph)
    kvec = [0 0 1]';    % up:    z = r cos(th)
    upvec = kvec;       % up: used for distinguishing hanging wall from foot wall

    % catalog convention for rake for dip = 90
    % plane 1: strike from 135 to 315 (half circle)
    % plane 2: strike from 225 to 405 (half circle)
    n1vec = unit([-1 -1 0]'); % northwest (used for 90-deg-dip faults only)
    n2vec = unit([-1  1 0]'); % northeast (used for 90-deg-dip faults only)

    % remove isotropic part
    % note: M is in up-south-east
    [Miso,Mdev,trM] = CMTdecom_iso(M);

    % number of moment tensors
    n = length(trM);
    disp(sprintf('CMT2dcfaultvec.m: %i input moment tensors',n));

    % KEY: convert M from up-south-east to south-east-up
    %Mdev = transform_MT(Bconvert',Mdev);
    Mdev = convert_MT(1,5,Mdev);

    % compute eigenbasis for the DC
    % check 1: is the eigenbasis orthogonal (U*U^t = I)?
    % check 2: is the eigenbasis right-handed (det U = 1)?
    isort = 1;
    [eigvals,U] = CMTdecom(Mdev,isort);
    %M, Mdev
    U = Ueigvec(U);     % convention for rotation matrix (= eigenbasis)
    %U = Udetcheck(U);

    % obtain DC representation
    idecom = 1;         % 1 = default
    if idecom==1
        MDC = 0*M;              % initialize
        for ii=1:n
            U0 = U(:,:,ii);
            MDC(:,ii) = transform_MT(U0,IDC);
        end

    elseif idecom==2
        % note: here we don't need the CLVD, just the DC
        [MDC,MCLVD,f] = CMTdecom_dcclvd(Mdev);

    elseif idecom==3
        MDC = Mdev;     % testing only
    end

    % double couple
    Mmat = Mvec2Mmat(MDC,1);
    
    %-----------------------
    % compute fault vectors

    % initialize
    k1 = zeros(3,n); k2 = zeros(3,n);
    d1 = zeros(3,n); d2 = zeros(3,n);
    n1 = zeros(3,n); n2 = zeros(3,n);

    % two different codes (defaul = 1)
    iversion = 1;
    
    if iversion == 1

        % four sets of candidate normals
        % note: eigenbasis is assumed to contain orthonormal basis vectors
        Pp1Pp3 = squeeze(U(:,1,:) + U(:,3,:));
        Pp1Mp3 = squeeze(U(:,1,:) - U(:,3,:));
        Mp1Pp3 = squeeze(-U(:,1,:) + U(:,3,:));
        Mp1Mp3 = squeeze(-U(:,1,:) - U(:,3,:));

        % sort them in order of increasing theta (dip)
        th1 = xyz2tp(Pp1Pp3);
        th2 = xyz2tp(Pp1Mp3);
        th3 = xyz2tp(Mp1Pp3);
        th4 = xyz2tp(Mp1Mp3);
        thall = [th1 th2 th3 th4];
        [thallsort,isortmat] = sort(thall,2);
        th1 = thallsort(:,1);
        th2 = thallsort(:,2);
        th3 = thallsort(:,3);
        th4 = thallsort(:,4);
        %disp('sorted:'); for ii=1:n, disp(thall(ii,isortmat(ii,:))); end

        % define p2up
        p2up = squeeze(U(:,2,:));
        i1 = find(p2up(3,:) > 0);
        i2 = find(p2up(3,:) == 0);
        i3 = find(p2up(3,:) < 0);
        if length([i1 i2 i3]) ~= n, error('unexpected p2up values'); end
        % flip p2 vectors that point down
        if ~isempty(i3), p2up(:,i3) = -p2up(:,i3); end
        % flip horizontal p2 vectors that point to -y
        if ~isempty(i2)
            j3 = find(p2up(2,i2) < 0);
            if ~isempty(j3), p2up(:,i2(j3)) = -p2up(:,i2(j3)); end
        end

        % default assignment of normal vector (valid for UNEQUAL DIPS)
        for ii=1:n
            isort = isortmat(ii,:);     % KEY: ordering
            % normal for plane 1
            if isort(1)==1
                n1(:,ii) = Pp1Pp3(:,ii);
            elseif isort(1)==2
                n1(:,ii) = Pp1Mp3(:,ii);
            elseif isort(1)==3
                n1(:,ii) = Mp1Pp3(:,ii);
            elseif isort(1)==4
                n1(:,ii) = Mp1Mp3(:,ii);
            end
            % normal for plane 2
            if isort(2)==1
                n2(:,ii) = Pp1Pp3(:,ii);
            elseif isort(2)==2
                n2(:,ii) = Pp1Mp3(:,ii);
            elseif isort(2)==3
                n2(:,ii) = Mp1Pp3(:,ii);
            elseif isort(2)==4
                n2(:,ii) = Mp1Mp3(:,ii);
            end
        end

        % check for case where both planes have SAME DIPS
        EPSVAL = 0;
        isamedip = find( abs(th1-th2) <= EPSVAL );
        %isamedip = find( round((th1-th2)*deg) == 0);

        if ~isempty(isamedip)
            disp(sprintf('CMT2dcfaultvec.m: %i/%i MTs having planes with same dip',length(isamedip),n));
            for jj=1:length(isamedip)   % loop over 
                ii = isamedip(jj);

                if abs(th1(ii) - pi/2) > EPSVAL     % non-vertical planes
                    if det([ n1(:,ii) Mmat(:,:,ii)*n1(:,ii) p2up(:,ii)]) > 0
                        n0 = n1(:,ii);
                        n1(:,ii) = n2(:,ii);
                        n2(:,ii) = n0;
                    end

                else                                % vertical planes
                    ncall = [Pp1Pp3(:,ii) Pp1Mp3(:,ii) Mp1Pp3(:,ii) Mp1Mp3(:,ii)];
                    signD = zeros(4,1);
                    for kk=1:4
                        nc = ncall(:,kk);
                        signD(kk) = sign(det([ nc Mmat(:,:,ii)*nc p2up(:,ii)]));
                    end
                    %signD
                    kpos = find(signD == 1);
                    kneg = find(signD == -1);
                    nc1 = ncall(:,kneg(1));    % candidate for normal 1
                    nc2 = ncall(:,kpos(1));    % candidate for normal 2
                    n1(:,ii) = pickvertfault(nc1,n1vec);
                    n2(:,ii) = pickvertfault(nc2,n2vec);
                end
            end
        end

        %n1, n2, thallsort
    
    else

        % fault normal
        for ii=1:n
            % DC eigenbasis
            U0 = U(:,:,ii);

            % two potential normal vectors: sum or difference
            % note: not unit vectors
            Se1e3 = unit(U0(:,1) + U0(:,3));  % sum of eigenvectors 
            De1e3 = unit(U0(:,1) - U0(:,3));  % difference of eigenvectors

            % fault normal 1
            dot1 = dot(Se1e3,upvec);
            if dot1 > 0
                n1(:,ii) = Se1e3;
            elseif dot1 < 0
                n1(:,ii) = -Se1e3;
            else
                dot1b = dot(Se1e3,n1vec);
                %disp(sprintf('%i/%i vertical fault (plane 1)',ii,n));
                %Se1e3, n1vec, dot1b
                if dot1b > 0
                    n1(:,ii) = Se1e3;
                elseif dot1b < 0
                    n1(:,ii) = -Se1e3;
                else
                    error('very special case');
                end
            end
            % fault normal 2
            dot2 = dot(De1e3,upvec);
            if dot2 > 0
                n2(:,ii) = De1e3;
            elseif dot2 < 0
                n2(:,ii) = -De1e3;
            else
                dot2b = dot(De1e3,n2vec);
                %disp(sprintf('%i/%i vertical fault (plane 2)',ii,n));
                %De1e3, n2vec, dot2b
                if dot2b > 0
                    n2(:,ii) = De1e3;
                elseif dot2b < 0
                    n2(:,ii) = -De1e3;
                else
                    error('very special case');
                end
            end
        end

    end   % iversion
    
    % slip vector and strike vector
    for ii=1:n
        % slip vector
        M0 = Mmat(:,:,ii);
        d1(:,ii) = M0*n1(:,ii);
        d2(:,ii) = M0*n2(:,ii);

        % strike vector
        k1(:,ii) = cross(upvec,n1(:,ii));
        k2(:,ii) = cross(upvec,n2(:,ii));
    end

    % convert all to unit vectors
    k1 = unit(k1); k2 = unit(k2);
    d1 = unit(d1); d2 = unit(d2);
    n1 = unit(n1); n2 = unit(n2);

%     % output fault vectors: plane 1
%     Xout1 = [k1 ; d1 ; n1];
%     [a,b] = size(Xout1);
%     if a~=9, error('Xout1 should be 9 x n'); end
%     if b~=n, error('Xout1 should be 9 x n'); end
% 
%     % output fault vectors: plane 2
%     Xout2 = [k2 ; d2 ; n2];
%     [a,b] = size(Xout2);
%     if a~=9, error('Xout2 should be 9 x n'); end
%     if b~=n, error('Xout2 should be 9 x n'); end
    
    Xout1 = fvec2fmat(k1,d1,n1);
    Xout2 = fvec2fmat(k2,d2,n2);
   
    % convert M from south-east-up to up-south-east
    MDC0 = MDC;                        % only used for output below
    %MDC = transform_MT(Bconvert,MDC);  % note Bconvert, not Bconvert' (as above)
    MDC = convert_MT(5,1,MDC);
    Xout3 = MDC;
    Xout4 = U;
    Xout5 = eigvals;
    
    if idisplay==1
        for ii = 1:n
            disp('-------------------------------------');
            %disp(sprintf('CMT (up-south-east): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',M(:,ii)'));
            disp(sprintf('Mdev (south-east-up): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',Mdev(:,ii)'));
            disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',Mdev(:,ii)'/max(abs(Mdev(:,ii)))));
            disp(sprintf('MDC  (south-east-up): %11.3e%11.3e%11.3e%11.3e%11.3e%11.3e',MDC0(:,ii)'));
            disp(sprintf('                    : %11.3f%11.3f%11.3f%11.3f%11.3f%11.3f',MDC0(:,ii)'/max(abs(MDC0(:,ii)))));

            disp('eigenvectors (U = [v1 v2 v3]):');
            Ux = U(:,:,ii)
            disp('eigenvalues:');
            disp(sprintf('%11.3e%11.3e%11.3e',eigvals(:,ii)));
            UUt = Ux*Ux'
            detU = det(Ux)
            %disp('index, strike, dip, rake:');
            %disp(sprintf('%6i(1)%6.1f%6.1f%6.1f',ii,kap1(ii),theta1(ii),sig1(ii)));
            %disp(sprintf('%6i(2)%6.1f%6.1f%6.1f',ii,kap2(ii),theta2(ii),sig2(ii)));
            disp('fault vectors (and magnitudes):');
            disp(sprintf('   n1: %8.4f%8.4f%8.4f : %8.4e',n1(:,ii),sqrt(n1(1,ii)^2+n1(2,ii)^2+n1(3,ii)^2)));
            disp(sprintf('   k1: %8.4f%8.4f%8.4f : %8.4e',k1(:,ii),sqrt(k1(1,ii)^2+k1(2,ii)^2+k1(3,ii)^2) ));
            disp(sprintf('   d1: %8.4f%8.4f%8.4f : %8.4e',d1(:,ii),sqrt(d1(1,ii)^2+d1(2,ii)^2+d1(3,ii)^2)));
            disp(sprintf('   n2: %8.4f%8.4f%8.4f : %8.4e',n2(:,ii),sqrt(n2(1,ii)^2+n2(2,ii)^2+n2(3,ii)^2)));
            disp(sprintf('   k2: %8.4f%8.4f%8.4f : %8.4e',k2(:,ii),sqrt(k2(1,ii)^2+k2(2,ii)^2+k2(3,ii)^2)));
            disp(sprintf('   d2: %8.4f%8.4f%8.4f : %8.4e',d2(:,ii),sqrt(d2(1,ii)^2+d2(2,ii)^2+d2(3,ii)^2)));
        end
    end

%-------------------------------------------------------
% OPTION 2: fault vectors to moment tensor

elseif itype==0
    % get fault vectors from input matrix
    % note: strike vector is not needed
    [~,dvec,nvec] = fmat2fvec(Xin);
    [~,n] = size(dvec);
%     [a,n] = size(Xin);
%     if a~=9, error('Xin must be 9 x n: 3 sets of 3-vectors k, d, n'); end
%     %kvec = Xin(1:3,:);
%     dvec = Xin(4:6,:);
%     nvec = Xin(7:9,:);
    
    % eigenvectors (initial convention -- see Ueigvec.m below)
    p1vec = unit(dvec + nvec);
    p3vec = unit(dvec - nvec);
    p2vec = -cross(p1vec,p3vec);    % unit vector by definition
    U = zeros(3,3,n);
    U(:,1,:) = p1vec;
    U(:,2,:) = p2vec;
    U(:,3,:) = p3vec;
    %U = Ueigvec(U);
    U = Udetcheck(U);
    p1vec = U(:,1,:);
    p2vec = U(:,2,:);
    p3vec = U(:,3,:);
    
    % transform the double couple IDC
    % NOTE: To obtain the original moment tensor instead, transform diag(eigvals).
    MDC = zeros(6,n);
    for ii=1:n
        U0 = U(:,:,ii);
        MDC(:,ii) = transform_MT(U0,IDC);
    end
    
%     % "standard" computation (no physical meaning)
%     % note: this assumes "unit" seismic moment; use Xout1 = M0 * Xout1
%     %       to assign the desired moment
%     Xout1 = zeros(6,n);
%     d1 = dvec(1,:); d2 = dvec(2,:); dvec(3,:);
%     n1 = nvec(1,:); n2 = nvec(2,:); nvec(3,:);
% 	Xout1(1,:) = 2*d1.*n1;          % M_11
% 	Xout1(2,:) = 2*d2.*n2;          % M_22
% 	Xout1(3,:) = 2*d3.*n3;          % M_33
% 	Xout1(4,:) = d1.*n2 + d2.*n1;   % M_12 = M_21
% 	Xout1(5,:) = d1.*n3 + d3.*n1;   % M_13 = M_31
% 	Xout1(6,:) = d2.*n3 + d3.*n2;   % M_23 = M_32
    
    % convert M from south-east-up to up-south-east
    MDC = convert_MT(5,1,MDC);
    % note Bconvert, not Bconvert' as above
    %MDC = transform_MT(Bconvert,MDC);  % double couple
    %Mdev = transform_MT(Bconvert,Mdev);  % deviatoric
    
    % output matrices
    Xout1 = MDC;                            % DC moment tensors
    Xout2 = fvec2fmat(p1vec,p2vec,p3vec);   % eigenvectors
    
%     Xout2 = [p1vec ; p2vec ; p3vec];
%     [a,b] = size(Xout2);
%     if a~=9, error('Xout2 should be 9 x n'); end
%     if b~=n, error('Xout2 should be 9 x n'); end
    
else
    error(sprintf('itype (%i) must be 0 or 1'),itype);
    
end

%==========================================================================

function n = pickvertfault(nvec,evec)

vdot = dot(nvec,evec);
%disp(sprintf('%i/%i vertical fault (plane 1)',ii,n));
%nvec, evec, vdot
if vdot > 0
    n = nvec;
elseif vdot < 0
    n = -nvec;
else
    error('very special case');
end

%==========================================================================
