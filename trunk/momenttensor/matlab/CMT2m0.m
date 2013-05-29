function M0 = CMT2m0(im0,M)
%CMT2M0 convert from moment tensor to scalar seismic moment
%
% INPUT
%   im0     =1 for Silver and Jordan (1982) formula
%           =2 for GCMT formula
%           =3 for old 'Caltech way'
%   M       6 x n set of moment tensors, M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%
% OUTPUT
%   M0      1 x n vector of seismic moments
%
% The units of M0 are the same as the elements of Mij, which should be
% Newton-meter (N-m), although it does not affect the function here.
%
% See Ekstrom email (11-Oct-2006) and corresponding Latex notes.
%
% calls Mdim.m
%
% Carl Tape, 02-Feb-2007
%

% make sure M is 6 x n
[M,n] = Mdim(M);

Mrr = M(1,:); Mtt = M(2,:); Mpp = M(3,:);
Mrt = M(4,:); Mrp = M(5,:); Mtp = M(6,:);

M0 = zeros(1,n);   % NOTE: row vector

for ii=1:n
    
    % if you need to compute the eigenvalues
    if or(im0==2,im0==3)
        % convention: r (up), theta (south), phi (east)
        Mcmt = zeros(3,3);
        Mcmt = [ Mrr(ii) Mrt(ii) Mrp(ii) ;
                 Mrt(ii) Mtt(ii) Mtp(ii) ;
                 Mrp(ii) Mtp(ii) Mpp(ii) ];

        [V, D] = eig(Mcmt);
        lams = diag(D)';
        isign = sign(lams);
        %(inv(V)*Mcmt*V - D) / norm(Mcmt)

        % adjust the order of eigenvectors
        % [lamsort, isort] = sort(abs(lams),'descend');
        % Vsort = V(:,isort)
        % Dsort = diag(lamsort.*isign(isort))
        % (inv(Vsort)*Mcmt*Vsort - Dsort) / norm(Mcmt)  % check
    end

    % formula to convert M to M0
    switch im0
        case 1
            % (1) Silver and Jordan (1982)
            M0(ii) = 1/sqrt(2) * sqrt( Mrr(ii)^2 + Mtt(ii)^2 + Mpp(ii)^2 ...
                + 2*(Mrt(ii)*Mrt(ii) + Mrp(ii)*Mrp(ii) + Mtp(ii)*Mtp(ii) ) );
        case 2
            % (2) GCMT -- double-couple moment
            M0(ii) = (abs(min(lams)) + abs(max(lams)) ) / 2;
        case 3
            % (3) 'Caltech way'
            M0(ii) = max(abs(lams));
    end
end

%==========================================================================
