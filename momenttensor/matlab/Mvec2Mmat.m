function Mout = Mvec2Mmat(Min,itype)
%MVEC2MMAT converts between two representation of a set of moment tensors:
%    (1) 6 x n
%    (2) 3 x 3 x n
%
% moment tensor M = [M11 M22 M33 M12 M13 M23]
% The variable names below are based on 1=r, 2=t, 3=p.
%
% Carl Tape, 10-Nov-2010
%

if itype==1     % 6 x n --> 3 x 3 x n
    % make sure M is 6 x n
    [Min,n] = Mdim(Min);

    Mrr = Min(1,:); Mtt = Min(2,:); Mpp = Min(3,:);
    Mrt = Min(4,:); Mrp = Min(5,:); Mtp = Min(6,:);

    Mout = zeros(3,3,n);
    for ii = 1:n
        Mout(:,:,ii) = [ Mrr(ii) Mrt(ii) Mrp(ii) ;
                        Mrt(ii) Mtt(ii) Mtp(ii) ;
                        Mrp(ii) Mtp(ii) Mpp(ii) ];
    end

else            % 3 x 3 x n --> 6 x n
    if length(Min(:)) == 9
        % assume Min is symmetric and take the six upper triangular entries
        Mout = [Min(1,1) Min(2,2) Min(3,3) Min(1,2) Min(1,3) Min(2,3)]';
    else
        [~,~,n] = size(Min);
        Mout = zeros(6,n);
    
        % this simply takes the six upper triangular entries, without
        % checking if the 3 x 3 matrix is symmetric or not
        Mout(1,:) = squeeze(Min(1,1,:))';
        Mout(2,:) = squeeze(Min(2,2,:))';
        Mout(3,:) = squeeze(Min(3,3,:))';
        Mout(4,:) = squeeze(Min(1,2,:))';
        Mout(5,:) = squeeze(Min(1,3,:))';
        Mout(6,:) = squeeze(Min(2,3,:))';
    end
end

%--------------------------------------------------------------------------

if 0==1
    n = 4; M = randn(6,n); Mmat = Mvec2Mmat(M,1); Mcheck = Mvec2Mmat(Mmat,0); M, Mcheck
    n = 1; M = randn(6,n); Mmat = Mvec2Mmat(M,1); Mcheck = Mvec2Mmat(Mmat,0); M, Mcheck
end

%==========================================================================