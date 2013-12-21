function Xout = dcfaultvec2faultpar(Xin,itype,idisplay)
%DCFAULTVEC2FAULTPAR converts between fault vectors and strike/dip/rake
%
% This is only meaningful for double couple moment tensors, for which the
% slip vector is in the fault plane (and therefore the rake angle is meaningful).
% The basis for all vectors is SOUTH-EAST-UP.
%
% itype = 1     convert fault vectors to fault parameters
%  INPUT:
%    Xin        9 x n fault vectors [k1 ; d1 ; n1]
%  OUTPUT:
%    Xout       n x 3 set of fault parameters [kap1 theta1 sig1]
%
% itype = 0     convert fault parameters to fault vectors
%  INPUT:
%    Xin        n x 3 set of fault parameters [kap1 theta1 sig1]
%  OUTPUT:
%    Xout       9 x n fault vectors [k1 ; d1 ; n1]
%
% Carl Tape, 31-Mar-2011
%

deg = 180/pi;
if nargin==2, idisplay=0; end

% OPTION 1: fault vectors (k, n, d) to fault parameters (strike, dip, rake)
if itype==1
    
    % get input fault vectors
    [a,n] = size(Xin);
    if a~=9, error('Xout should be 9 x n'); end
    disp(sprintf('dcfaultvec2faultpar.m: %i sets of input fault vectors',n));
    k1 = Xin(1:3,:);
    d1 = Xin(4:6,:);
    n1 = Xin(7:9,:);
    
    % initialize
    kap1 = zeros(n,1);
    theta1 = zeros(n,1);
    sig1 = zeros(n,1);

    % dip and strike
    % dip: theta is [0, 90] because we "force" the normal to point up
    % strike: kappa is unrestricted, but we force it to [0, 360]
    % NOTE: instead of xyz2tp.m, we should use rotmat.m (see below)
    [n1th,n1ph] = xyz2tp(n1);
    theta1 = n1th;
    kap1 = wrapTo2Pi(pi/2 - n1ph);

    % rake
    for ii=1:n
        % sign will return -1, 0, or 1
        % e.g., =0 if strike vector is opposite from slip vector, in which case
        %    we can let acos below return 180
        sig1 = sign(det([k1(:,ii) d1(:,ii) n1(:,ii)]));
        if sig1==0, sig1 = 1; end
        % compute rake
        dot1 = dot(k1(:,ii),d1(:,ii));
        % ugly correction to assure that rake is not imaginary
        if abs(dot1) > 1, dot1 = sign(dot1); end
        % note: acos operation to get sigma assumes UNIT vectors k and d
        sig1(ii) = sig1 * acos( dot1 );
    end

    % convert to degrees
    kap1   = kap1*deg;
    theta1 = theta1*deg;
    sig1   = sig1*deg;
    
    % output matrix is n x 3
    Xout = [kap1(:) theta1(:) sig1(:)];

%-------------------------------------------------------
% OPTION 2: fault parameters (strike, dip, rake) to fault vectors (k, n, d)

elseif itype==0
    % NOTE: xyz-coordinates are such that the coordinate vectors i,j,k point
    %       SOUTH, EAST, UP (see CMT2dcfaultpar.m)
    upvec = [0 0 1]';
    northvec = [-1 0 0]';
    
    % get input fault parameters
    [n,b] = size(Xin);
    if b~=3, error('Xin should be n x 3'); end
    disp(sprintf('dcfaultvec2faultpar.m: %i sets of input fault parameters',n));
    kap1 = Xin(:,1);
    theta1 = Xin(:,2);
    sig1 = Xin(:,3);
    
    % compute fault vectors
    k1 = zeros(3,n);
    d1 = zeros(3,n);
    n1 = zeros(3,n);
    for ii=1:n
        if 1==1
            % most efficient operations (uses rotmat.m)
            
            % fault normal vector (from strike and dip angles)
            R = rotmat(-kap1(ii),3) * rotmat(-theta1(ii),1);
            n1(:,ii) = R*upvec;

            % slip vector
            R = rotmat(-kap1(ii),3) * rotmat(-theta1(ii),1) * rotmat(sig1(ii),3);
            d1(:,ii) = R*northvec;

            % strike vector
            k1(:,ii) = cross(upvec,n1(:,ii));
            
        else
            % most conceptual operations (uses rotmat_gen.m):
            % To get the strike vector, rotate the north vector through angle -kap about the up vector.
            % To get the normal vector, rotate the up vector through angle theta about the strike vector.
            % To get the slip vector, rotate the strike vector through angle sigma about the normal vector.
            
            % strike vector
            R = rotmat(-kap1(ii),3);
            k1(:,ii) = R*northvec;
            
            % fault normal vector
            R = rotmat_gen(k1(:,ii),theta1(ii));
            n1(:,ii) = R*upvec;
            
            % slip vector
            R = rotmat_gen(n1(:,ii),sig1(ii));
            d1(:,ii) = R*k1(:,ii);
        end
    end
    
    % ensure that all are unit vectors
    k1 = unit(k1);
    d1 = unit(d1);
    n1 = unit(n1);
    
    % output 
    Xout = [k1 ; d1 ; n1];
    [a,b] = size(Xout);
    if a~=9, error('Xout should be 9 x n'); end
    if b~=n, error('Xout should be 9 x n'); end
  
else
    error(sprintf('itype (%i) must be 0 or 1'),itype);   
    
end

if idisplay==1
    for ii=1:n
        disp('-------------------------------------');
        disp('index, strike, dip, rake:');
        disp(sprintf('%6i%6.1f%6.1f%6.1f',ii,kap1(ii),theta1(ii),sig1(ii)));
        disp('fault vectors (and magnitudes):');
        disp(sprintf('   k: %8.4f%8.4f%8.4f : %8.4e',k1(:,ii), sqrt(k1(1,ii)^2+k1(2,ii)^2+k1(3,ii)^2) ));
        disp(sprintf('   d: %8.4f%8.4f%8.4f : %8.4e',d1(:,ii), sqrt(d1(1,ii)^2+d1(2,ii)^2+d1(3,ii)^2) ));
        disp(sprintf('   n: %8.4f%8.4f%8.4f : %8.4e',n1(:,ii), sqrt(n1(1,ii)^2+n1(2,ii)^2+n1(3,ii)^2) ));
    end
end

%==========================================================================
% EXAMPLES

if 0==1
    kap1 = 350;
    theta1 = 80;
    sig1 = -20;
    Xin = [kap1 theta1 sig1]
    Xout = dcfaultvec2faultpar(Xin,0,1);
    Xcheck = dcfaultvec2faultpar(Xout,1,1);
end

%==========================================================================
