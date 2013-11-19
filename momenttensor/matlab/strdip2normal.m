function Xout = strdip2normal(Xin,itype,idisplay)
%STRDIP2NORMAL converts between fault normal vector and strike/dip
%
% The basis for all vectors is SOUTH-EAST-UP.
%
% itype = 0     convert fault vectors to fault parameters
%  INPUT:
%    Xin        3 x n fault normal vectors
%  OUTPUT:
%    Xout       n x 2 strikes and dips
%
% itype = 1     convert fault parameters to fault vectors
%  INPUT:
%    Xin        n x 2 strikes and dips
%  OUTPUT:
%    Xout       3 x n fault normal vectors
%
% Copied from faultvec2faultpar.m on 13-Aug-2013
%
% Carl Tape, 31-Mar-2011
%

deg = 180/pi;
if nargin==2, idisplay=0; end

% OPTION 1: fault normal vector to strike and dip
if itype==0
    
    % get input fault vectors
    [a,n] = size(Xin);
    if a~=3, error('Xout should be 3 x n'); end
    n1 = Xin;
    
    % dip and strike
    % dip: theta is [0, 90] because we "force" the normal to point up
    % strike: kappa is unrestricted, but we force it to [0, 360]
    % NOTE: instead of xyz2tp.m, we should use rotmat.m (see below)
    [n1th,n1ph] = xyz2tp(n1);
    theta1 = n1th;
    kap1 = wrapTo2Pi(pi/2 - n1ph);

    % convert to degrees
    kap1   = kap1*deg;
    theta1 = theta1*deg;
    
    % output matrix is n x 2
    Xout = [kap1(:) theta1(:)];

%-------------------------------------------------------
% OPTION 2: strike and dip to fault normal vector

elseif itype==1
    % NOTE: xyz-coordinates are such that the coordinate vectors i,j,k point
    %       SOUTH, EAST, UP (see CMT2faultpar.m)
    upvec = [0 0 1]';
    northvec = [-1 0 0]';
    
    % get input fault parameters
    [n,b] = size(Xin);
    if b~=2, error('Xin should be n x 2'); end
    disp(sprintf('faultvec2faultpar.m: %i sets of input fault parameters',n));
    kap1   = Xin(:,1);
    theta1 = Xin(:,2);
    
    % compute fault vectors
    n1 = zeros(3,n);
    for ii=1:n
        % fault normal vector (from strike and dip angles)
        R = rotmat(-kap1(ii),3) * rotmat(-theta1(ii),1);
        n1(:,ii) = R*upvec;
    end
    % unit vectors
    n1 = unit(n1);
    
    % output 
    Xout = n1;
    [a,b] = size(Xout);
    if a~=3, error('Xout should be 3 x n'); end
    if b~=n, error('Xout should be 3 x n'); end
  
else
    error(sprintf('itype (%i) must be 0 or 1'),itype);   
    
end

if idisplay==1
    for ii = 1:n
        disp('-------------------------------------');
        disp('index, strike, dip:');
        disp(sprintf('%6i%6.1f%6.1f',ii,kap1(ii),theta1(ii)));
        disp('fault normal vector (magnitude):');
        disp(sprintf('   n: %8.4f%8.4f%8.4f : %8.4e',n1(:,ii), sqrt(n1(1,ii)^2+n1(2,ii)^2+n1(3,ii)^2) ));
    end
end

%==========================================================================
% EXAMPLE

if 0==1
    kappa1 = 350;
    theta1 = 80;
    Xin = [kappa1 theta1]
    Xout = strdip2normal(Xin,1,1);
    Xcheck = strdip2normal(Xout,0,1)
end

%==========================================================================
