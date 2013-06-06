function [xi,ixi,q,trALL,imaxtr] = U2q(U,itrace,idisplay)
%U2Q convert rotation matrix to unit quaternion
%
% INPUT
%   U         3 x 3 x n set of rotation matrices
%   itrace    =1 to use quaternions for orthogonalization
%   idisplay  optional: display details (=1)
%
% OUTPUT   
%   xi     n x 1 set of kagan angles
%   ixi    index into subset zone (=1 for w; =2 for x, =3 for y, =4 for z)
%   q      structure with the following:
%            qw     4 x n set of (near) unit quaternions
%            qx     4 x n set of (near) unit quaternions
%            qy     4 x n set of (near) unit quaternions
%            qz     4 x n set of (near) unit quaternions
% 
% Carl Tape, 8/12/2012
%

if ~exist('idisplay','var'), idisplay = 0; end 

[~,~,n] = size(U);

u11 = squeeze(U(1,1,:));
u12 = squeeze(U(1,2,:));
u13 = squeeze(U(1,3,:));
u21 = squeeze(U(2,1,:));
u22 = squeeze(U(2,2,:));
u23 = squeeze(U(2,3,:));
u31 = squeeze(U(3,1,:));
u32 = squeeze(U(3,2,:));
u33 = squeeze(U(3,3,:));

% reshape to row vectors
u11 = u11(:)'; u12 = u12(:)'; u13 = u13(:)';
u21 = u21(:)'; u22 = u22(:)'; u23 = u23(:)';
u31 = u31(:)'; u32 = u32(:)'; u33 = u33(:)';

qw = zeros(4,n);
qx = zeros(4,n);
qy = zeros(4,n);
qz = zeros(4,n);

w = 0.5 * sqrt(1 + u11 + u22 + u33);
x = 0.5 * sqrt(1 + u11 - u22 - u33);
y = 0.5 * sqrt(1 - u11 + u22 - u33);
z = 0.5 * sqrt(1 - u11 - u22 + u33);

qw(1,:) = w;
qw(2,:) = (u32-u23)./(4*w);
qw(3,:) = (u13-u31)./(4*w);
qw(4,:) = (u21-u12)./(4*w);

qx(1,:) = (u32-u23)./(4*x);
qx(2,:) = x;
qx(3,:) = (u12+u21)./(4*x);
qx(4,:) = (u13+u31)./(4*x);

qy(1,:) = (u13-u31)./(4*y);
qy(2,:) = (u12+u21)./(4*y);
qy(3,:) = y;
qy(4,:) = (u23+u32)./(4*y);

qz(1,:) = (u21-u12)./(4*z);
qz(2,:) = (u13+u31)./(4*z);
qz(3,:) = (u23+u32)./(4*z);
qz(4,:) = z;

if itrace==1
   Xpi = diag([ 1 -1 -1]);
   Ypi = diag([-1  1 -1]);
   Zpi = diag([-1 -1  1]);
   trALL = NaN(4,n);
   for kk=1:n
       U0 = U(:,:,kk);
       trALL(:,kk) = [sum(diag(U0))     ; sum(diag(U0*Xpi)) ;
                      sum(diag(U0*Ypi)) ; sum(diag(U0*Zpi)) ];
   end
   [~,imaxtr] = max(trALL);
end

% compute the kagan angle (xi0 in TapeTape2013)
%[val,ixi] = max(abs(qw));
[val,ixi] = max([w ; x ; y; z]);
xi = 2*acos(val)*180/pi;
xi = xi(:);

% rotation angle (xi in TapeTape2013)
% note: this is currently not returned as an output argument
xirot = 2*acos(w)*180/pi;

if and(any(~isreal([w x y z])),itrace==0)
    disp('WARNING (U2q.m): imaginary entry in w, x, y, z');
    %error('imaginary entry in w, x, y, z');
end

if idisplay==1
    stcol = {'w','x','y','z'};
    stfmt = '%24.12f%20.12f%20.12f';
    for kk=1:n
        U0 = U(:,:,kk);
        disp(sprintf('input U (%i/%i):',kk,n));
        for ii=1:3, disp(sprintf(stfmt,U0(ii,:))); end
        disp('UT*U:'); UTU = U0'*U0;
        for ii=1:3, disp(sprintf(stfmt,UTU(ii,:))); end
        disp('det(U):'); det(U0)
        if itrace==1
            disp('computing trace to pick the best choice for orthogonalization:');
            disp(sprintf('  tr(U)     = %10.4f ',trALL(1,kk)));
            disp(sprintf('  tr(U*Xpi) = %10.4f ',trALL(2,kk)));
            disp(sprintf('  tr(U*Ypi) = %10.4f ',trALL(3,kk)));
            disp(sprintf('  tr(U*Zpi) = %10.4f ',trALL(4,kk)));
            disp(sprintf('  --> max is for %s',stcol{imaxtr(kk)}));
        end
        if any(~isreal([w(kk) x(kk) y(kk) z(kk)]))
            disp('WARNING: imaginary entry in w, x, y, z');
            disp('w = '); w(kk)
            disp('x = '); x(kk)
            disp('y = '); y(kk)
            disp('z = '); z(kk)
            disp('quaternions (REAL PART ONLY):');
        else
            disp('quaternions:');
        end
        %whos qw qx qy qz
        qdisp = [qw(:,kk) qx(:,kk) qy(:,kk) qz(:,kk)]';
        disp(sprintf('%15s%9s%9s%9s','qw','qx','qy','qz'));
        for ii=1:4, disp(sprintf('   q%i %9.4f%9.4f%9.4f%9.4f',ii,qdisp(ii,:))); end
        disp(sprintf('%13s%9s%9s%9s','w','x','y','z'));
        disp(sprintf('%15.4f%9.4f%9.4f%9.4f',w(kk),x(kk),y(kk),z(kk)));
        disp(sprintf('max(|w|,|x|,|y|,|z|) = %.4f',val(kk)));
        disp(sprintf('      rotation angle = %.3f deg (w)',xirot(kk)));
        disp(sprintf('principal axis angle = %.3f deg (%s)',xi(kk),stcol{ixi(kk)}));
        disp('--------------------');
    end
end

% cell array for quaternion output
q = {qw,qx,qy,qz};

%==========================================================================
% EXAMPLES

if 0==1
    clear, close all, clc
    % example bases from Kagan (1991)
    % U1 and U2 are NOT quite orthogonal
    % optional: transform from south-east-up to north-west-up
    P = [-1 0 0 ; 0 -1 0 ; 0 0 1];
    U1 = P * U2pa([24 120 41 232],0); % not orthogonal
    U2 = P * U2pa([55 295 17 51],0);  % not orthogonal
    U1o1 = Uorth(U1,1);
    U2o1 = Uorth(U2,1);
    U1o2 = Uorth(U1,2);
    U2o2 = Uorth(U2,2);
    disp('INPUT U1:'); U2q(U1,0,1);
    disp('INPUT U1o1 (svd):'); U2q(U1o1,0,1);
    disp('INPUT U1o2 (tape):'); U2q(U1o2,0,1);
    disp('INPUT U2:'); U2q(U2,1,1);
    disp('INPUT U2o1 (svd):'); U2q(U2o1,0,1);
    disp('INPUT U2o2 (tape):'); U2q(U2o2,0,1);
    disp('now consider the matrix U = U1o2T * U2o2 (tape):');
    U = U1o2'*U2o2; U2q(U,0,1);
    %xi = U2xi(U1o2,U2o2)  % check
    
    %---------------------------------------
    % testing orthogonalization
    disp('INPUT U1:');
    [xi,ixi,q,trALL,imaxtr] = U2q(U1,1,1);
    qpick = q{imaxtr};
    qpick = qpick/norm(qpick);
    
    %---------------------------------------
    % random rotation matrices
    n = 1000;
    U = randomU(n);
    xi = U2q(U,0,0); figure; plot_histo(xi,[0:5:120]);
    xlabel('xi, principal axis angle'); title('random rotation matrices');
    
    %---------------------------------------
    % try several sets
    % here we use each GCMT basis to represent U = U1'*U2
    isub = 1:10000;
    [otime,tshift,hdur,lat,lon,dep,M] = readCMT(isub);
    n = length(otime);
    [~,U] = CMTdecom(M);
    U = Udetcheck(U);   % ensure that U are rotation matrices
    xi = U2q(U,0); figure; plot_histo(xi,[0:5:120]);
    xlabel('xi, principal axis angle');
end

%==========================================================================
