% 
% TT2012kaganAppE.m
%
% Reproduce the results in Appendix E of Tape and Tape (2012c)
%   "Angle between principal axis triples"
% These are the same two events used in Kagan (1991).
%
% Carl Tape, 10/22/2013
%

close all, clear, clc

%-----------------------------------------------------------------------

% % read events from GCMT catalog
% ax3 = [140 150 -10 10 -10 700];
% oran = [datenum(1977,1,6) datenum(1977,1,7)];
% [otime1,~,~,lat1,lon1,dep1,M1,M01,Mw1,eid1] = readCMT(oran,ax3,[]);
% oran = [datenum(1980,9,26) datenum(1980,9,27)];
% [otime2,~,~,lat2,lon2,dep2,M2,M02,Mw2,eid2] = readCMT(oran,ax3,[]);

% Appendix E of TapeTape2012 "Angle betweeen principal axis triples"
% Kagan (1991) comparison events from GCMT catalog
eid1 = 'C010677A'; eid2 = 'C092680B';   % New Guinea
fac1 = 1e19; fac2 = 1e18;
if 0==1     % if you have access to the full GCMT catalog
    % read full GCMT catalog
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid,elabel,...
        str1,dip1,rk1,str2,dip2,rk2,lams,pl1,az1,pl2,az2,pl3,az3] = readCMT;
    %eid1 = 'B060586B'; eid2 = 'B062486F';  % south Pacific
    ieid1 = find(strcmp(eid1,eid)==1);
    ieid2 = find(strcmp(eid2,eid)==1);
    % event 1
    disp(sprintf('%s %s',eid{ieid1},datestr(otime(ieid1),29)));
    disp(sprintf('str/dip/rk (1) = %i/%i/%i',str1(ieid1),dip1(ieid1),rk1(ieid1)));
    disp(sprintf('str/dip/rk (2) = %i/%i/%i',str2(ieid1),dip2(ieid1),rk2(ieid1)));
    disp(sprintf('pl1/az1 = %i/%i, pl2/az2 = %i/%i, pl3/ax3 = %i/%i',...
        pl1(ieid1),az1(ieid1),pl2(ieid1),az2(ieid1),pl3(ieid1),az3(ieid1)));
    M1 = M(:,ieid1);
    % event 2
    disp(sprintf('%s %s',eid{ieid2},datestr(otime(ieid2),29)));
    disp(sprintf('str/dip/rk (1) = %i/%i/%i',str1(ieid2),dip1(ieid2),rk1(ieid2)));
    disp(sprintf('str/dip/rk (2) = %i/%i/%i',str2(ieid2),dip2(ieid2),rk2(ieid2)));
    disp(sprintf('pl1/az1 = %i/%i, pl2/az2 = %i/%i, pl3/ax3 = %i/%i',...
        pl1(ieid2),az1(ieid2),pl2(ieid2),az2(ieid2),pl3(ieid2),az3(ieid2)));
    M2 = M(:,ieid2);
else
    % New Guinea
    M1 = [-0.4120    0.0840    0.3280    0.3980   -1.2390    1.0580]'*1e19;
    M2 = [ 5.0540   -2.2350   -2.8190   -0.4760    5.4200    5.5940]'*1e18;
end
% event 1
disp('///////////////// EVENT 1 /////////////////');
Mvec2Mmat(M1/fac1,1)
disp(sprintf('multiplication factor is %.0e N-m',fac1));
[lam,U] = CMTdecom(M1); U = Udetcheck(U); U14 = Ufour(U);
% note that FRAME 4 will match the U1 listed in TapeTape Appendix E
for kk=1:4
    disp(sprintf('========= FRAME %i ===========',kk));
    U1 = U14(:,:,kk);
    U1, det(U1), diag(lam), U1*diag(lam)*U1', Mvec2Mmat(M1,1)
end
% event 2
disp('///////////////// EVENT 2 /////////////////');
Mvec2Mmat(M2/fac2,1)
disp(sprintf('multiplication factor is %.0e N-m',fac2));
[lam,U] = CMTdecom(M2); U = Udetcheck(U); U24 = Ufour(U);
for kk=1:4
    disp(sprintf('========= FRAME %i ===========',kk));
    U2 = U24(:,:,kk);
    U2, det(U2), diag(lam), U2*diag(lam)*U2', Mvec2Mmat(M2,1)
end
    
% now consider the 'difference' matrix U = U1^-1 U2
% the values for U and q should match the results in TapeTape2012, App. E
U1 = U14(:,:,4);
U2 = U24(:,:,4);
U12 = U1'*U2;
[omega,xi0] = CMT2omega_xi0(M1,M2,0,1);
[xi0,ixi0,q] = U2xi0(U12,1,1);     % more info

% four possible quaternions
Uall = Ufour(U12);
[~,~,qall] = U2xi0(Uall,0,1);

disp('========== FULL RESULTS FROM APPENDIX E OF TAPE AND TAPE (2012) ===========');
disp(sprintf('first moment tensor with scale factor (%.3e) removed:',fac1))
disp('M1:'); disp(Mvec2Mmat(M1/fac1,1))
U1, U2
disp('U1^T * U2:'); disp(U12);
disp('one of four possible quaternions for U1^T*U2:'); disp(qall(:,1))
xi0
disp('============================================================================');

% checking output on 10/23/2013

% ========== FULL RESULTS FROM APPENDIX E OF TAPE AND TAPE (2012) ===========
% first moment tensor with scale factor (1.000e+19) removed:
% M1:
%    -0.4120    0.3980   -1.2390
%     0.3980    0.0840    1.0580
%    -1.2390    1.0580    0.3280
% U1 =
%     0.4041    0.6425    0.6511
%    -0.4565    0.7585   -0.4651
%    -0.7927   -0.1093    0.5998
% U2 =
%    -0.8217    0.4860   -0.2978
%    -0.2373   -0.7667   -0.5965
%    -0.5183   -0.4195    0.7453
% U1^T * U2:
%     0.1871    0.8789   -0.4388
%    -0.6512   -0.2235   -0.7252
%    -0.7355    0.4215    0.5305
% one of four possible quaternions for U1^T*U2:
%     0.6112
%     0.4691
%     0.1213
%    -0.6259
% xi0 =
%   102.5062
% ============================================================================
