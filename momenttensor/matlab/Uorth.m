function [Uout,dtUin,dtUout] = Uorth(Uin,itype,idisplay)
%UORTH turn a nearly orthogonal rotation matrix into a rotation matrix
% 
% INPUT
%   Uin       3 x 3 x n set of bases
%   dtUin     determinant of input matrix
%   dtUout    determinant of output matrix
%
% Carl Tape, 08/10/2012
%

if ~exist('idisplay','var'), idisplay = 0; end 
if ~any(itype==[1 2 3]), error('itype = 1, 2, 3'); end
stypes = {'svd','tape','quaternion'};
stfmt = '%24.12f%20.12f%20.12f';

if idisplay==1
   disp('============');
   disp(sprintf('USING THE ORTHOGONALIZATION METHOD OF %s ',stypes{itype})); 
   disp('============');
end

[~,~,n] = size(Uin);
Uout = Uin;  % initialize
dtUin = NaN(n,1);
dtUout = NaN(n,1);
for ii=1:n
   U0 = Uin(:,:,ii); 
   switch itype
       case 1
           % svd
           [U,S,V] = svd(U0);
           Uout(:,:,ii) = U*V';
       case 2
           % TapeTape2013 suggestion
           t = Uout(:,1,ii);
           b = Uout(:,2,ii);
           p = cross(t,b);
           Uout(:,:,ii) = [t b/norm(b) p/norm(p)];
       case 3
           % quaternions
           [~,~,q,~,imaxtr] = U2q(U0,1,idisplay);
           qpick = q{imaxtr};   % q is a cell array
           Uout(:,:,ii) = q2U(qpick);
           %Uout(:,:,ii) = Ueigvec(Uout(:,:,ii));
   end
   dtUin(ii) = det(U0);
   dtUout(ii) = det(Uout(:,:,ii));
   
   if idisplay==1
        Unew = Uout(:,:,ii);
        disp('--------------------');
        disp(sprintf('input U (%i/%i):',ii,n));
        for jj=1:3, disp(sprintf(stfmt,U0(jj,:))); end
        disp('UT*U:'); oU = U0'*U0;
        for jj=1:3, disp(sprintf(stfmt,oU(jj,:))); end
        disp('det(U):'); det(U0)
        Udiff = U0 - Unew;
        disp(sprintf('ortho U using %s:',stypes{itype}));
        for jj=1:3, disp(sprintf(stfmt,Unew(jj,:))); end
        disp('UT*U:'); oU = Unew'*Unew;
        for jj=1:3, disp(sprintf(stfmt,oU(jj,:))); end
        disp('det(U):'); det(Unew)
        Udiff = U0 - Unew;
        disp('U - Unew:');
        for jj=1:3, disp(sprintf(stfmt,Udiff(jj,:))); end
        disp('norm(U-Unew):');
        norm(Udiff)
   end
end

%==========================================================================
% EXAMPLES

if 0==1
    format long
    % http://www.mathworks.com/matlabcentral/newsreader/view_thread/239674
    T = [0.9988 -0.0482 -0.0165
         0.0479 0.9988 -0.0167
         0.0173 0.0159 0.9998 ];
    [Tout,dtUin,dtUout] = Uorth(T,1)
    T,Tout
    norm(T-Tout)
    
    % example in TapeTape2013, Appendix E
    % transform from south-east-up to north-west-up
    clc, clear
    P = [-1 0 0 ; 0 -1 0 ; 0 0 1];
    U1 = P * U2pa([24 120 41 232],0);
    U2 = P * U2pa([55 295 17 51],0);
    Uorth(U1,1,1);
    Uorth(U1,2,1);
    Uorth(U1,3,1);
    
    % EXAMPLE: compare different orthogonalization methods by applying them
    % to GCMT catalog bases constructed from integer rounded trend-azimuth
    % angles of the basis vectors
    % plunge: 0 to 90
    % azimuth: 0 to 360
    clear, clc, close all
    [otime,tshift,hdur,lat,lon,dep,M,M0,Mw,eid,elabel,...
       str1,dip1,rk1,str2,dip2,rk2,lams,pl1,az1,pl2,az2,pl3,az3] = readCMT;
    n = length(otime);
    % U in south-east-up convention
    Uazpl = U2pa([pl1 az1 pl2 az2 pl3 az3],0);
    % ensure that all U are rotation matrices
    Uazpl = Udetcheck(Uazpl);
    % three types of orthogonal versions of U
    Uazplo1 = Uorth(Uazpl,1);
    Uazplo2 = Uorth(Uazpl,2);
    Uazplo3 = Uorth(Uazpl,3);
    % put all the U in the green so you can directly difference the U's
    % warning: high EPSVAL this will force more matrices toward identity
    EPSVAL  = 0.02;
    Uazpl   = Ueigvec(Uazpl,EPSVAL);
    Uazplo1 = Ueigvec(Uazplo1,EPSVAL);
    Uazplo2 = Ueigvec(Uazplo2,EPSVAL);
    Uazplo3 = Ueigvec(Uazplo3,EPSVAL);
    %% pick a subset
    npair = 10000;
    ivec1 = randi(n,npair,1);
    U1    = Uazpl(:,:,ivec1); 
    U1o1  = Uazplo1(:,:,ivec1); 
    U1o2  = Uazplo2(:,:,ivec1); 
    U1o3  = Uazplo3(:,:,ivec1); 
    stypes = {'svd','tape','quat'};
    normlim = 0.02; nedge = [0:0.001:normlim]; nmx = 0.24;
    anglelim = 1.5; aedge = [0:0.1:anglelim]; amx = 0.34;
    for kk=1:5
        if kk==1, UA=U1; UB=U1o1; slab1='U1'; slab2='Usvd'; end
        if kk==2, UA=U1; UB=U1o2; slab1='U1'; slab2='Utape'; end
        if kk==3, UA=U1; UB=U1o3; slab1='U1'; slab2='Uquat'; end
        if kk==4, UA=U1o1; UB=U1o2; slab1='Usvd'; slab2='Utape'; end
        if kk==5, UA=U1o1; UB=U1o3; slab1='Usvd'; slab2='Uquat'; end
        % compute various differences
        Udiff_L2 = norm_mat(UA-UB,2);
        Udiff_Linf = norm_mat(UA-UB,Inf);
        [omega_diff,xi_diff] = CMT2omega_xi(UA,UB,0);
        xirot_diff = U2xi(UA,UB);

        figure; nr=3; nc=2;
        subplot(nr,nc,1); plot_histo(Udiff_L2,nedge); ylim([0 nmx]);
        xlabel(['|| ' slab1 ' - ' slab2 ' ||_2']);
        %title(sprintf('GCMT catalog, %i pairs',npair));
        title(sprintf('max = %.4f (%i/%i > %.3f)',...
            max(Udiff_L2),sum(Udiff_L2>normlim),n,normlim));
        subplot(nr,nc,2); plot_histo(Udiff_Linf,nedge); ylim([0 nmx]);
        xlabel(['|| ' slab1 ' - ' slab2 ' ||_{\infty}']);
        title(sprintf('max = %.4f (%i/%i > %.3f)',...
            max(Udiff_Linf),sum(Udiff_Linf>normlim),n,normlim));
        subplot(nr,nc,3); plot_histo(xirot_diff,aedge); ylim([0 amx]);
        xlabel(['xi ( ' slab1 ', ' slab2 ' )']);
        title(sprintf('XI: max = %.4f (%i/%i > %.3f)',...
            max(xirot_diff),sum(xirot_diff>anglelim),n,anglelim));
        subplot(nr,nc,4); plot_histo(xi_diff,aedge); ylim([0 amx]);
        xlabel(['xi_0 ( ' slab1 ', ' slab2 ' )']);
        title(sprintf('XI_0: max = %.4f (%i/%i > %.3f)',...
            max(xi_diff),sum(xi_diff>anglelim),n,anglelim));
        subplot(nr,nc,5); plot_histo(omega_diff,aedge); ylim([0 amx]);
        xlabel(['omega ( ' slab1 ', ' slab2 ' )']);
        title(sprintf('OMEGA: max = %.4f (%i/%i > %.3f)',...
            max(omega_diff),sum(omega_diff>anglelim),n,anglelim));
        % print figure
        %orient tall; print(gcf,'-dpsc',sprintf('Uortho%2.2i',kk));
        %for file in `ls Uortho0*.ps` ; do ps2pdf $file ; done ; pdcat -r Uortho0*pdf Uortho.pdf
    end
    
%     ibig = find(Udiff_L2 > 0.1);
%     ipick = ibig(1);
%     Udiff_L2(ipick)
%     Utest   = U1(:,:,ipick)
%     Utesto1 = Uorth(Utest,1,1);
%     Utesto1 = Ueigvec(Utesto1,EPSVAL);
%     Utesto2 = Uorth(Utest,1,1);
%     Utesto2 = Ueigvec(Utesto2,EPSVAL);    
%     Utesto3 = Uorth(Utest,3,1);
%     %Utesto3 = Ueigvec(Utesto3,EPSVAL);
%     Utest,Utesto1,Utesto2,Utesto3
%     norm_mat(Utesto1-Utesto3,2)
    
end

%==========================================================================