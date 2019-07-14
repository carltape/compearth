%
% run_genlsq.m
%
% genlsq = generalized least squares
% This program tests several optimization methods presented in Albert
% Tarantola's book Inverse Problem Theory (2nd Ed., 2005) in Section 6.22.
% The program was written in July 2007 following Tarantola's sabbatical
% visit to Caltech. It is designed for teaching purposes.
%
% MORE INFORMATION:
% See INFO/optimization_guide.pdf
%
% NOTE: If there are constants in the Matlab in-line functions, then these
% constants are assigned at the time the function is INITIALIZED, not at
% the time the function is called.
%
% THIS PROGRAM CALLS:
%    optimization_method.m          all the iterative algorithms
%    forward_epicenter.m
%    forward_crescent.m
%    plot_histo.m                   plot a particular style of histogram
%
% Carl Tape (ctape@alaska.edu), 2010-03-20
%

clc, clear, close all
format compact
format short

% optimization method options
stlabels = {
    'none (stop here)',
    'Newton',
    'quasi-Newton',
    'steepest descent',
    'conjugate gradient',
    'conjugate gradient (polynomial line search)',
    'variable metric (matrix version)',
    'variable metric (vector version)',
    'square-root variable metric (matrix version)',
    'square-root variable metric (vector version)'
    };
stlabels2 = {'none','newton','quasi','steepest','cg','cgpoly','vmmatrix','vmvector','srvmmatrix','srvmvector'}';
nmethod0 = length(stlabels)-1;
%nmethod0 = 6;   % only first 6 are used here

%==========================================================================
% USER INPUT : choose optimization method

% set this =1 if you want to use the default values or =0 if you want
% to test different experiments
idefault = 0;

disp('  ');
disp('TYPE A NUMBER AFTER EACH PROMPT AND HIT ENTER:');
disp('Optimization problem set-up:');
iforward = input('Forward problem (1 = epicenter; 2 = epicenter-cresent): ');
%iforward = input('Forward problem (1 = epicenter; 2 = epicenter-cresent; 3 = leaf): ');
%iforward = 1;
if idefault==0
    nsamples = input('Select the number of samples for the distributions (1000): ');
    irandom_initial_model  = input('Type 1 for random initial model or 0 for fixed: ');
    irandom_target_model  = input('Type 1 for random target model or 0 for fixed: ');
    idata_errors  = input('Data errors (0 = none, 1 = random, 2 = fixed): ');
    ifig = input('Type 1 to plot figures or 0 to not: ');
else
    nsamples = 1000;
    irandom_initial_model = 0;
    irandom_target_model = 0;
    idata_errors = 2;
    ifig = 1;
end
%inormalization = input('Type 1 for normalization factors for Cd and Cm: ');
inormalization = 1;

% print figures to EPS files in directory pdir
bprint = false;
pdir = strcat(pwd,'/');
%pdir = '/home/carltape/latex/notes/tomo/figures_optim/';
if ~exist(pdir,'dir'), error('pdir does not exist'); end

%==========================================================================
% FORWARD PROBLEM    

% KEY: call the forward model
if iforward == 1, forward_epicenter; end
if iforward == 2, forward_epicenter_crescent; end
%if iforward == 3, forward_leaf; end

stnsamples = [num2str(nsamples) ' samples'];
stlabS = {'Sd(m^n)','Sm(m^n)','S(m^n) = Sd + Sm'};
stc = {'kp-','ro-','bV-','gs-','m^-','kd-','rV-','bs-','g^-','mo-'};

% predictions for prior and initial models (not necessary)
dprior = d(mprior);
dinitial = d(minitial);

if ifig==1
    if 0==1
        % plot Gaussian random samples (randn_vecs_m) based on randn (not rand)
        figure; nr=4; nc=2;
        sigma = 1;
        edges = [-4*sigma: sigma/2 : 4*sigma];
        for kk=1:nparm
            etemp = randn_vecs_m(kk,:);
            subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.301]); grid on;
            title({['Model parameter ' num2str(kk) ],
                ['mean = ' sprintf('%.5f',mean(etemp)) ...
                '; std = ' sprintf('%.5f',std(etemp)) ]});
        end
        subplot(2,1,2); hold on;
        plot(randn_vecs_m,'.-');
        title(sprintf('Gaussian distribution: %i samples',nsamples));
        xlim([0.5 nparm+0.5]);
        set(gca,'xtick',[1:nparm]);
        %set(gca,'xticklabel',mlabs);
        xlabel('model parameter'); ylabel('model parameter value'); grid on;
        if bprint, orient tall; print(gcf,'-depsc',[pdir 'randn_new']); end
    end

    %-----------------
    % plot different histograms of properties of the prior model covariance samples
    
    % display distributions for each model parameter (nparm ROWS of cov_samples_m)
    figure; nr=2; nc=2;
    for kk=1:nparm
        sigma = sigma_prior(kk);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_m(kk,:);
        subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
        title({['mprior samples: Model parameter ' num2str(kk) ' (' mlabs{kk} ')'],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    if bprint, print(gcf,'-depsc',[pdir 'mprior1']); end
    
    if 0==1
        % display model samples (nsamples COLUMNS of cov_samples_m)
        % -- this is NOT a very useful representation
        figure; nr=2; nc=1;
        subplot(nr,nc,1); plot(cov_samples_m,'.-'); xlim([0.5 nparm+0.5]);
        set(gca,'xtick',[1:nparm],'xticklabel',mlabs); grid on;
        xlabel('model parameter'); ylabel('deviation from mean model ');
        title([stnsamples ' of the prior model distribution']);
        subplot(nr,nc,2); hold on;
        plot(mprior_samples,'.-');
        plot(mprior,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k','markeredgecolor','w');
        plot(mtarget,'ro-','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
        title(' BLACK = prior model;  RED = target model');
        xlim([0.5 nparm+0.5]); set(gca,'xtick',[1:nparm],'xticklabel',mlabs);
        xlabel('model parameter'); ylabel('model parameter value'); grid on;
        if bprint, print(gcf,'-depsc',[pdir 'mprior2']); end
    end
    
    %-----------------
    % plot different histograms of properties of the data covariance samples

    % display distributions for each data index (ndata ROWS of cov_samples_d)
    figure; nr=4; nc=3;
    for ii=1:ndata
        sigma = sigma_obs(ii);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_d(ii,:);
        subplot(nr,nc,ii); plot_histo(etemp,edges); ylim([0 0.4]);
        title({['Data index ' num2str(ii)],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    orient tall, wysiwyg
    if bprint, print(gcf,'-depsc',[pdir 'CD1']); end
    
    if 0==1
        % display some length-ndata samples (nsamples COLUMNS of cov_samples_d)
        figure; plot(cov_samples_d,'.-');
        xlim([0.5 ndata+0.5]); set(gca,'xtick',[1:ndata]); grid on;
        xlabel('Data index'); ylabel('Deviation from target data ');
        title([stnsamples ' of the data distribution']);
        if bprint, print(gcf,'-depsc',[pdir 'CD2a']); end
    end
    
    figure; hold on;
    plot(dobs_samples,'.-');
    p1 = plot(dprior,'bo-','linewidth',2,'markersize',10,'markerfacecolor','b','markeredgecolor','w');
    p2 = plot(dinitial,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k','markeredgecolor','w');
    p3 = plot(dtarget,'ro--','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    p4 = plot(dobs,'ro-','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    legend([p1 p2 p3 p4],'g(mprior)','g(minitial)','g(mtarget)','g(mtarget) + errors',...
        'location','northwest');
    %title(' BLACK = d(mprior);  RED DASHED = d(mtarget);  RED = d(mtarget) + errors');
    xlim([0.5 ndata+0.5]); set(gca,'xtick',[1:ndata]);
    xlabel('Data index'); ylabel('Prediction value, g(m)'); grid on;
    if bprint, print(gcf,'-depsc',[pdir 'CD2b']); end
end

%---------------------------------------------
% MISFIT FUNCTION : least squares, Tarantola (2005), Eq. 6.251
% (This calls the function d to compute the predictions.)

if 1==1
    % data misfit
    Sd = @(m,dobs,icobs) ( 0.5* ...
          (d(m)-dobs)' * icobs * (d(m)-dobs) );
    % model misfit (related to regularization)
    Sm = @(m,mprior,icprior) ( 0.5* ...
          (m - mprior)' * icprior * (m - mprior) );
    % total misfit
    S = @(m,dobs,mprior,icobs,icprior) ( Sd(m,dobs,icobs) + Sm(m,mprior,icprior) );
    
else
    S = @(m,dobs,mprior,icobs,icprior) ( 0.5*( ...
          (d(m)-dobs)' * icobs * (d(m)-dobs) ...
        + (m - mprior)' * icprior * (m - mprior) ) );
end

% misfit for prior model and target model
%S_mtarget = S(mtarget,dobs,mprior,icobs,icprior);
%S_mprior  = S(mprior,dobs,mprior,icobs,icprior);

%///////////////////////////////
% grid search to plot misfit function
% (In many problems, the reason for using an iterative gradient-based method 
% is that we do not have the computational resources needed to evaluate the
% misfit function over the entire model space.)

if iforward==2
    % regular grid of epicenters
    numx = 50;
    [xvec,yvec,numy,a,b,dx] = gridvec(xmin,xmax,numx,ymin,ymax);
    ng = length(xvec);
    sd0 = NaN(ng,1);
    sm0 = NaN(ng,1);
    for xx=1:length(xvec)
        mg = [xvec(xx) yvec(xx)]';
        sd0(xx) = Sd(mg,dobs,icobs);
        sm0(xx) = Sm(mg,mprior,icprior);
    end
    s0 = sd0 + sm0;
    
    % plot misfit function and its two parts: S = Sm + Sd
    X = reshape(xvec,a,b);
    Y = reshape(yvec,a,b);
    figure; nr=1; nc=3; cmin = -1e-5;
    subplot(nr,nc,1); pcolor(X,Y,reshape(s0,a,b)); shading flat;
    caxis([cmin max(s0)]); colorbar; axis equal, axis(axepi); title('S(m)');
    subplot(nr,nc,2); pcolor(X,Y,reshape(sd0,a,b)); shading flat;
    caxis([cmin max(s0)]); colorbar; axis equal, axis(axepi); title('Sd(m)');
    subplot(nr,nc,3); pcolor(X,Y,reshape(sm0,a,b)); shading flat;
    caxis([cmin max(sm0)]); colorbar; axis equal, axis(axepi); title('Sm(m)');
    % for this file: ps2eps -f genlsq_misfit_f2.ps
    if bprint, print(gcf,'-dpsc',sprintf('%sgenlsq_misfit_f%i',pdir,iforward)); end
end

%///////////////////////////////

%========================================================================
% USER INPUT : choose optimization method

disp('  ');
disp('Optimization methods:');
for xx=0:nmethod0, disp(['    ' num2str(xx) ' : ' stlabels{xx+1}]); end
disp('TYPE A NUMBER AFTER EACH PROMPT AND HIT ENTER:');
disp('   IF YOU WANT MULTIPLE METHODS, LIST THE NUMBERS IN BRACKETS (e.g., [1 4 5])');
disp(['   IF YOU WANT ALL METHODS, USE [1:' num2str(nmethod0) ']']);
imethod_vec = input(['Select your optimization method(s) (1-' num2str(nmethod0) '): ']);
if imethod_vec==0, break; end
nmethod = length(imethod_vec);

niter = input('Select the number of iterations (10): ');
if any([idata_errors irandom_initial_model irandom_target_model]==1)
    nrun = input(['Select the number of different runs for each inversion method (1 <= nrun <= ' num2str(nsamples) '): ']);
else
    nrun = 1;
end
stnrun = sprintf('%4.4irun',nrun);
    
%========================================================================

Sd_mat = zeros(niter+1,nmethod);
Sm_mat = zeros(niter+1,nmethod);
S_mat = zeros(niter+1,nmethod);
Sd_array = zeros(niter+1,nmethod,nrun);
Sm_array = zeros(niter+1,nmethod,nrun);
S_array = zeros(niter+1,nmethod,nrun);
% note: these are used for summary plots if nmethod=1 and nrun > 1
mpost_models = zeros(nparm,nrun);
mtarget_models = zeros(nparm,nrun); 
minitial_models = zeros(nparm,nrun);    

% LOOP OVER DIFFERENT RUNS (VARYING INITIAL MODEL, TARGET MODEL, OR DATA ERRORS)
for irun = 1:nrun
    disp('  '); disp(['Run ' num2str(irun) ' out of ' num2str(nrun)]);
    
    % compute a set of observations for a randomly picked target model,
    % and add in the appropriate errors to the target data
    %if or(irandom==1, nmodel > 1)
    %    mtarget = mprior_samples(:,irun); 
    %    dtarget = d(mtarget);
    %    eobs = cov_samples_d(:,irun);
    %    dobs = dtarget + eobs;
    %end
    
    % pick initial model
    if irandom_initial_model==1
        iinitial = irun;            % draw from the START of mprior samples
        minitial = mprior_samples(:,iinitial); 
    end
    minitial_models(:,irun) = minitial;
    
    % pick target model
    if irandom_target_model==1 
        itarget = nsamples+1-irun;  % draw from the END of mprior samples
        mtarget = mprior_samples(:,itarget);
        dtarget = d(mtarget);
    end
    mtarget_models(:,irun) = mtarget;
    
    % pick errors for target model
    if idata_errors==1
        eobs = cov_samples_d(:,irun);
    end
    dobs = dtarget + eobs;
    
    % LOOP OVER DIFFERENT METHODS
    for zz = 1:nmethod
        imethod = imethod_vec(zz);

        % preconditional F
        F0 = eye(nparm);

        % initial model
        %mnew = mprior;     % prior model
        mnew = minitial;
        Sd_0 = Sd(mnew,dobs,icobs);
        Sm_0 = Sm(mnew,mprior,icprior);
        S_0 = S(mnew,dobs,mprior,icobs,icprior);
        stS0 = sprintf(' S(m0) = %.3f = %.3f(D) + %.3f(M)',S_0,Sd_0,Sm_0);

        % initialize arrays
        iter_vec = [0:niter]';
        Sd_vec = zeros(niter+1,1);
        Sm_vec = zeros(niter+1,1);
        S_vec = zeros(niter+1,1);
        
        % misfit for initial model
        Sd_vec(1) = Sd_0;
        Sm_vec(1) = Sm_0;
        S_vec(1) = S_0;
        
        % KEY: perform the iterative inversion
        disp('==========================================================');
        disp(['Running ' stlabels{imethod+1} '...']); 
        
        disp('==========================================================');
        optimization_method;
        %optimization_method_ex;
        
        %==================================================================

        % NEW: incorporate normalization factors into Cm and Cd
        % OLD: normalize misfit vector to the initial value
        %S_vec = S_vec/S_vec(1);
        
        % limits for y-axis on convergence plots
        ylims = 10.^[-1 2];
        if iforward==2, ylims = 10.^[-2 0.5]; end
        %ylims = [6 10];

        %--------------------------
        % posterior model and covariance matrix
        % (Regardless of the method, we can compute the posterior covariance
        % matrix, since we have the matrix of Frechet derivatives evaluated for
        % the final model.)

        % posterior model
        mpost = mnew;
        dpost = d(mpost);

        % posterior covariance matrix (e.g., Tarantola Eq. 3.53)
        % note: cpost0 does not include normalization factors Cdfac and Cmfac
        Gpost = G(mpost);
        cpost0 = inv(Gpost'*icobs0*Gpost + icprior0);
        %cpost0 = (cpost0 + cpost0')/2;    % force to be symmetric

        if exist('Fhat','var')
            disp('compare Fhat with Cpost:');
            Fhat
            cpost = inv(Gpost'*icobs*Gpost + icprior)
        end
        
        % store mpost
        %if and(nmodel > 1,nmethod==1)
        if nmethod==1
            mpost_models(:,irun) = mpost;
        end
        
        % using matlab functions or not
        if 0==1
            % posterior model uncertainties
            for kk=1:nparm, sigma_post(kk) = cpost0(kk,kk)^(1/2); end

            % a posteriori model correlations
            rho_post = zeros(nparm,nparm);
            for k1=1:nparm
                for k2=1:nparm
                    rho_post(k1,k2) = (cpost0(k1,k2)/(sigma_post(k1)*sigma_post(k2)));
                end
            end
        else
            sigma_post = sqrt( diag(cpost0) );
            rho_post = corrcov(cpost0);
        end
        
        % a priori model correlations (for comparison)
        rho_prior = corrcov(cprior0);

        cpost0, rho_post

        % SEE optim_hw.m FOR MORE DETAILS
        % posterior data covariance matrix (e.g., Tarantola Eq. 3.44)
        %cpost0_D = Gpost*cpost0*Gpost';
        %rho_post_D = corrcov(cpost0_D);     % posterior correlation matrix
        %rho_prior_D = corrcov(cobs0);       % prior, for comparison
        
        %format long
        disp(sprintf('model summary (%i iterations):',niter));
        disp('    prior    initial  posterior   target');
        disp([mprior minitial mpost mtarget]);
        disp(sprintf('data summary (%i observations):',ndata));
        disp('    prior    initial   posterior   target   actual');
        disp([dprior dinitial dpost dtarget dobs]);

        % Cholesky decomposition to obtain the square-root of cpost0
        % NOTE: for large problems, this is not possible due to poor
        %       conditioning of cpost0 or the inability to compute cpost0
        Lpost = chol(cpost0,'lower');
        
        % samples of the posterior distribution
        if or(imethod >= 1, imethod <= nmethod0)
            
            mcov_samples = zeros(nparm,nsamples);
            randn_vecs_m = randn(nparm,nsamples);
            if or(imethod >= 1, imethod <= nmethod0)
                mcov_samples  = Lpost * randn_vecs_m;
            elseif or(imethod==8, imethod==9)
                % SQUARE-ROOT VARIABLE METRIC METHOD
                for kk = 1:nsamples
                    chi = randn_vecs_m(:,kk);
                    mcov_samples(:,kk) = srvm_Shat_chi(chi,niter-1,S0hat,nu,a,w,0);
                end  
            end
            mpost_samples = repmat(mpost,1,nsamples) + mcov_samples;

            % compare the standard deviation with sigma_post
            std_samples = std(mpost_samples');

            % compare posterior model distribution with prior
            disp('  ');
            disp(' Compare model uncertainties : ');
            disp(['             model parameter : ' sprintf('%13s', mlabs{:}) ]);
            disp(['                       units : ' sprintf('%13s', ulabs{:}) ]);
            disp(['                 sigma_prior = ' sprintf('%13.5f', sigma_prior) ]);
            disp(['                  sigma_post = ' sprintf('%13.5f', sigma_post) ]);
            disp(['   std(' sprintf('%6.0f', nsamples) ' mpost_samples) = ' sprintf('%13.5f', std_samples) ]);
            disp(['    sigma_post / sigma_prior = ' sprintf('%13.5f', sigma_post./sigma_prior) ]);
            disp('  ');
            
            if ifig==1
                ftag = [stlabels2{imethod+1} '_' sprintf('run%4.4i',irun) ];
                
                if 0==1
                    % this representation is not very useful
                    figure; plot(mcov_samples,'.-'); xlim([0.5 nparm+0.5]);
                    set(gca,'xtick',[1:nparm],'xticklabel',mlabs); grid on;
                    xlabel('model parameter'); ylabel('deviation from mean ');
                    title({stlabels{imethod},[num2str(nsamples) ' samples of the posterior']});
                    if bprint, print(gcf,'-depsc',[pdir 'mpost1a_' stlabels2{imethod+1}]); end
                    
                    figure; hold on; plot(mpost_samples,'.-');
                    p1 = plot(mprior,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k','markeredgecolor','w');
                    p2 = plot(mtarget,'ro-','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
                    p3 = plot(mpost,'bo-','linewidth',2,'markersize',10,'markerfacecolor','b','markeredgecolor','w');
                    legend([p1 p2 p3],'mprior','mtarget','mpost','location','northwest')
                    xlim([0.5 nparm+0.5]); set(gca,'xtick',[1:nparm],'xticklabel',mlabs);
                    xlabel('model parameter'); ylabel('model parameter value'); grid on;
                    if bprint, print(gcf,'-depsc',[pdir 'mpost1b_' ftag]); end
                end

                % display distributions for each model parameter (nparm ROWS of cov_samples_m)
                figure; nr=2; nc=2;
                for kk=1:nparm
                    sigma = sigma_post(kk);
                    edges = [-4*sigma: sigma/2 : 4*sigma];
                    etemp = mcov_samples(kk,:);
                    subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
                    stl1 = ['mpost samples for ' stlabels2{imethod+1} ' method'];
                    stl2 = ['Model parameter ' num2str(kk) ' (' mlabs{kk} ')'];
                    stl3 = ['mean = ' sprintf('%.5f',mean(etemp)) ...
                        '; std = ' sprintf('%.5f',std(etemp)) ];
                    if kk==1, title({stl1,stl2,stl3});
                    else title({stl2,stl3}); end
                end
                if bprint, print(gcf,'-depsc',[pdir 'mpost2_' ftag]); end

                % epicenter problems only
                if or(iforward==1, iforward==2)
                    % opts is set in forward_epicenter
                    plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost);
                    % plot the cpost0 samples and re-plot the two markers
                    plot(mpost_samples(1,:),mpost_samples(2,:),'c.');
                    plot(mpost(1),mpost(2),'o','markersize',10,'markerfacecolor','c','markeredgecolor','w');
                    plot(mtarget(1),mtarget(2),'o','markersize',10,'markerfacecolor','r','markeredgecolor','w');
                    title({stlabels{imethod+1},
                        sprintf('samples of prior (blue) and posterior (cyan), run %i out of %i',irun,nrun)});
                    %orient tall, wysiwyg
                    if bprint, print(gcf,'-depsc',[pdir 'mpost_' ftag '_epi']); end
                end
                
                % convergence curve
                ya = floor(min(log10([Sm_vec ; Sd_vec ; S_vec])));
                yb = ceil(max(log10([Sm_vec ; Sd_vec ; S_vec])));
                ylims = 10.^[ya yb];
                stit = [stlabels{imethod+1} ': ' num2str(niter) ' iterations, ' ...
                    sprintf('run %i out of %i',irun,nrun) ];
                figure; hold on;
                plot(iter_vec, log10(Sd_vec),'r.-',iter_vec, log10(Sm_vec),'b.-',iter_vec, log10(S_vec),'k.-',...
                    'linewidth',2,'markersize',20);
                legend(stlabS); xlim([-0.5 niter+0.5]); ylim(log10(ylims));
                set(gca,'xtick',[-1:niter+1]); grid on;
                xlabel('n, iteration'); ylabel(' log10[ S(m^n) ], misfit function'); title(stit);
                if bprint, print(gcf,'-depsc',[pdir 'converge_' ftag]); end
            end
            
        else
            error('check imethod');
        end

        %--------------------------
        % misfit function values

        disp(sprintf('%8s%16s%16s%16s','iter','Sd','Sm','S = Sm + Sd'));
        for nn=1:niter+1, disp(sprintf('%8i%16.10f%16.10f%16.10f',iter_vec(nn),Sd_vec(nn),Sm_vec(nn),S_vec(nn))); end

%         if and(ifig==1, nmethod==1)
%            figure; semilogy(iter_vec, S_vec,'b.','markersize',24);
%            xlim([-0.5 niter+0.5]); ylim(ylims);
%            set(gca,'xtick',[-1:niter+1]); grid on;
%            xlabel('n, iteration'); ylabel(' S(m^n), misfit function');
%            title([stlabels{imethod+1} ' method with ' num2str(niter) ' iterations']);
%         end
        
        % plot figure for a SINGLE method for a SINGLE run
        if ifig==1

        end

        Sd_mat(:,zz) = Sd_vec;
        Sm_mat(:,zz) = Sm_vec;
        S_mat(:,zz) = S_vec;
        
    end  % for zz=1:imethod
    
    % store all misfit curves
    Sd_array(:,:,irun) = Sd_mat;
    Sm_array(:,:,irun) = Sm_mat;
    S_array(:,:,irun) = S_mat;

    % superimpose convergence plots for difference methods
    if and(ifig==1, nmethod > 1)
        
        figure; hold on;
        for zz=1:nmethod
            plot(iter_vec, log10(S_mat(:,zz)), stc{zz},'markersize',12);
        end
        legend(stlabels{imethod_vec+1});
        xlim([-0.5 niter+0.5]); ylim(log10(ylims));
        set(gca,'xtick',[-1:niter+1]); grid on;
        xlabel('n, iteration'); ylabel('log10[ S(m^n) ], misfit function');
        title(stS0);
        %orient tall, wysiwyg
        if bprint, print(gcf,'-depsc',[pdir 'converge_Nmethod_' stnrun]); end
        %ylim([-0.1 1.3]); print(gcf,'-depsc',['converge_Nmethod_' stnmodel]);
        
%         % variable metric methods (check that ALL FOUR ARE the same)
%         if imethod==9
%             stc = {'kp-','ro-','bV-','gx-','m^-'};
%             figure; hold on;
%             for zz=5:8
%                 plot(iter_vec, log10(S_mat(:,zz)), stc{zz-4},'markersize',12);
%             end
%             legend(stlabels{6},stlabels{7},stlabels{8},stlabels{9});
%             xlim([-0.5 niter+0.5]); ylim(log10(ylims));
%             set(gca,'xtick',[-1:niter+1]); grid on;
%             xlabel('n, iteration'); ylabel('log10[ S(m^n) ], misfit function');
%             title(stS0);
%             if bprint, print(gcf,'-depsc',[pdir 'variable_metric']); end
%         end
    end
    
%     % plot epicenters
%     if and(and(ifig==1, iforward==1),nmethod==1)
%         plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost);
%         title(sprintf('run %i out of %i',irun,nrun));
%         orient tall, wysiwyg
%         if bprint, print(gcf,'-depsc',[pdir 'mpost_epi']); end
%     end

end  % for irun

%======================================================================
% SUMMARY FIGURES

if and(nrun==1, nmethod==1)
    if ifig==0, disp('use ifig=1 to compare one method for a single run'); end
    
%     % plot figure for a SINGLE method for a SINGLE run
%     stit = [stlabels{imethod+1} ': ' num2str(niter) ' iterations, ' num2str(nrun) ' runs'];
% 
%     figure; hold on;
%     plot(iter_vec, log10(Sd_array),'r.-',iter_vec, log10(Sm_array),'b.-',iter_vec, log10(S_array),'k.-',...
%         'linewidth',2,'markersize',20);
%     legend(stlabS); xlim([-0.5 niter+0.5]); ylim(log10(ylims));
%     set(gca,'xtick',[-1:niter+1]); grid on;
%     xlabel('n, iteration'); ylabel(' log10[ S(m^n) ], misfit function'); title(stit);
%     if bprint, print(gcf,'-depsc',[pdir 'converge_' stlabels2{imethod+1} '_1run']); end
    
elseif and(nrun > 1, nmethod==1)
    % plot figure for a SINGLE method for multiple runs
    stit = sprintf('%s: %i iterations, %i runs',stlabels{imethod+1},niter,nrun);
    
    % compute mean curves
    sumcurve1 = zeros(niter+1,1); sumcurve2 = zeros(niter+1,1); sumcurve3 = zeros(niter+1,1);
    for xx=1:nrun
        temp1 = log10(Sd_array(:,1,xx));
        temp2 = log10(Sm_array(:,1,xx));
        temp3 = log10(S_array(:,1,xx));
        sumcurve1 = sumcurve1 + temp1;
        sumcurve2 = sumcurve2 + temp2;
        sumcurve3 = sumcurve3 + temp3;
    end
    
    figure; nr=2; nc=2;
    for xx=1:10
        temp1 = log10(Sd_array(:,1,xx));
        temp2 = log10(Sm_array(:,1,xx));
        temp3 = log10(S_array(:,1,xx));
        subplot(nr,nc,1); hold on; plot(iter_vec, temp1,'m');
        subplot(nr,nc,2); hold on; plot(iter_vec, temp2,'c');
        subplot(nr,nc,3); hold on; plot(iter_vec, temp3,'color',0.5*[1 1 1]);
        %subplot(nr,nc,4); hold on; plot(iter_vec, temp1, 'r',iter_vec, temp2, 'b',iter_vec, temp3, 'k');
    end
    subplot(nr,nc,1); plot(iter_vec, sumcurve1/nrun,'r','linewidth',4);
    subplot(nr,nc,2); plot(iter_vec, sumcurve2/nrun,'b','linewidth',4);
    subplot(nr,nc,3); plot(iter_vec, sumcurve3/nrun,'k','linewidth',4);
    subplot(nr,nc,4); hold on;
    p1 = plot(iter_vec, sumcurve1/nrun,'r','linewidth',4);
    p2 = plot(iter_vec, sumcurve2/nrun,'b','linewidth',4);
    p3 = plot(iter_vec, sumcurve3/nrun,'k','linewidth',4);
    legend([p1 p2 p3],stlabS);
    % labels and axis limits for each subplot
    for xx=1:4
        subplot(nr,nc,xx);   
        xlim([-0.5 niter+0.5]); ylim(log10(ylims));
        set(gca,'xtick',[-1:niter+1]); grid on;
        xlabel('n, iteration');
        if xx==4, ylabel(sprintf('log10[ %s ]',stlabS{3})); title(stit);
        else ylabel(sprintf('log10[ %s ]',stlabS{xx})); end
    end
    orient tall, wysiwyg
    if bprint
        pfile = sprintf('%sconverge_%s_%s_derr%i',pdir,stlabels2{imethod+1},stnrun,idata_errors);
        print(gcf,'-depsc',pfile);
    end

    % modified version of plot_epicenter.m
    % NOTE: we are plotting mpost for each run, but NOT the associated
    %       cpost0 samples for each run
    if or(iforward==1, iforward==2)
        figure; hold on;
        for xx=1:nrun
            % line segment connecting initial (x,y) to post (x,y)
            plot([minitial_models(1,xx) mpost_models(1,xx)],[minitial_models(2,xx) mpost_models(2,xx)],'k');
        end
        p0 = plot(mprior_samples(1,:),mprior_samples(2,:),'.');
        p1 = plot(minitial_models(1,:),minitial_models(2,:),'o','markersize',10,'markerfacecolor','k','markeredgecolor','w');
        p2 = plot(mpost_models(1,:),mpost_models(2,:),'o','markersize',10,'markerfacecolor','g','markeredgecolor','w');
        pP = plot(mprior(1),mprior(2),'o','markersize',10,'markerfacecolor','b','markeredgecolor','w');
        pT = plot(mtarget_models(1,:),mtarget_models(2,:),'o','markersize',10,'markerfacecolor','r','markeredgecolor','w');
        plot(xrec,yrec,'kV','markersize',10,'Linewidth',2);
        for ii=1:ndata, text(xrec(ii),yrec(ii),num2str(ii)); end
        legend([p0(1) p1(1) p2(1) pP pT],'Cprior sample','m00',sprintf('m%2.2i',niter),'mprior','mtarget');
        axis equal; axis(axepi); grid on;
        %set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
        xlabel('X distance (km)'); ylabel('Y distance (km)'); title(stit);
        orient tall, wysiwyg
        if bprint, print(gcf,'-depsc',[pfile '_epi']); end
    end

elseif and(nrun > 1, nmethod > 1)
    % compare the convergence results for all nrun data sets
    % NOTE: here we do NOT separate S = Sd + Sm parts
    figure; nc=2; nr=ceil(nmethod/2);
    meancurves = zeros(niter+1,nmethod);
    for zz = 1:nmethod
        subplot(nr,nc,zz); hold on;
        sumcurve = zeros(niter+1,1);
        for xx=1:nrun
            temp = log10(S_array(:,zz,xx));
            plot(iter_vec, temp, 'k');
            sumcurve = sumcurve + temp;
        end
        meancurves(:,zz) = sumcurve/nrun;
        plot(iter_vec, meancurves(:,zz),'r','linewidth',2);
        xlim([-0.5 niter+0.5]); ylim(log10(ylims));
        set(gca,'xtick',[-1:niter+1]); grid on;
        xlabel('n, iteration'); ylabel('log10[ S(m^n) ], misfit function');
        title([stlabels{zz+1} ', ' num2str(nrun) ' runs']);
    end
    orient tall, wysiwyg
    if bprint, print(gcf,'-depsc',[pdir 'converge_Nmethod_' stnrun '_all']); end

    figure; hold on;
    for zz=1:nmethod, plot(iter_vec, meancurves(:,zz), stc{zz},'markersize',12); end
    legend(stlabels(imethod_vec+1));
    xlim([-0.5 niter+0.5]); ylim(log10(ylims));
    set(gca,'xtick',[-1:niter+1]); grid on;
    xlabel('n, iteration'); ylabel('log10[ S(m^n) ], misfit function');
    title(['MEANS for each method, ' num2str(nrun) ' runs']);
    if bprint, print(gcf,'-depsc',[pdir 'converge_Nmethod_' stnrun]); end
    
else
    if ifig==0, disp('use ifig=1 to compare multiple methods for a single run'); end
end

%==========================================================================
