%
% [dlon,dlat,d,dsig,ax0,slabel,ulabel,iall] = get_1D_dataset_carl(ropt,dopt)
% Carl Tape, 01-Feb-2011
%
% This loads a set of discrete (1D) observations on the sphere
% that we will use to estimate a smooth continuous field.
%
% calls xxx
% called by sphereinterp.m
%

function [dlon,dlat,d,dsig,ax0,slabel,ulabel,iall] = get_1D_dataset_carl(ropt,dopt)

% GEOGRAPHIC REGION 
switch ropt
    case 1, rlabel = 'socal'; ax0 = [-122 -113 30 38]; szone = '11S';
    %case 2, rlabel = 'cal'; ax0 = [-128.05 -112.95 30.95 43.05]; szone = '11S';
    case 2, rlabel = 'cal'; ax0 = [-128.05 -111.95 29.95 43.05]; szone = '11S';
    case 3, rlabel = 'maricopa'; ax0 = [-119.7 -118.5 34.5 35.5]; szone = '11S';    
    case 4, rlabel = 'nenana'; ax0 = [-151.05 -147.45 63.45 65.55]; szone = '6W';
    case 5, rlabel = 'alaska'; ax0 = [190 230 54 72]; szone = '5V';        
end

%========================================================
% LOAD OBSERVATIONS

if dopt == 1
    %---------------------------
    % load Moho points from all available sources
    % 1: Crust2.0 global modelcd CV
    % 2: EARS -- automated receiver functions
    % 3: Magistrale et al. (2000) -- for Salton trough
    % 4: Bryant and Jones (1992) -- for western Transverse Ranges
    % 5: Chulick and Mooney (2002)
    % 6: Yan and Clayton (2007)
    % 7: Gilbert (2010)
    % 8: digitized points from SJB gocad project

    % USER INPUT: refers to type of combination of data sets (see below)
    idata = 1;
    
    stmlabs = {'crust2','EARS','salton','bryant-jones','chulick-mooney','yan-clayton','gilbert','gocad'};
    
    % 1: Crust2.0 global model
    [xlon_crust2,ylat_crust2,zdep_crust2] = read_moho_crust2;
    zdep_sigma_crust2 = zdep_crust2*0.1;
    ncrust2 = length(xlon_crust2);

%     if 0==1
%         [X,Y,Z] = griddataXB(xlon_crust2,ylat_crust2,zdep_crust2,200,'cubic');
%         figure; pcolor(X,Y,Z); shading flat; colorbar
%         xlabel('Longitude'); ylabel('Latitude');
%         title('Crust 2.0 -- Basin, Laske, Masters, 2000');
%     end

    % 2: EARS -- automated receiver functions
    [xlon_EARS,ylat_EARS,zdep_EARS,zdep_sigma_EARS] = read_moho_EARS;
    nEARS = length(xlon_EARS);

    % 3: Magistrale et al. (2000) -- for Salton trough
    [xlon_cvm4,ylat_cvm4,zdep_cvm4] = read_moho_Magistrale2000;
    isalton = find( zdep_cvm4 == 22);
    xlon_salton = xlon_cvm4(isalton);
    ylat_salton = ylat_cvm4(isalton);
    zdep_salton = zdep_cvm4(isalton);
    %zdep_sigma_salton = zdep_sigma_cvm4(isalton);
    nsalton = length(xlon_salton);

    % 4: Bryant and Jones (1992) -- for western Transverse Ranges
    [xlon_wtrb,ylat_wtrb,zdep_wtrb] = read_moho_BryantJones1992;
    nwtrb = length(xlon_wtrb);

    % 5: Chulick and Mooney (2002)
    [xlon_mooney,ylat_mooney,zdep_mooney] = read_moho_ChulickMooney2002;
    nmooney = length(xlon_mooney);

    % 6: Yan and Clayton (2007)
    [xlon_yan,ylat_yan,zdep_yan,zdep_sigma_yan] = read_moho_YanClayton2007;
    nyan = length(xlon_yan);

    % 7: Gilbert et al. (2010)
    [xlon_gilbert,ylat_gilbert,zdep_gilbert,zdep_sigma_gilbert] = read_moho_Gilbert2010;
    ngilbert = length(xlon_gilbert);

    % 8: digitized cross sections in SJB gocad project
    [xlon_gocad,ylat_gocad,zdep_gocad,zdep_sigma_gocad] = read_moho_gocad;
    ngocad = length(xlon_gocad);

    %---------------------------
    % modify/assign uncertainty estimates

    zdep_sigma_salton  = zdep_salton*0.2;
    zdep_sigma_cvm4    = zdep_cvm4*0.2;
    zdep_sigma_wtrb    = zdep_wtrb*0.2;
    zdep_sigma_mooney  = 1 + zdep_mooney*0.1;
    zdep_sigma_gilbert = 1 + zdep_sigma_gilbert;

    %---------------------------

    % KEY: specify the combination of datasets that you want
    %  1   Crust2 + Gilbert + YanClayton + ChulickMooney + 22kmSalton + gocad_digitized_profiles
    %  2   Crust2 + Gilbert
    % 10   only Crust2
    % 11   all datasets minus Gilbert (CVM6p2)

    % make one big dataset
    switch idata
        case 1
            % default choice: Crust2 + Gilbert + YanClayton + ChulickMooney
            %                        + 22kmSalton + gocad_digitized_profiles
            % This was used for created the surface for the SCEC CVM
            % (Southern California Earthquake Center Community Velocity
            % Model), using scales q=2 to q=8.
            xlon_all = [xlon_crust2 ; xlon_salton ; xlon_mooney ; xlon_yan ; xlon_gilbert ; xlon_gocad];
            ylat_all = [ylat_crust2 ; ylat_salton ; ylat_mooney ; ylat_yan ; ylat_gilbert ; ylat_gocad];
            zdep_all = [zdep_crust2 ; zdep_salton ; zdep_mooney ; zdep_yan ; zdep_gilbert ; zdep_gocad];
            zdep_sigma_all = [zdep_sigma_crust2 ; zdep_sigma_salton ; ...
                zdep_sigma_mooney ; zdep_sigma_yan ; zdep_sigma_gilbert ; zdep_sigma_gocad];
            iall = [1*ones(ncrust2,1) ; 3*ones(nsalton,1) ; ...
                5*ones(nmooney,1) ; 6*ones(nyan,1)  ; 7*ones(ngilbert,1) ; 8*ones(ngocad,1)];
        case 2
            xlon_all = [xlon_crust2 ; xlon_gilbert];
            ylat_all = [ylat_crust2 ; ylat_gilbert];
            zdep_all = [zdep_crust2 ; zdep_gilbert];
            zdep_sigma_all = [zdep_sigma_crust2 ; zdep_sigma_gilbert ];
            iall = [1*ones(ncrust2,1) ; 7*ones(ngilbert,1) ];
        case 3
            xlon_all = [xlon_crust2 ; xlon_yan];
            ylat_all = [ylat_crust2 ; ylat_yan];
            zdep_all = [zdep_crust2 ; zdep_yan];
            zdep_sigma_all = [zdep_sigma_crust2 ; zdep_sigma_yan ];
            iall = [1*ones(ncrust2,1) ; 6*ones(nyan,1) ];
        case 10
            xlon_all = [xlon_crust2];
            ylat_all = [ylat_crust2];
            zdep_all = [zdep_crust2];
            zdep_sigma_all = [zdep_sigma_crust2];
            iall = ones(ncrust2,1);
        case 11
            xlon_all = [xlon_crust2 ; xlon_EARS ; xlon_salton ; xlon_wtrb ; xlon_mooney ; xlon_yan];
            ylat_all = [ylat_crust2 ; ylat_EARS ; ylat_salton ; ylat_wtrb ; ylat_mooney ; ylat_yan];
            zdep_all = [zdep_crust2 ; zdep_EARS ; zdep_salton ; zdep_wtrb ; zdep_mooney ; zdep_yan];
            zdep_sigma_all = [zdep_sigma_crust2 ; zdep_sigma_EARS ; zdep_sigma_salton ; zdep_sigma_wtrb ; ...
                zdep_sigma_mooney ; zdep_sigma_yan ];
            iall = [1*ones(ncrust2,1) ; 2*ones(nEARS,1) ; 3*ones(nsalton,1) ; 4*ones(nwtrb,1) ; ...
                5*ones(nmooney,1) ; 6*ones(nyan,1) ];
    end

    % NOTE: increase sigma value if it is less than sigminfrac times moho depth
    sigminfrac = 0.02;
    ilowerror = find(zdep_sigma_all <= sigminfrac*zdep_all);
    nlowerror = length(ilowerror);
    disp(sprintf('increasing zdepth uncertainty for %i points :',nlowerror));
    for ii=1:nlowerror
        jj = ilowerror(ii);
        fprintf('%6i%6i%10.4f%10.4f%10.4f%10.4f%10.4f\n',ii,iall(jj),xlon_all(jj),ylat_all(jj),...
            zdep_all(jj),zdep_sigma_all(jj),sigminfrac*zdep_all(jj));
    end
    zdep_sigma_all(ilowerror) = sigminfrac*zdep_all(ilowerror);

    % remove NaN values if any exist
    iokay = find(~isnan(zdep_all));
    xlon_all       = xlon_all(iokay);
    ylat_all       = ylat_all(iokay);
    zdep_all       = zdep_all(iokay);
    zdep_sigma_all = zdep_sigma_all(iokay);
    iall           = iall(iokay);

    % quick plot showing the range of uncertainties for each dataset
    figure; plot(iall, zdep_sigma_all,'.')
    xlim([min(iall)-1 max(iall)+1]);
    xlabel('Index of dataset'); ylabel('Uncertainty in Moho depth');

    % more informative display -- histogram for each dataset
    figure; nr=4; nc=2;
    for ii=1:8
        itemp = find(iall==ii);
        subplot(nr,nc,ii);
        N = plot_histo(zdep_sigma_all(itemp),[0:1:10]);
        ylim([0 0.5]); xlabel('Uncertainty in Moho depth, km');
        title(stmlabs{ii});
    end

    %[X,Y,Z] = griddataXB(xlon_cal,ylat_cal,zdep_cal,200,'cubic');
    %figure; pcolor(X,Y,Z); shading flat; colorbar
    %xlabel('Longitude'); ylabel('Latitude'); caxis([10 60])

    dlon = xlon_all;
    dlat = ylat_all;
    d = zdep_all;
    dsig = zdep_sigma_all;
    dlabel = sprintf('moho%2.2i',idata);
    ulabel = 'zdep, km';
    
elseif dopt==2
    [dlon,dlat,zdep,zdep_sigma] = read_basement_cal;
    d = zdep;
    dsig = zdep_sigma;
    dlabel = 'USGSxtalbasement';
    ulabel = 'zdep, km';
    
elseif dopt==3
    [dlon,dlat,zdep,zdep_sigma] = read_basement_maricopa;
    d = zdep;
    dsig = zdep_sigma;
    dlabel = 'maricopa_basement';
    ulabel = 'zdep, km';
    
elseif dopt==4
    [dlon,dlat,zdep,zdep_sigma] = read_basement_sjb; 
    d = zdep;
    dsig = zdep_sigma;
    dlabel = 'base_tertiary';
    ulabel = 'zdep, km';
    
elseif dopt==5
    % isostatic residual gravity anomaly (column 5)
    [dlon,dlat,dg,gsig] = read_grav_akusgs(ax0);
    d = dg(:,5);
    dsig = gsig;
    dlabel = 'USGSgrav';
    ulabel = 'dgrav, mgal';
    
elseif dopt==6
    ifile = '/home/carltape/PROJECTS/alaska/tomo_models/alaska_moho_depth_data.dat';
    [dlon,dlat,d,dsig] = textread(ifile,'%f%f%f%f');
    dlabel = 'moho';
    ulabel = 'zdep, kml';   
end

% subset
inds = getsubset(dlon,dlat,ax0);
dlon = dlon(inds);
dlat = dlat(inds);
d = d(inds);
dsig = dsig(inds);
if dopt==1, iall = iall(inds); end

if isempty(dlon), error('get_1D_dataset.m: zero observations within specified region'); end

% label for files: geographic region, then data set
slabel = [rlabel  '_' dlabel];

ifig = 1;
if ifig==1
    if dopt==1
        figure; 
        socal_basemap(ax0,0,0,0);
        plot(dlon,dlat,'.');
        plot(xlon_crust2,ylat_crust2,'o');
        axis(ax0);
        for ii=1:length(dlon)
           text( dlon(ii), dlat(ii), sprintf('%i:%.0f',iall(ii),d(ii)));
        end
    end
    
    % plot data and uncertainties
    figure; nr=2; nc=1; msize = 6^2;
    subplot(nr,nc,1); hold on;
    scatter(dlon,dlat,msize,d,'filled');
    %scatter(dlon,dlat,msize,'ko');   % black outline
    colorbar; axis(ax0); grid on;
    xlabel('Longitude'); ylabel('Latitude');
    title(sprintf('%s observations (%i)',slabel,length(d)),'interpreter','none');
    
    subplot(nr,nc,2); hold on;
    scatter(dlon,dlat,msize,dsig,'filled');
    %scatter(dlon,dlat,msize,'ko');   % black outline
    colorbar; axis(ax0); grid on;
    xlabel('Longitude'); ylabel('Latitude');
    title(sprintf('%s uncertainties (%i)',slabel,length(d)),'interpreter','none');

end

%========================================================
