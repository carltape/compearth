%
% plot_wavelets.m
% Pablo Muse
% 
% This outputs a set of figures for our 2009 GPS paper.
%

close all, clear
clc, format short, format compact

%------------

%scales = [2 3 4 5 6 7]; %scale
scales = [0:8]; %scale
N = length(scales);

if 0==1
    figure;
    [phi,theta] = sphgrid(512);
    for i = 1:N

        j = scales(i);

        subplot(2,3,i);
        aj = 1/2^j;
        maxffj = (1-1/(1.25)^2)/aj;
        
        ff_j = dogsph(aj,theta,phi);  %./ maxffj;
        
        yashow(ff_j,'spheric','relief');
        %yashow(ff_j,'spheric');
        
        grid off; axis off
        %axis([-1 1 -1 1 -1 1])
        %axis manual
        camzoom(2.0)
        colormap([.8 .8 .8])
        view(109,34)
        light('Position',[0 1 1],'Style','infinite');
      %  text(-2.3,-1.6,['Scale a = 1/2^',num2str(j)],'FontSize',10)

    end
    
    break
end


%-----------------------------------------------------------

icompute = 0;

if icompute == 1

    % KEY: determines the sampling of the curves (512 default)
    Nx = 1024;
    
    A1 = cell(1,N); Q1 = cell(1,N);
    [phi,theta] = sphgrid(2*Nx);
    vecthetas = theta(:,1);
    vecff_all = zeros(2*Nx,N);
    disp('computing spatial wavelets...');
    for i = 1:N
        i
        j = scales(i);
        ff = dogsph(1/2^j,theta,phi);
        vecff = ff(:,1);
        vecff_all(:,i) = vecff;
        A1{i} = ['a = 2^{-',num2str(j),'}'];
        Q1{i} = ['q = ' num2str(j) ];
    end
    
    A2 = cell(1,N); Q2 = cell(1,N);
    fmat_plot = zeros(Nx,N);
    disp('computing spectral wavelets...');
    for i = 1:N
        i
        j = scales(i);
        ff = dogsph(1/2^j,theta,phi);
        fmat = fst(ff);
        fmat_plot(:,i) = abs(fmat(1,:)).^2;   % NOTE: squared value
        A2{i} = ['a = 2^{-',num2str(j),'}'];
        Q2{i} = ['q = ' num2str(j) ];
    end
    
    save('wavelets_plots','scales','vecff_all','fmat_plot','theta','phi',...
        'vecthetas','A1','A2','Q1','Q2');
else
    load('wavelets_plots');
end

% limits for the axes
% reduce the number of cross-sections for plotting
%imin = 3; imax = 6; degmax = 25; Lmax = 100;
imin = 1; imax = 2; degmax = 180; Lmax = 20;

% reduce the number of cross-sections for plotting
vecff_all_red = vecff_all(:,imin:imax);
fmat_plot_red = fmat_plot(:,imin:imax);

figure; nr=2; nc=1;

% axes limits
zstretch = 0.1;
zmax0 = max(vecff_all_red(:)); zmin0 = min(vecff_all_red(:));
zrange = zmax0 - zmin0;
zmax = zmax0 + zstretch/2*zrange;
zmin = zmin0 - zstretch/2*zrange;

subplot(nr,nc,1); hold on;
plot(vecthetas*180/pi,vecff_all_red','LineWidth',2)
%for i = 1:N
%    vecff = vecff_all(:,i);
%    plot(vecthetas*180/pi,vecff,'LineWidth',2)
%end
legend(Q1(imin:imax));
axis([0 degmax zmin zmax])
xlabel('Angular distance from North Pole, degrees');
ylabel('Amplitude');
plot(vecthetas*180/pi,zeros(size(vecthetas)),'k:')

%-----------------------------------

% axes limits
zstretch = 0.1;
zmax0 = max(fmat_plot_red(:));
zmin0 = 0;
zrange = zmax0 - zmin0;
zmax = zmax0 + zstretch/2*zrange;
zmin = zmin0 - zstretch/2*zrange;

subplot(nr,nc,2); hold on;
plot(fmat_plot_red,'LineWidth',2);
%plot(fmat_plot_red,'LineWidth',2);
%for i = 1:N
%    fmat = fmat_all(:,:,i);
%    plot(abs(fmat(1,:)),'LineWidth',2)
%end
legend(Q2(imin:imax));
axis([0 Lmax zmin zmax])
xlabel('Index of spherical harmonic coefficient');
ylabel('Spectral power');

%===================================================================
