% function plot_okada_fields(E_stations,N_stations,lambda,mu,L,nL,W,nW,Z,strike,dip,rake,slip3D,slipelev)
%
% plots a very dense (almost continuous) estimation of the displacements
% using okada model.
%
% Called by socal_gps_syn.m
%
% Pablo Muse, 10/23/07 (adapted from okadanevG.m by Francisco Ortega)
%
function plot_okada_fields(E_stations,N_stations,lambda,mu,L,nL,W,nW,Z,strike,dip,rake,slip3D,slipelev)


minE=min(E_stations);
maxE=max(E_stations);

delta = 0.2*abs(max(N_stations)-min(N_stations));
minN = min(N_stations) - delta;
maxN = max(N_stations) + delta;

pasoE=abs(maxE-minE)/400;
pasoN=abs(maxN-minN)/400;

[e,n]=meshgrid(minE:pasoE:maxE,minN:pasoN:maxN);
v=zeros(size(e));
un=zeros(size(e));


archivo='slip_constante.txt';
rakeslipchiIguales(rake,slip3D,slipelev,nL,nW,archivo);
[un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,Z,strike,dip,archivo,e,n,v);

%%%%%%%%% Define coordinates of the plane fault before displacement.
planoX(1)=-L/2; planoY(1)=-(W/2)*cos(deg2rad(dip));
planoX(2)=-L/2; planoY(2)=(W/2)*cos(deg2rad(dip)); 
planoX(3)=L/2;  planoY(3)=(W/2)*cos(deg2rad(dip)); 
planoX(4)=L/2;  planoY(4)=-(W/2)*cos(deg2rad(dip));
planoX(5)=-L/2; planoY(5)=-(W/2)*cos(deg2rad(dip));
% del 6 al 10 es para graficar el plano en tama�o real, dando una
% referencia de la magnitud del dip.
planoX(6)=-L/2; planoY(6)=(W/2); 
planoX(7)=L/2;  planoY(7)=(W/2); 
planoX(8)=L/2;  planoY(8)=-(W/2);
planoX(9)=-L/2;  planoY(9)=-(W/2);
planoX(10)=-L/2; planoY(10)=-(W/2)*cos(deg2rad(dip));
%rotacion de las coordenadas del plano.divido por 1000 para pasar a km.
planoN=(cos(deg2rad(strike))*planoX + sin(deg2rad(strike))*planoY)/1000.0;
planoE=(sin(deg2rad(strike))*planoX - cos(deg2rad(strike))*planoY)/1000.0;
% Ahora defino las coordenadas del vector slip2D para graficarlo
slip2D=slip3D*cos(deg2rad(slipelev));
slipX(1)=0; slipY(1)=0;
slipX(2)=slip2D*cos(deg2rad(rake)); slipY(2)=slip2D*sin(deg2rad(rake));
%rotacion
slipN=cos(deg2rad(strike))*slipX + sin(deg2rad(strike))*slipY;
slipE=sin(deg2rad(strike))*slipX - cos(deg2rad(strike))*slipY;
%Ahora defino los ejes del sistema de coordenadas de okadanevG XYZ, para
%graficarlos.
ejeX(1)=1.75*L/1000.0; ejeY(1)=-(W/2000.0)*cos(deg2rad(dip));
ejeX(2)=-L/2000.0;             ejeY(2)=-(W/2000.0)*cos(deg2rad(dip));
ejeX(3)=-L/2000.0;             ejeY(3)=1.1*L/1000.0;
%rotacion
ejesN=cos(deg2rad(strike))*ejeX + sin(deg2rad(strike))*ejeY;
ejesE=sin(deg2rad(strike))*ejeX - cos(deg2rad(strike))*ejeY;



LineWidth=1.3;
figure
subplot(1,3,1); imagesc(e(1,:)/1000,n(:,1)/1000,un);colorbar;axis xy;%axis equal;
title('Velocity field: v_n'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',LineWidth);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',LineWidth);
plot(E_stations/1000.0,N_stations/1000.0,'+r');
hold off;
subplot(1,3,2); imagesc(e(1,:)/1000,n(:,1)/1000,ue);colorbar;axis xy;%axis equal;
title('Velocity field: v_e'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',LineWidth);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',LineWidth);
plot(E_stations/1000.0,N_stations/1000.0,'+r');
hold off;
subplot(1,3,3);imagesc(e(1,:)/1000,n(:,1)/1000,-uv);colorbar;axis xy;%axis equal;
title('Velocity field: v_u +=up'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',LineWidth);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',LineWidth);
plot(E_stations/1000.0,N_stations/1000.0,'+r');
hold off;
