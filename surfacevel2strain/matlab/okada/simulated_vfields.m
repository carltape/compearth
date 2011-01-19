dir_gps     = '/home/muse/gmt/gps_data/';
[dlon,dlat,ve,vn,se,sn,ren,name] = ...
    textread([dir_gps 'socal/socal_vfield_3p0_psvelo.dat'],...
    '%f%f%f%f%f%f%f%s');

utm_zone = 11;
[xrec,yrec,I_ZONE]=utm2ll(dlon,dlat,utm_zone,1);

% change scales to km
xrec = xrec*0.001;
yrec = yrec*0.001;

% create grid on the observation domain
max_xrec = max(xrec); min_xrec = min(xrec);
max_yrec = max(yrec); min_yrec = min(yrec);

stepx = (max_xrec - min_xrec)/50;
stepy = (max_yrec - min_yrec)/50;

[e, n] = creaPuntosObs(min_xrec,stepx,max_xrec - stepx,...
    min_yrec,stepy,max_yrec-stepy,'sampling_domain.txt');
    
xcenter_zone = 0.5*(min_xrec + max_xrec)
ycenter_zone = 0.5*(min_yrec + max_yrec)

% create file of stations for okada code

archivo=fopen('stations_coord.txt', 'w+');
fprintf(archivo,'%% Coordenadas de estaciones socal \n');
fprintf(archivo,'%% NumeroEstacion, Nestacion(km)  Eestacion(km)\n');
fprintf(archivo,'%%\n');

nstations = length(xrec);

for j=1:nstations
        fprintf(archivo,'%f %f %f \n',j,yrec(j),xrec(j));
end
fclose(archivo);