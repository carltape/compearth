function [e,n] = creaPuntosObs(minE,pasoE,maxE,minN,pasoN,maxN,archivoSALIDA)
%function [e,n] = creaPuntosObs(minE,pasoE,maxE,minN,pasoN,maxN,archivoSALIDA)
% Todas las unidades en (km)
% archivoSALIDA es un string.
[e,n]=meshgrid(minE:pasoE:maxE,minN:pasoN:maxN);
[n1,n2]=size(n);
archivo=fopen(archivoSALIDA, 'w+');
fprintf(archivo,'%% NumeroEstacion, Nestacion(km)  Eestacion(km)\n');
fprintf(archivo,'%%\n');

for j=1:n2
    for i=1:n1
        progreso(i + (j-1)*n1,n1*n2);
        num=(j-1)*n1 + i;
        NN=n(i,j);
        EE=e(i,j);
        fprintf(archivo,'%f %f %f \n',num,NN,EE);
              
    end
end
fclose(archivo);

