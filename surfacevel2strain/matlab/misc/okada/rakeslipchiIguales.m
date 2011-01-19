function rakeslipchiIguales(rake,slip3D,chi,nL,nW,nombrearchivo);
%function rakeslipchiIguales(rake,slip3D,chi,nL,nW,archivo);
% Crea el archivo necesario para okadaparchesNEV, con los parametros
% ingresados iguales para todos los parches.
archivo=fopen(nombrearchivo, 'w+');
fprintf(archivo,'%% rake slip3D chi Li Wj\n');

for j=1:nW
    
    for i=1:nL
        fprintf(archivo,'%d\t',rake);
        fprintf(archivo,'%d\t',slip3D);
        fprintf(archivo,'%d\t',chi);
        fprintf(archivo,'%i\t %i\n',i,j);
    end
end
st=fclose(archivo); 