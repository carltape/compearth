function parchesH(L,nL,W,nW,nombrearchivo)
%function parchesH(L,nL,W,nW,nombrearchivo)
%
% parchesH crea un archivo de nombre "nombrearchivo" de configuracion de los parches necesario para
% la inversion usando el programa GenInversoP.
%
% L [km] = Largo del plano de falla.
% W [km] = Ancho del plano de falla.
% nL = numero de parches en la direccion L.
% nW = numero de parches en la direccion W.
% nombrearchivo = string con el nombre del archivo en que se desea guardar
% la configuracion de los parches, Ej: 'parches.dat'.
% 
% La numeracion de los parches sera como se ilustra en el esquema siguiente
%
% |---------|----------|---------|---------|---------|---------|
% |         |          |         |         |         |         |
% |    13   |     14   |    15   |    16   |   17    |    18   |
% |---------|----------|---------|---------|---------|---------|
% |         |          |         |         |         |         |
% |    7    |     8    |    9    |    10   |   11    |    12   |
% |---------|----------|---------|---------|---------|---------|
% |         |          |         |         |         |         |
% |    1    |     2    |    3    |    4    |   5     |    6    |
% |---------|----------|---------|---------|---------|---------|
%
% |<---------------------------- L --------------------------->|
% 
% En el archivo de salida hay 4 columnas, la primera es i (numero del
% parche),la segunda AL1, la tercera AW1, la cuarta AL2 y la quinta AW2
% (Okada Y. 1992).

archivo=fopen(nombrearchivo, 'w+');
fprintf(archivo,'%% k i j AL1i AW1i AL2i AW2i\n');

for j=1:nW
    
    for i=1:nL
        k=(j-1)*nL + i;
        AL1=(i-1)*L/nL;
        AW1=(j-1)*W/nW;
        AL2=i*L/nL;
        AW2=j*W/nW;

        fprintf(archivo,'%d\t %d\t %d\t %d\t %d\t %d\t %d\n',k,i,j,AL1,AW1,AL2,AW2);
    end
end
st=fclose(archivo);