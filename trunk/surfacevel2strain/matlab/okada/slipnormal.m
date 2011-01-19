function slipnormal(rake,slip3D,muX,sigmaX,muY,sigmaY,chi,L,nL,W,nW,nombrearchivo);
%function slipnormal(rake,slip3D,muX,sigmaX,muY,sigmaY,chi,L,nL,W,nW,nombrearchivo);
% Crea el archivo necesario para okadaparchesNEV, con los parametros
% ingresados iguales para todos los parches, salvo la magnitud del
% desplazamiento que sigue una tendencia normal en x e y (se debe
% especificar las esperanzas muX, muY y las varianzas sigmaX y sigma Y.
% Sistema de coordenadas X e Y de Okada.
% El vector slip3d se calcula en el centro de cada parche.
archivo=fopen(nombrearchivo, 'w+');
fprintf(archivo,'%% rake slip3D chi Li Wj\n');

for j=1:nW
    for i=1:nL
        u1(i)=(i - 0.5)*L/nL;
        u2(j)=(j - 0.5)*W/nW; %Coordenadas del centro de cada parche.
        slip(i,j)=slip3D*sigmaX*sigmaY*2*pi*normpdf(u1(i),muX,sigmaX)*normpdf(u2(j),muY,sigmaY);
    end
end
promedio_slip=0;
for j=1:nW
    
    for i=1:nL
        promedio_slip=promedio_slip+slip(i,j)/nW/nL;
    end
end
slip=slip*slip3D/promedio_slip; %Desplazamiento promedio=slip3D.
for j=1:nW
    
    for i=1:nL
        fprintf(archivo,'%d\t',rake);
        fprintf(archivo,'%d\t',slip(i,j));
        fprintf(archivo,'%d\t',chi);
        fprintf(archivo,'%i\t %i\n',i,j);
    end
end

st=fclose(archivo);
figure;
imagesc(u2/1000,u1/1000,slip);colorbar;axis equal;axis xy;
xlabel('W [km]');ylabel('L [km]');title('Distribucion de la dislocacion en el plano de falla (m)');