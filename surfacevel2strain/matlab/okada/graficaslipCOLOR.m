function graficaslipCOLOR(slipSol,rake,parches,ValorSlipColorMin,ValorSlipColorMax)
%function graficaslipCOLOR(slipSol,paramFoco,parches,ValorSlipColorMin,ValorSlipColorMax)
%Grafica la solucion.
% slipSol [m] : vector solucion de la inversion del slip en los parches.
% paramFoco: nombre del archivo con los parametros focales usado en la
%            inversion.
% parches : nombre del archivo con la configuracion de parches usada en la
%           inversion.

% Comprobacion
if ValorSlipColorMin==ValorSlipColorMax
    ValorSlipColorMin=ValorSlipColorMin - 0.01;
    ValorSlipColorMax=ValorSlipColorMax + 0.01;
end
%Obtencion de parametros necesarios para graficar.
FA=0.5;
P=load(parches);
%Verificacion
if min(size(slipSol))>1 
    error('slipSol debe ser un vector');
end
if length(slipSol)~=max(P(:,1))
    error('El vector slipSol debe tener la misma longitud que el numero total de parches');
end
% dibujo de los parches

Lc=zeros(max(P(:,1)),1);
Wc=zeros(max(P(:,1)),1);
Lf=zeros(max(P(:,1)),1);
Wf=zeros(max(P(:,1)),1);
for i=1:max(P(:,1)) %para la cantidad de parches.
    %Calculo del centro de cada parche para graficar el vector slip.
    Lc(i)=(P(i,2) + P(i,4))/2.0;
    Wc(i)=(P(i,3) + P(i,5))/2.0;
    %Calculo de la posicion final de cada parche.
    Lf(i)=Lc(i) + FA*slipSol(i)*cos(deg2rad(rake));
    Wf(i)=Wc(i) + FA*slipSol(i)*sin(deg2rad(rake));
end;



%Ahora con colores.

%primero debo eliminar los elementos repetidos de Lc y Wc (los almaceno en
%X para Lc e Y para Wc.
X=0;
Y=0;
Kx=0; %indice para x
Ky=0; %indice para y
for i=1:length(Lc) %Para limpiar X.
    if length(find(X==Lc(i)))==0 %si no encuentro un elemento de Lc en X lo agrego.
        Kx=Kx+1;
        X(Kx)=Lc(i);
    end
end

for i=1:length(Wc) %Para limpiar Y.
    if length(find(Y==Wc(i)))==0 %si no encuentro un elemento de Wc en Y lo agrego.
        Ky=Ky+1;
        Y(Ky)=Wc(i);
    end
end
    
% Ahora hay que ordenar X e Y con el comando SORT
X=sort(X);
Y=sort(Y);

% Ahora hay que hacer un meshgrid de X e Y
[XI,YI]=meshgrid(X,Y);

%Ahora creo la matriz Z usando el mesh grid  recorriendo las matrices X e Y
%y buscando en los vectores Lc y Wc si encuentro el par ordenado X Y, si
%esta lo agrego en la posicion correspondiente en la matriz Z y si no esta
%dejo como true una variable booleana para realizar el comando interp2 para
%rellenarla.
slipSolI = GRIDDATA(Lc,Wc,slipSol,XI,YI);
%impongo el codigo de colores en el rango deseado.
% El comando siguiente es para dejar fijo el color en los graficos con
% imagesc.
figure;
imagesc(XI(1,:),YI(:,1),slipSolI)
axis xy
set(gca,'cLimMode','manual');
set(gca,'clim',[ValorSlipColorMin ValorSlipColorMax]);
colorbar
hold on;

for i=1:max(P(:,1)) %para la cantidad de parches.
    %Ploteo de los parches
    plot([P(i,2) P(i,4) P(i,4) P(i,2) P(i,2)],[P(i,3) P(i,3) P(i,5) P(i,5) P(i,3)],'k')
    plot(Lc(i),Wc(i),'ok');
    %Grafico el slip en cada parche.
    plot([Lc(i) Lf(i)],[Wc(i) Wf(i)],'k')
end;
title('Distribucion de la dislocacion en el plano de falla [m]');
xlabel('L[km]');
ylabel('W[km]');
hold off;