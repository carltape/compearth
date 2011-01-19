function [unsum,uesum,uvsum,estado]=okadaparchesNEVp(lambda,mu,PARCHES,profPL,strike,dip,rake,opening,slipParches,e,n,v)
%function [unsum,uesum,uvsum,uxsum,uysum,uzsum,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,profPL,strike,dip,rakeslipchiFILE,e,n,v)
%
% * e,n,v: coordenadas, en el sistema geografico local(cuyo origen es la proyeccion en superficie del plano de falla) de los puntos en los cuales 
%          se quiere calcular los corrimientos. (m).
% * estado : numero de puntos singulares encontrados durante el calculo. 

if nargin ~= 12
    msgbox('EL NUMERO DE ARGUMENTOS DEBE SER 12');
    error('EL NUMERO DE ARGUMENTOS DEBE SER 12');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Determinacion de los parametros del medio %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parametro del medio para Okada %%%%%
alpha=(lambda+mu)/(lambda+2*mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Determinacion de los parametros de la fuente de Okada %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(slipParches)
    
U1(k)=slipParches(k)*cos(deg2rad(opening))*cos(deg2rad(rake));
U2(k)=slipParches(k)*cos(deg2rad(opening))*sin(deg2rad(rake));
U3(k)=slipParches(k)*sin(deg2rad(opening));

end

W=max(PARCHES(:,7)) - min(PARCHES(:,5));
L=max(PARCHES(:,6)) - min(PARCHES(:,4));
prof = profPL + W*sin(deg2rad(dip));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Transformacion de coordenadas de los puntos en que se mediran los corrimientos %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verificacion Preliminar %%%%%
input_size=size(n);
if prod(size(e)) ~= prod(input_size) | prod(size(v)) ~= prod(input_size);
	error('n,e y v deben tener iguales dimensiones!');
end;

%%%%% Cambio de base: sistema geografico ---> sistema de Okada %%%%%
x= cos(deg2rad(strike))*n + sin(deg2rad(strike))*e+ L/2.0;;
y= sin(deg2rad(strike))*n - cos(deg2rad(strike))*e + (W/2.0)*cos(deg2rad(dip));
z=-v;

%%%%% Definicion de las dimensiones de los datos de salida %%%%%
uxsum=zeros(input_size);
uysum=zeros(input_size);
uzsum=zeros(input_size);
uxxsum=zeros(input_size);
uyxsum=zeros(input_size);
uzxsum=zeros(input_size);
uxysum=zeros(input_size);
uyysum=zeros(input_size);
uzysum=zeros(input_size);
uxzsum=zeros(input_size);
uyzsum=zeros(input_size);
uzzsum=zeros(input_size);
ux=zeros(input_size);
uy=zeros(input_size);
uz=zeros(input_size);
uxx=zeros(input_size);
uyx=zeros(input_size);
uzx=zeros(input_size);
uxy=zeros(input_size);
uyy=zeros(input_size);
uzy=zeros(input_size);
uxz=zeros(input_size);
uyz=zeros(input_size);
uzz=zeros(input_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Llamada a la rutina MEX Okada %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:length(slipParches) % for para recorrer los parches.
        al1=PARCHES(k,4);
        al2=PARCHES(k,6);
        aw1=PARCHES(k,5);
        aw2=PARCHES(k,7);     
        j=0;
for i=1:prod(input_size);
	[ux(i),uy(i),uz(i),uxx(i),uyx(i),uzx(i),uxy(i),uyy(i),uzy(i),uxz(i),uyz(i),uzz(i),iret]=dc3d(alpha,x(i),y(i),z(i),prof,dip,al1,al2,aw1,aw2,U1(k),U2(k),U3(k));
	if iret == 1;
		j=j+1;
		if nargout ~= 30;
			disp(sprintf('Point asked at geographical coordinates (%g,%g,%g) and Okada coordinates (%g,%g,%g) is singular...',n(i),e(i),v(i),x(i),y(i),z(i)));
		end;
		indexes(j)=i;
    end;
end;
estado=j;
if j ~= 0;
	if nargout == 30;
		disp(sprintf('okada92f found %g singular points among the points where it is supposed to work',j));
		indexes=indexes';
	else
		disp(sprintf('okada92f removed %g singular points among the points where it is supposed to work',j));
		i=1:prod(input_size);
		for k=1:j;
			int_indexes=find(i ~= indexes(k));
			i=i(int_indexes);
		end;
		x=x(i);
		y=y(i);
		z=z(i);
		ux=ux(i);
		uy=uy(i);
		uz=uz(i);
		uxx=uxx(i);
		uyx=uyx(i);
		uzx=uzx(i);
		uxy=uxy(i);
		uyy=uyy(i);
		uzy=uzy(i);
		uxz=uxz(i);
		uyz=uyz(i);
		uzz=uzz(i);
		clear indexes int_indexes k;
	end;
end;
uxsum=uxsum + ux;
uysum=uysum + uy;
uzsum=uzsum + uz;

end %for k

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generacion de la salida %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cambio de sistema de coordenadas de Okada a %%%%%
%%%%%    Geografico para los desplazamientos      %%%%%



unsum= cos(deg2rad(strike)).*uxsum + sin(deg2rad(strike)).*uysum;
uesum= sin(deg2rad(strike)).*uxsum - cos(deg2rad(strike)).*uysum;
uvsum=-uzsum;
