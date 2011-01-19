function [un,ue,uv,ux,uy,uz,estado]=okadapanchonev(lambda,mu,L,W,profPL,strike,dip,rake,slip3D,chi,e,n,v)
%function [un,ue,uv,ux,uy,uz,estado]=okadapanchonev(lambda,mu,L,W,profPL,strike,dip,rake,slip3D,chi,e,n,v)
%
% * lambda, mu : constantes de LAME del medio (Pa).
% * L: Largo de la superficie de ruptura (m).
% * W: Ancho de la superficie de ruptura (m).
% * profPL: profundidad del centro de la superficie de falla (m).
% * strike, dip, rake (Angulos en grados sexagesimales).
% * slip3D: dislocacion en [m] medida en el eje que define el vector de
%           movimiento 3D (considerando U1, U2 y U3 de okada).
% * chi: angulo que forma el vector de movimiento con su proyeccion en el
%        plano de falla.
% * e,n,v: coordenadas, en el sistema geografico, de los puntos en los cuales 
%          se quiere calcular los corrimientos. (m).
% * estado : numero de puntos singulares encontrados durante el calculo. 

if nargin ~= 13
    msgbox('EL NUMERO DE ARGUMENTOS DEBE SER 13');
    error('EL NUMERO DE ARGUMENTOS DEBE SER 13');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Determinacion de los parametros del medio %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% parametro del medio para Okada %%%%%
alpha=(lambda+mu)/(lambda+2*mu);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Determinacion de los parametros de la fuente de Okada %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U1=slip3D*cos(deg2rad(chi))*cos(deg2rad(rake));
U2=slip3D*cos(deg2rad(chi))*sin(deg2rad(rake));
U3=slip3D*sin(deg2rad(chi));
al1=0;
al2=L;
aw1=0;
aw2=W;
prof = profPL + W*sin(deg2rad(dip));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Transformacion de coordenadas de los puntos en que se mediran los corrimientos %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Verificacion Preliminar %%%%%
saleplano=prof - (W)*sin(deg2rad(dip));
%if saleplano <= 0 
    %msgbox('El plano de falla se sale de la superficie libre: modificar W o PROF');
    %error('El plano de falla se sale de la superficie libre: modificar W o PROF');
%end;  
input_size=size(n);
if prod(size(e)) ~= prod(input_size) | prod(size(v)) ~= prod(input_size);
	error('n,e y v deben tener iguales dimensiones!');
end;

%%%%% Cambio de base: sistema geografico ---> sistema de Okada %%%%%
x= cos(deg2rad(strike))*n + sin(deg2rad(strike))*e+ L/2.0;;
y= sin(deg2rad(strike))*n - cos(deg2rad(strike))*e + (W/2.0)*cos(deg2rad(dip));
z=-v;

%%%%% Definicion de las dimensiones de los datos de salida %%%%%
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
j=0;
for i=1:prod(input_size);
	[ux(i),uy(i),uz(i),uxx(i),uyx(i),uzx(i),uxy(i),uyy(i),uzy(i),uxz(i),uyz(i),uzz(i),iret]=dc3d(alpha,x(i),y(i),z(i),prof,dip,al1,al2,aw1,aw2,U1,U2,U3);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Generacion de la salida %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Cambio de sistema de coordenadas de Okada a %%%%%
%%%%%    Geografico para los desplazamientos      %%%%%



un= cos(deg2rad(strike))*ux + sin(deg2rad(strike))*uy;
ue= sin(deg2rad(strike))*ux - cos(deg2rad(strike))*uy;
uv=-uz;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Vaciando de la memoria las variables intermedias %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
