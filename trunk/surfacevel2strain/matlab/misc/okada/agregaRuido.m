function agregaRuidoArchivo(ArchivoObs,ruidoN,ruidoE,ruidoV,archivoSalida)
%function agregaRuidoArchivo(ArchivoObs,ruidoN,ruidoE,ruidoV,archivoSalida)

disp('Leyendo Archivo de datos...');

datos=load(ArchivoObs);
Un=datos(:,5);
Ue=datos(:,6);
Uv=datos(:,7);

disp('Añadiendo ruido a los datos...');
%se agrega ruido a Un
rand('state',sum(100*clock));
UnRuido=Un + ruidoN*(2*rand(size(Un)) - 1);

%se agrega ruido a Ue
rand('state',sum(500*clock));
UeRuido=Ue + ruidoE*(2*rand(size(Ue)) - 1);


%se agrega ruido a Uv
rand('state',sum(256*clock));
UvRuido=Uv + ruidoV*(2*rand(size(Uv)) - 1);

disp(['Guardando datos con ruido en archivo : ' archivoSalida]);

datos(:,5)=UnRuido;
datos(:,6)=UeRuido;
datos(:,7)=UvRuido;

c=['save ' archivoSalida ' datos -ASCII'];
eval(c);


disp('Proceso terminado...');



