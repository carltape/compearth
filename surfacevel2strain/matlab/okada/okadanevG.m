function varargout = okadanevG(varargin)
% OKADANEVG M-file for okadanevG.fig
%      OKADANEVG, by itself, creates a new OKADANEVG or raises the existing
%      singleton*.
%
%      H = OKADANEVG returns the handle to a new OKADANEVG or the handle to
%      the existing singleton*.
%
%      OKADANEVG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OKADANEVG.M with the given input arguments.
%
%      OKADANEVG('Property','Value',...) creates a new OKADANEVG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before okadanevG_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to okadanevG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help okadanevG

% Last Modified by GUIDE v2.5 29-Jun-2004 00:51:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @okadanevG_OpeningFcn, ...
                   'gui_OutputFcn',  @okadanevG_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before okadanevG is made visible.
function okadanevG_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to okadanevG (see VARARGIN)

% Choose default command line output for okadanevG
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes okadanevG wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = okadanevG_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function strikebox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function strikebox_Callback(hObject, eventdata, handles)
handles.hand_strike=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function dipbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function dipbox_Callback(hObject, eventdata, handles)
handles.hand_dip=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rakebox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function rakebox_Callback(hObject, eventdata, handles)
handles.hand_rake=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slip3Dbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function slip3Dbox_Callback(hObject, eventdata, handles)
handles.hand_slip3D=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function lbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lbox_Callback(hObject, eventdata, handles)
handles.hand_L=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function wbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function wbox_Callback(hObject, eventdata, handles)
handles.hand_W=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function profbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function profbox_Callback(hObject, eventdata, handles)
handles.hand_prof=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function lambdabox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function lambdabox_Callback(hObject, eventdata, handles)
handles.hand_lambda=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function mubox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function mubox_Callback(hObject, eventdata, handles)
handles.hand_mu=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slipelevbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function slipelevbox_Callback(hObject, eventdata, handles)
handles.hand_slipelev=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function minNbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function minNbox_Callback(hObject, eventdata, handles)
handles.hand_minN=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function maxNbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function maxNbox_Callback(hObject, eventdata, handles)
handles.hand_maxN=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pasoNbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function pasoNbox_Callback(hObject, eventdata, handles)
handles.hand_pasoN=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function minEbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function minEbox_Callback(hObject, eventdata, handles)
handles.hand_minE=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function maxEbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function maxEbox_Callback(hObject, eventdata, handles)
handles.hand_maxE=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function pasoEbox_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function pasoEbox_Callback(hObject, eventdata, handles)
handles.hand_pasoE=get(hObject,'String');
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function Az_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function Az_Callback(hObject, eventdata, handles)
handles.hand_Az=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function elevacion_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function elevacion_Callback(hObject, eventdata, handles)
handles.hand_elevacion=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nLbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nLbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

function nLbox_Callback(hObject, eventdata, handles)
handles.hand_nL=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function nWbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nWbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nWbox_Callback(hObject, eventdata, handles)
handles.hand_nW=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function mhuLbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mhuLbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mhuLbox_Callback(hObject, eventdata, handles)
handles.hand_mhuL=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function mhuWbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mhuWbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function mhuWbox_Callback(hObject, eventdata, handles)
handles.hand_mhuW=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sigmaLbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaLbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function sigmaLbox_Callback(hObject, eventdata, handles)
handles.hand_sigmaL=get(hObject,'String');
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function sigmaW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sigmaW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function sigmaW_Callback(hObject, eventdata, handles)
handles.hand_sigmaW=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function gaussbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in gaussbox.
function gaussbox_Callback(hObject, eventdata, handles)
handles.hand_gauss=get(hObject,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function nombreArchivo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nombreArchivo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function nombreArchivo_Callback(hObject, eventdata, handles)
% hObject    handle to nombreArchivo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nombreArchivo as text
%        str2double(get(hObject,'String')) returns contents of nombreArchivo as a double
handles.hand_nombreArchivo=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function datosbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datosbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


function datosbox_Callback(hObject, eventdata, handles)
% hObject    handle to datosbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of datosbox as text
%        str2double(get(hObject,'String')) returns contents of datosbox as a double
handles.hand_datosSalida=get(hObject,'String');
guidata(hObject,handles);





% --- Executes during object creation, after setting all properties.
function estacionesbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to estacionesbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function estacionesbox_Callback(hObject, eventdata, handles)
% hObject    handle to estacionesbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of estacionesbox as text
%        str2double(get(hObject,'String')) returns contents of estacionesbox as a double
handles.hand_estaciones=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Ncpbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ncpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Ncpbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ncpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ncpbox as text
%        str2double(get(hObject,'String')) returns contents of Ncpbox as a double
handles.hand_Ncp=get(hObject,'String');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Ecpbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ecpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Ecpbox_Callback(hObject, eventdata, handles)
% hObject    handle to Ecpbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ecpbox as text
%        str2double(get(hObject,'String')) returns contents of Ecpbox as a double
handles.hand_Ecp=get(hObject,'String');
guidata(hObject,handles);





% --- Executes on button press in boton_ejecutar.
function boton_ejecutar_Callback(hObject, eventdata, handles)
%Primero obtengo las variables de la estructura handles (estan como
%String's), convirtiendolas al tipo de datos necesario.

strike=str2double(handles.hand_strike);
dip=str2double(handles.hand_dip);
rake=str2double(handles.hand_rake);
slip3D=str2double(handles.hand_slip3D);
L=str2double(handles.hand_L);
W=str2double(handles.hand_W);
prof=str2double(handles.hand_prof);
slipelev=str2double(handles.hand_slipelev);
lambda=str2double(handles.hand_lambda);
mu=str2double(handles.hand_mu);
Az=str2double(handles.hand_Az);
elevacion=str2double(handles.hand_elevacion);

nL=str2double(handles.hand_nL);
nW=str2double(handles.hand_nW);
mhuL=str2double(handles.hand_mhuL);
mhuW=str2double(handles.hand_mhuW);
sigmaL=str2double(handles.hand_sigmaL);
sigmaW=str2double(handles.hand_sigmaW);
gauss=str2double(handles.hand_gauss);
nombreArchivo=handles.hand_nombreArchivo;
archivoSalida=handles.hand_datosSalida;

estaciones=handles.hand_estaciones;
Ncp=str2double(handles.hand_Ncp);
Ecp=str2double(handles.hand_Ecp);

%Cargo en memoria la matriz con las coordenadas PLANIMETRICAS de las
%estaciones. 
%formato del archivo : Numestacion  Nestacion(km)  Eestacion(km)
CoordEst=load(estaciones);
NumEst=CoordEst(:,1);
Nest=CoordEst(:,2);
Eest=CoordEst(:,3);

%Conversion de las unidades, ya que las longitudes deben ir en metros
L=L*1000;
W=W*1000;
prof=prof*1000;
Nest=Nest*1000;
Eest=Eest*1000;
Ncp=Ncp*1000;
Ecp=Ecp*1000;

% Traslacion de los puntos de observacion, para que el sistema de
% coordenadas en el cual estan descritos coincida con el necesario para el
% software (origen en la proyeccion en superficie del centro del plano de
% falla)
Nest=Nest-Ncp;
Eest=Eest-Ecp;



%Calculo los corrimientos teoricos
e=Eest;
n=Nest;
v=zeros(size(e));
un=zeros(size(e));

if(gauss==0)
    archivo='slip_constante.txt';
    rakeslipchiIguales(rake,slip3D,slipelev,nL,nW,archivo);
    [un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end
if(gauss==1) 
archivo='slip_normal.txt';
slipnormal(rake,slip3D,mhuL*L,sigmaL*L,mhuW*W,sigmaW*W,slipelev,L,nL,W,nW,archivo);
[un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end   
if(gauss==2) 
archivo=nombreArchivo;
[un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end   


%Ahora proyecto sobre L.O.S.
%ui=.33*ue+.94*uv-0.07*un;
losN=cos(deg2rad(elevacion))*cos(deg2rad(Az));
losE=cos(deg2rad(elevacion))*sin(deg2rad(Az));
losV=-sin(deg2rad(elevacion));
losE=0.37; losN=-0.09; losV=-0.91; %ANTOFAGASTA
uLOS = (losE*ue + losV*uv + losN*un);
UULOS=uLOS;
%Guardo resultados en archivo "archivoSalida" los corrimientos teoricos
%calculados en los puntos de observacion (N,E,V).
%Formato del archivo " #punto N(km) E(km) V(km) Un(m) Ue(m) Uv(m) Ulos(m) "
n=n+Ncp;
e=e+Ecp;

[n1,n2]=size(n);
archivo=fopen(archivoSalida, 'w+');
fprintf(archivo,'%% strike= %f, dip= %f, rake= %f, chi= %f(angulo de opening del slip3D) \n',strike,dip,rake,slipelev);
fprintf(archivo,'%% L[km]= %f, W[km]= %f, slip3D[m]= %f, Prof[km]= %f\n',L/1000.0,W/1000.0,slip3D,prof/1000.0);
fprintf(archivo,'%% Parametros de Lame: lambda[Pa]= %f, mhu[Pa]= %f\n',lambda,mu);
fprintf(archivo,'%%\n');



fprintf(archivo,'%% #punto N(km) E(km) V(km) Un(m) Ue(m) Uv(m) Ulos(m)\n');
fprintf(archivo,'%%\n');

for j=1:n2
    for i=1:n1
        num=(j-1)*n1 + i;
        NN=n(i,j)/1000.0;
        EE=e(i,j)/1000.0;
        VV=v(i,j)/1000.0;
        UUN=un(i,j);
        UUE=ue(i,j);
        UUV=uv(i,j);
        UULOS=uLOS(i,j);
        fprintf(archivo,'%f %f %f %f %f %f %f %f \n',num,NN,EE,VV,UUN,UUE,UUV,UULOS);
              
    end
end
fclose(archivo);

ccc=['save ' archivoSalida '.mat n e v un ue uv UULOS -MAT'];
eval(ccc);


%Ahora para graficar los campos de corrimientos teoricos.
minE=min(Eest);
maxE=max(Eest);

minN=min(Nest) - 0.2*abs(max(Nest)-min(Nest));
maxN=max(Nest) + 0.2*abs(max(Nest)-min(Nest));

pasoE=abs(maxE-minE)/400;
pasoN=abs(maxN-minN)/400;

[e,n]=meshgrid(minE:pasoE:maxE,minN:pasoN:maxN);
v=zeros(size(e));
un=zeros(size(e));

if(gauss==0)
    archivo='slip_constante.txt';
    rakeslipchiIguales(rake,slip3D,slipelev,nL,nW,archivo);
    [un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end
if(gauss==1) 
archivo='slip_normal.txt';
slipnormal(rake,slip3D,mhuL*L,sigmaL*L,mhuW*W,sigmaW*W,slipelev,L,nL,W,nW,archivo);
[un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end   
if(gauss==2) 
archivo=nombreArchivo;
[un,ue,uv,ux,uy,uz,estado]=okadaparchesNEV(lambda,mu,L,nL,W,nW,prof,strike,dip,archivo,e,n,v);
end   
%%%%%%%%% Defino coordenadas del plano de falla antes del deslizamiento para graficarlo.
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



GrosorLinea=1.3;
figure
subplot(2,2,1); imagesc(e(1,:)/1000,n(:,1)/1000,un);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "un"'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
plot(Eest/1000.0,Nest/1000.0,'+r');
hold off;
subplot(2,2,2);,imagesc(e(1,:)/1000,n(:,1)/1000,ue);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "ue"'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
plot(Eest/1000.0,Nest/1000.0,'+r');
hold off;
subplot(2,2,3);imagesc(e(1,:)/1000,n(:,1)/1000,-uv);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "uv" +=up'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
plot(Eest/1000.0,Nest/1000.0,'+r');
hold off;

%Ahora proyecto sobre L.O.S.
%ui=.33*ue+.94*uv-0.07*un;
%losN=cos(deg2rad(elevacion))*cos(deg2rad(Az));
%losE=cos(deg2rad(elevacion))*sin(deg2rad(Az));
%losV=-sin(deg2rad(elevacion));

%uLOS = (losE*ue + losV*uv + losN*un);

subplot(2,2,4);imagesc(e(1,:)/1000,n(:,1)/1000,mod(uLOS,0.028));colorbar;axis equal;axis xy;
title('L.O.S.'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
plot(Eest/1000.0,Nest/1000.0,'+r');
hold off;



figure
subplot(2,2,1); imagesc(e(1,:)/1000,n(:,1)/1000,un);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "un"'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;
subplot(2,2,2);,imagesc(e(1,:)/1000,n(:,1)/1000,ue);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "ue"'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;
subplot(2,2,3);imagesc(e(1,:)/1000,n(:,1)/1000,-uv);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "uv" +=up'); xlabel('E(km)');ylabel('N (km)');
hold on;
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;

%Ahora proyecto sobre L.O.S.
%ui=.33*ue+.94*uv-0.07*un;
losN=cos(deg2rad(elevacion))*cos(deg2rad(Az));
losE=cos(deg2rad(elevacion))*sin(deg2rad(Az));
losV=-sin(deg2rad(elevacion));

uLOS = (losE*ue + losV*uv + losN*un);

subplot(2,2,4);imagesc(e(1,:)/1000,n(:,1)/1000,mod(uLOS,0.028));colorbar;axis equal;axis xy;
title('L.O.S.'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;

figure;imagesc(e(1,:)/1000,n(:,1)/1000,uLOS);colorbar;axis equal;axis xy;
title('L.O.S. Desarrollado'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;


figure
% Aqui se grafica 
subplot(2,2,1); imagesc(e(1,:)/1000,n(:,1)/1000,ux);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "ux"'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;

subplot(2,2,2);,imagesc(e(1,:)/1000,n(:,1)/1000,uy);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "uy"'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;
subplot(2,2,3);imagesc(e(1,:)/1000,n(:,1)/1000,uz);colorbar;axis equal;axis xy;
title('Campo de desplazamientos "uz" +=up'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;

%Ahora proyecto sobre L.O.S.
%ui=.33*ue+.94*uv-0.07*un;
losX=cos(deg2rad(elevacion))*cos(deg2rad(Az-strike));
losY=-cos(deg2rad(elevacion))*sin(deg2rad(Az-strike));
losZ=sin(deg2rad(elevacion));


uLOSxyz = (losX*ux + losZ*uz + losY*uy);

subplot(2,2,4);imagesc(e(1,:)/1000,n(:,1)/1000,mod(uLOSxyz,.028));colorbar;axis equal;axis xy;
title('L.O.S.'); xlabel('E(km)');ylabel('N (km)');
hold on; 
plot(planoE,planoN,'k-','LineWidth',GrosorLinea);
plot(slipE,slipN,'m.-');
plot(ejesE,ejesN,'k-','LineWidth',GrosorLinea);
hold off;

figure;
disp_fig(n/1000,e/1000,v/1000,0,0,0,un,ue,uv,1,1,1,200,'N','E','V','OKADANEVG');



