function disp_fig(n,e,v,n_step,e_step,v_step,un,ue,uv,factn,facte,factv,fringe_length,ncomment,ecomment,vcomment,titlecomment)
% DISP_FIG, a matlab utility function to display the displacements generated by okada92f.
%
%	Usage : disp_fig(n,e,v,n_step,e_step,v_step,un,ue,uv,factn,facte,factv,fringe_length,ncomment,ecomment,vcomment,
%		titlecomment), where:
%               - n = abscissa (North direction) of the grid where the displacements have to be represented (m) ;
%               - e = ordinate (East direction) of the grid where the displacements have to be represented (m) ;
%               - v = 3rd coordinate (descending vertical direction) of the grid where the displacements have to be
%                     represented (m) ;
%               - n_step = step in the North direction of the grid where the displacements have to be represented (m) ;
%               - e_step = step in the East direction of the grid where the displacements have to be represented (m) ;
%               - v_step = step in the descending vertical direction of the grid where the displacements have to be
%                          represented (m) ;
%               - un = displacement in the North direction (m) ;
%               - ue = displacement in the East direction (m) ;
%               - uv = displacement in descending vertical direction (m) ;
%		- factn = possible amplification to be applied to the North displacements ;
%               - facte = possible amplification to be applied to the East displacements ;
%               - factv = possible amplification to be applied to the descending vertical displacements ;
%		- fringe_length = length of the fringe for the vertical displacements (in m) ;
%		- ncomment = comment to be applied to the North axis ;
%		- ecomment = comment to be applied to the East axis ;
%		- vcomment = comment to be applied to the vertical axis ;
%		- titlecomment = general comment to be applied to the figure.
%
%       Version : 1.0b first written 11/06/1997, revised 04/02/1998
%
%       Author : Patrick Pinettes (pinettes@ipgp.jussieu.fr), IPGP Sismologie.

if n_step ~= 0;
        nn=round((max(n)-min(n))/n_step+1);
else
        nn=1;
end;
clear n_step;
if e_step ~= 0;
        ne=round((max(e)-min(e))/e_step+1);
else
        ne=1;
end;
clear e_step;
if isnan(v_step);
        if v_step ~= 0;
                nv=round((max(v)-min(v))/v_step+1);
        else
                nv=1;
        end;
else
        nv=1;
end;
clear v_step;
n=n+un*factn;
e=e+ue*facte;
v=v+uv*factv;
clear un ue uv factn facte;
max_color=max(max(v));
min_color=min(min(v));
%n=reshape(n,nn,ne,nv);
%e=reshape(e,nn,ne,nv);
%v=reshape(v,nn,ne,nv);
clear nn ne nv;
delta=max_color-min_color;
max_color=max_color+0.1*delta;
min_color=min_color-0.1*delta;
clear delta;
fringe_length=fringe_length*factv;
clear factv;
nb=round(16*((max_color-min_color)/fringe_length+1));
clear fringe_length;
colormap(sar(nb));
clear nb;
Display=uimenu(gcf,'Label','Display');
ColorMap=uimenu(Display,'Label','Color map');
uimenu(ColorMap,'Label','SAR', 'Callback','colormap(sar)');
uimenu(ColorMap,'Label','Rainbow', 'Callback','colormap(rainbow)');
uimenu(ColorMap,'Label','Saturated', 'Callback','colormap(hsv)');
uimenu(ColorMap,'Label','Jet','Callback','colormap(jet)');
uimenu(ColorMap,'Label','Hot','Callback','colormap(hot)');
uimenu(ColorMap,'Label','Cool','Callback','colormap(cool)');
uimenu(ColorMap,'Label','Spring','Callback','colormap(spring)');
uimenu(ColorMap,'Label','Summer','Callback','colormap(summer)');
uimenu(ColorMap,'Label','Autumn','Callback','colormap(autumn)');
uimenu(ColorMap,'Label','Winter','Callback','colormap(winter)');
uimenu(ColorMap,'Label','Bone','Callback','colormap(bone)');
uimenu(ColorMap,'Label','Pinks','Callback','colormap(pink)');
uimenu(ColorMap,'Label','Copper','Callback','colormap(copper)');
uimenu(ColorMap,'Label','Prism','Callback','colormap(prism)');
uimenu(ColorMap,'Label','White','Callback','colormap(white)');
Shading=uimenu(Display,'Label','Shading');
uimenu(Shading,'Label','Mesh', ...
	'Callback','axis_children=get(gca,''Children'');for i=1:length(axis_children);if strcmp(get(axis_children(i),''Type''),''surface'') == 1;set(axis_children(i),''EdgeColor'',''flat'',''FaceColor'',[0 0 0]);end;end;clear i;clear axis_children;');
uimenu(Shading,'Label','Surf', ...
	'Callback','axis_children=get(gca,''Children'');for i=1:length(axis_children);if strcmp(get(axis_children(i),''Type''),''surface'') == 1;set(axis_children(i),''EdgeColor'',[0 0 0],''FaceColor'',''flat'');end;end;clear i;clear axis_children;');
uimenu(Shading,'Label','Smoothed', ...
	'Callback','axis_children=get(gca,''Children'');for i=1:length(axis_children);if strcmp(get(axis_children(i),''Type''),''surface'') == 1;set(axis_children(i),''EdgeColor'',''none'',''FaceColor'',''flat'');end;end;clear i;clear axis_children;');
uimenu(Shading,'Label','Interpolated', ...
	'Callback','axis_children=get(gca,''Children'');for i=1:length(axis_children);if strcmp(get(axis_children(i),''Type''),''surface'') == 1;set(axis_children(i),''EdgeColor'',''none'',''FaceColor'',''interp'');end;end;clear i;clear axis_children;');
Hiding=uimenu(Display,'Label','Hiding');
uimenu(Hiding,'Label','On','Callback','hidden on');
uimenu(Hiding,'Label','Off','Callback','hidden off');
mesh(e,n,v);
clear n e v;
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
set(gca,'ZDir','reverse');
ylabel(ncomment);
clear ncomment;
xlabel(ecomment);
clear ecomment;
zlabel(vcomment);
clear ecomment;
title(titlecomment)
clear titlecomment;
caxis('manual');
caxis([min_color max_color]);
clear min_color max_color;
grid on;
hold off;
rotate3d on;
hidden off;