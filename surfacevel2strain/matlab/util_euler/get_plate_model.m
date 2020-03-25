function [exyz,names,name_labs,dir_bounds,ssfx,smod] = get_plate_model(imodel,bdisplay)
%GET_PLATE_MODEL return a plate model in the form of euler poles and plate names
%
% EXAMPLE: get_plate_model(1,true);
%
% called by platemodel2gps.m, platemodel2conv_vel.m
%
% Carl Tape, 2006-05-03
%

if nargin==1, bdisplay=false; end

% string labels for each plate model
mod_labs    = {'oneill','nuvel1A_nnr','revel','bird','gripp_hs3','bird_gripp','bird_morgan','bird_nnr','morvel_morgan','morvel_nnr'};

% directories for plate boundaries for each model
% NOTE: some models use the same plate boundaries
dir_labs    = {'oneill/data_mike_13/','nuvel1A/','bird/','bird/','nuvel1A/','bird/','bird/','bird/','morvel_nnr/','morvel_nnr/'};

% suffixes for plate boundary files for each model
suffix_labs = {'.xy','.180to180xy','.xy','.xy','.180to180xy','.xy','.xy','.xy','.xy','.xy'};

% imodel is input variable
if imodel > 10
    smod = sprintf('morvel_nnr_becker2015_%2.2i',imodel-10);
    sbnd = 'morvel_nnr/';
    ssfx = '.xy';
else
    smod = mod_labs{imodel};
    sbnd = dir_labs{imodel};
    ssfx = suffix_labs{imodel};
end

dir_plates  = '/home/carltape/gmt/plates/';
dir_models  = [dir_plates 'plate_models/' smod '/'];
dir_bounds  = [dir_plates 'plate_boundaries/' sbnd];

% load euler vectors and plate names (these are made in run_platemodel2gps.m)
ww = [dir_models smod '_euler_poles.dat'];
if ~exist(ww,'file'), error('file not found: %s',ww); end
[wx,wy,wz,name_labs,names] = textread(ww,'%f%f%f%s%s');
exyz = [wx wy wz]';     % euler poles
nump = length(wx);      % number of plates

% display info on euler poles
if bdisplay
    % convert euler poles (wx,wy,wz) --> (lat,lon,omg)
    outvec1 = euler_convert(exyz,1);
    elat = outvec1(1,:);
    elon = outvec1(2,:);

    % rotation rate, deg/Myr
    omegs = sqrt( exyz(1,:).^2 + exyz(2,:).^2 + exyz(3,:).^2 );

    % maximum surface velocity on a plate, in mm/yr
    % vmax will be at the points that are del=90 from the euler pole
    % vmax = omegs/mmyr2degMyr;
    
    % display info on euler poles
    disp('  ');
    disp(sprintf('Here are %i euler poles from the %s model:',nump,smod));
    disp('             label     elon        elat    deg/Myr    name'); 
    for ii=1:nump  
        disp(sprintf('  Plate #%2.2i : %s%12.4f%12.4f%10.4f   %s',...
            ii,name_labs{ii},elon(ii),elat(ii),omegs(ii),names{ii}));
    end
    disp('-----------------------------');
end

%==========================================================================