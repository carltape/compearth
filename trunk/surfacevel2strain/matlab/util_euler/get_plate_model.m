%
% get_plate_model.m
% Carl Tape, 03-May-2006
%
% This script is a header file that assigns directory names for reading and
% writing data.
%
% calls xxx
% called by platemodel2gps.m, platemodel2conv_vel.m
%

% string labels for each plate model
mod_labs    = {'oneill','nuvel1A_nnr','revel','bird','gripp_hs3','bird_gripp','bird_morgan','bird_nnr'};

% directories for plate boundaries for each model
% NOTE: some models use the same plate boundaries
dir_labs    = {'oneill/data_mike_13/','nuvel1A/','bird/','bird/','nuvel1A/','bird/','bird/','bird/'};

% suffixes for plate boundary files for each model
suffix_labs = {'.xy','.180to180xy','.xy','.xy','.180to180xy','.xy','.xy','.xy'};

% imodel is input variable
smod = mod_labs{imodel};
sbnd = dir_labs{imodel};
ssfx = suffix_labs{imodel};

dir_plates  = '/home/carltape/gmt/plates/';
dir_models  = [dir_plates 'plate_models/' smod '/'];
dir_bounds  = [dir_plates 'plate_boundaries/' sbnd];

% load euler vectors and plate names
ww = [dir_models smod '_euler_poles.dat'];
[wx,wy,wz,name_labs,names] = textread(ww,'%f%f%f%s%s');
exyz = [wx wy wz]';     % euler poles
nump = length(wx);      % number of plates

%=================================================================
