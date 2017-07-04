function [otime,tshift,hdur,lat,lon,dep,M,eid,elabel] = read_CMTSOLUTION(filename,nlines_cmt,nspace_between_entries)
%READ_CMTSOLUTION read a CMTSOLUTION format file (or a concatenated version)
%
% This file reads a concatenated CMTSOLUTION file and outputs the data.
% The CMT tshift parameter is added to the origin time, which is the
% time used for the initial CMT solution.
%
% INPUT
%   filename                 file containing all CMTSOLUTION files concatenated together
%   nlines_cmt               optional: number of lines per CMT solution (13)
%   nspace_between_entries   optional: number of spaces between CMTSOLUTION blocks (0 or 1)
%
% OUTPUT
%   otime   
%   tshift  
%   hdur    
%   lat     
%   lon     
%   dep     
%   M       6 x n moment tensors in CMT convention
%           M = [Mrr Mtt Mpp Mrt Mrp Mtp]; r=up, theta=south, phi=east
%           WARNING: units here are N-m, not dyne-cm as in CMTSOLUTION
%   eid     
%   elabel  
%
% EXAMPLE:
%   filename = '/home/carltape/PROJECTS/SPECFEM/CMTSOLUTION_finite_test';
%   [otime,tshift,hdur,lat,lon,dep,M,eid,elabel] = read_CMTSOLUTION(filename,13,0);
%
% See the wrapper function read_CMTSOLUTION_finite.m
%
% Carl Tape, 06/26/2007
%

% default optional arguments
if nargin==1
    nlines_cmt = 13;
    nspace_between_entries = 0;
end

% read in concatenated CMTSOLUTION files
lines = textread(filename,'%s','delimiter','\n');
nlines = length(lines);

% number of events
enum = (nlines + nspace_between_entries) / (nlines_cmt + nspace_between_entries);
if mod(enum,1) ~= 0
    disp('read_CMTSOLUTION.m: mismatch of expected number of lines');
    nlines, nspace_between_entries, nlines_cmt, nspace_between_entries
    error('enum = %.2f should be an integer',enum);
end
disp([' File : ' filename ]);
disp([' Number of events : ' num2str(enum) ]);

% initialize vectors
otime   = zeros(enum,1); tshift  = zeros(enum,1); hdur    = zeros(enum,1);
lat     = zeros(enum,1); lon     = zeros(enum,1); dep     = zeros(enum,1);
Mrr     = zeros(enum,1); Mtt     = zeros(enum,1); Mpp     = zeros(enum,1);
Mrt     = zeros(enum,1); Mrp     = zeros(enum,1); Mtp     = zeros(enum,1);

% initialize strings
eid = repmat(cellstr(''),enum,1); 
elabel = repmat(cellstr(''),enum,1); 

for kk = 1:enum
    if mod(kk,100)==0, disp(sprintf('%i/%i',kk,enum)); end
    in1 = (kk-1)*(nlines_cmt + nspace_between_entries);
    
    % first full line of CMTSOLUTION, as a string
    ltemp = lines{in1+1};
    
    if isempty(ltemp), error('the CMTSOLUTION line is length zero'); end
    
    % debug
    %disp(sprintf('kk = %i in1 = %i LINE: %s',kk,in1,ltemp));
    
    % replace W with space for more recent CMT solutions
    % SWEQ2006
    % PDEW2006
    if or(strcmp(ltemp(4),'W'),strcmp(ltemp(4),'Q')), ltemp(4) = ' '; end
    
    [j1,yr,mo,dy,hr,min,sec,lat_pde(kk),lon_pde(kk),dep_pde(kk),j5,j6,name1,name2,name3,name4,name5,name6] = ...
       strread(ltemp,'%s%f%f%f%f%f%f%f%f%f%f%f%s%s%s%s%s%s');
    % note: if there is no space after [event name:], then the following
    % line will crash, since it expects three fields to be present
    [j1,j2,eid(kk)] = strread(lines{in1+2},'%s%s%s');
    [j1,j2,tshift(kk)] = strread(lines{in1+3},'%s%s%f');
    [j1,j2,hdur(kk)] = strread(lines{in1+4},'%s%s%f');
    [j1,lat(kk)] = strread(lines{in1+5},'%s%f');
    [j1,lon(kk)] = strread(lines{in1+6},'%s%f');
    [j1,dep(kk)] = strread(lines{in1+7},'%s%f');
    [j1,Mrr(kk)] = strread(lines{in1+8},'%s%f');
    [j1,Mtt(kk)] = strread(lines{in1+9},'%s%f');
    [j1,Mpp(kk)] = strread(lines{in1+10},'%s%f');
    [j1,Mrt(kk)] = strread(lines{in1+11},'%s%f');
    [j1,Mrp(kk)] = strread(lines{in1+12},'%s%f');
    [j1,Mtp(kk)] = strread(lines{in1+13},'%s%f');

    % NOTE: time shift is added to the PDE origin time; thus the stored
    %       origin time is the centroid time obtained from the CMT inversion
    otime(kk) = datenum(yr,mo,dy,hr,min, sec + tshift(kk));
    
    elabel{kk} = [char(name1) ' ' char(name2) ' ' char(name3) ...
             ' ' char(name4) ' ' char(name5) ' ' char(name6)];
end

eid = eid(:);
elabel = elabel(:);

% moment tensor elements
M = [Mrr Mtt Mpp Mrt Mrp Mtp];

% convert moment tensor from dyne-cm to N-m
M = 1e-7 * M;

% convert to 6 x enum
M = M';

%==========================================================================
