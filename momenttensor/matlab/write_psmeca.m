function write_psmeca(filename,otime,lat,lon,dep,M,eid,slabel)
% WRITE_PSMECA write psmeca file for GMT plotting
%
% INPUT
%   filename    (files will be appended with _year, _Mw, _date, etc)
%   otime       n x 1 origin time (Matlab format)
%   lat         n x 1 latitude
%   lon         n x 1 longitude
%   dep         n x 1 depth, in km
%   M           6 x n array of moment tensors in GCMT convention (up-south-east)
%                   M = [Mrr Mtt Mpp Mrt Mrp Mtp], units N-m
%   eid         OPTIONAL: event IDs
%   slabel      OPTIONAL: names for events (e.g., CENTRAL_ALASKA)
%
% This inputs a set of moment tensors in units of N-m and outputs four
% different files for plotting in GMT using psmeca, which assumes units of
% dyne-cm. The output files differ only in what label to use plotting above
% the beach balls.
%
% slabel determines the text label for each moment tensor
% If slabel is present, then two files are written: (1) no labels (2) slabel
% If slabel is not present, then five files are written:
%   (1) no label (2) eid (3) Mw (4) year (5) date
%
% calls CMT2m0.m
%
% Carl Tape, 02-Feb-2011
%

n = length(otime);

% make sure M is 6 x n
[M,n1] = Mdim(M);

% check input argument dimensions
if n1~=n, whos M otime, error('dimension mismatch (M otime)'); end
if length(lat)~=n, whos lat otime, error('dimension mismatch (lat otime)'); end
if length(lon)~=n, whos lon otime, error('dimension mismatch (lon otime)'); end
if length(dep)~=n, whos dep otime, error('dimension mismatch (dep otime)'); end

% convert moment tensor from N-m to dyne-cm
M = 1e7 * M;

% exponent for computing magnitude in psmeca
M0 = CMT2m0(1,M);
iexp_all = floor(log10(M0));

% for labeling the moment magnitude
Mw = m02mw(1,CMT2m0(1,M*1e-7));   % M0 must be in N-m

% controls the labels for the beach balls
if nargin < 8
    % slabel not provided: write 5 different files
    imin = 1; imax = 5;
    % eid is not provided
    if nargin == 6
        eid = strtrim(cellstr(num2str([1:n]')));
    end
else
    % slabel IS provided: one file with no label (5), one file with slabel (6)
    imin = 5; imax = 6;
    if length(slabel)~=n, whos slabel otime, error('dimension mismatch (slabel otime)'); end
end

if length(eid)~=n, whos eid otime, error('dimension mismatch (eid otime)'); end

%--------------------

for ilab = imin:imax
    
    if ilab==0, ext = ''; end
    if ilab==1, ext = '_eid'; end
    if ilab==2, ext = '_Mw'; end
    if ilab==3, ext = '_year'; end
    if ilab==4, ext = '_date'; end
    if ilab==5, ext = ''; end 
    if ilab==6, ext = '_custom'; end 
    
    disp(['write_psmeca.m : extension is ' ext]);

    % write to file for GMT plotting
    file1 = [filename '_psmeca' ext];
    fid = fopen(file1,'w');
    for ii = 1:n

        % title for beach ball
        switch ilab
            case 1, cmtlabel = char(eid{ii});
            case 2, cmtlabel = sprintf('%.2f',Mw(ii));  
            case 3, cmtlabel = datestr(otime(ii),10); 
            case 4, cmtlabel = datestr(otime(ii),29); 
            case 5, cmtlabel = '  ';
            case 6, cmtlabel = slabel{ii};    
            %case 5, cmtlabel = sprintf('%4i-M%.1f-%i',year(otime(ii)),Mw(ii),isource(ii)); 
        end

        fac = 10^-iexp_all(ii);

        % originally had 16 char for final string; changed to open
        if ilab==5
            fprintf(fid,'%14.6f%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%4i\n',...
                lon(ii), lat(ii), dep(ii),...
                M(1,ii)*fac, M(2,ii)*fac, M(3,ii)*fac,...
                M(4,ii)*fac, M(5,ii)*fac, M(6,ii)*fac,...
                iexp_all(ii));
        else
            fprintf(fid,'%14.6f%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%4i%14.6e%14.6e %s\n',...
                lon(ii), lat(ii), dep(ii),...
                M(1,ii)*fac, M(2,ii)*fac, M(3,ii)*fac,...
                M(4,ii)*fac, M(5,ii)*fac, M(6,ii)*fac,...
                iexp_all(ii),...
                lon(ii), lat(ii), cmtlabel);
        end
    end
    fclose(fid);
end

% write a list of event IDs
file2 = [filename '_eid'];
fid = fopen(file2,'w');
for ii=1:n, fprintf(fid,'%s\n',char(eid{ii})); end
fclose(fid);

disp(' writing psmeca file for GMT plotting...');
disp([' output file : ' file1]);
disp([' output file : ' file2]);
disp([' number of CMT solutions : ' num2str(n)]);

%======================================================
