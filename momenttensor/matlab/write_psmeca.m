function write_psmeca(filename,otime,lat,lon,dep,M,eid,slabel,slabeltag)
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
%   slabel      OPTIONAL: custom labels for events
%   slabeltag   OPTIONAL: file tag for custom text
%
% This inputs a set of moment tensors in units of N-m and outputs four
% different files for plotting in GMT using psmeca, which assumes units of
% dyne-cm. The output files differ only in what label to use plotting above
% the beachballs.
%
% slabel determines the text label for each moment tensor
% If slabel is present, then two files are written: (1) no labels (2) slabel
% If slabel is not present, then five files are written:
%   (1) no label (2) eid (3) Mw (4) year (5) date
%
% calls CMT2m0.m
%
% Carl Tape, 2011-02-02
%

n = length(lat);

% make sure M is 6 x n
[M,n1] = Mdim(M);

% check input argument dimensions
if n1~=n, whos M lat, error('dimension mismatch (M lat)'); end
if length(lon)~=n, whos lon lat, error('dimension mismatch (lon lat)'); end
if length(dep)~=n, whos dep lat, error('dimension mismatch (dep lat)'); end

% convert moment tensor from N-m to dyne-cm
M = 1e7 * M;

% exponent for computing magnitude in psmeca
M0 = CMT2m0(1,M);       % dyne-cm
iexp_all = floor(log10(M0));

% for labeling the moment magnitude
Mw = m02mw(1,CMT2m0(1,M*1e-7));   % M0 must be in N-m

% controls the labels for the beachballs
if nargin < 9
    % slabel not provided: write 5 different files
    imin = 1; imax = 5;
    % eid is not provided
    if nargin == 6
        eid = strtrim(cellstr(num2str([1:n]')));
    end
else
    % slabel IS provided: one file with no label (5), one file with slabel (6)
    imin = 5; imax = 6;
    if length(slabel)~=n
        n, whos slabel lat, error('dimension mismatch (slabel lat)');
    end
    if isempty(slabeltag);
        scustom = 'custom';
    else
        scustom = slabeltag;
    end
end

% if no origin times are specified, then just plot one output file
%if isempty(otime), imin=5; imax=5; end

if isempty(eid), eid = strtrim(cellstr(num2str([1:n]'))); end

%if length(eid)~=n
%    n, whos eid otime, error('dimension mismatch (eid otime)');
%end

%--------------------

nfile = imax-imin+1;
disp(sprintf('looping from %i to %i for %i files',imin,imax,nfile));

for ilab = imin:imax            % loop over different text labels
    
    if ilab==0, ext = ''; end
    if ilab==1, ext = '_eid'; end
    if ilab==2, ext = '_Mw'; end
    if and(ilab==3,~isempty(otime)) ext = '_year'; end
    if and(ilab==4,~isempty(otime)) ext = '_date'; end
    if ilab==5, ext = ''; end 
    if ilab==6, ext = strcat('_',scustom); end 
    
    disp(['write_psmeca.m : extension is ' ext]);

    % write to file for GMT plotting
    file1 = [filename '_psmeca' ext];
    fid1 = fopen(file1,'w');
    for ii = 1:n                % loop over events

        % title for beachball
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

        if ilab==5
            stfmt = '%14.6f%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%4i\n';
            fprintf(fid1,stfmt,...
                lon(ii), lat(ii), dep(ii),...
                M(1,ii)*fac, M(2,ii)*fac, M(3,ii)*fac,...
                M(4,ii)*fac, M(5,ii)*fac, M(6,ii)*fac,...
                iexp_all(ii));
        else
            % originally had 16 char for final string; changed to open
            stfmt = '%14.6f%14.6f%14.6f%14.6e%14.6e%14.6e%14.6e%14.6e%14.6e%4i%14.6e%14.6e %s\n';
            %stfmt = '%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%14.6f%4i%14.6f%14.6f %s\n';
            fprintf(fid1,stfmt,...
                lon(ii), lat(ii), dep(ii),...
                M(1,ii)*fac, M(2,ii)*fac, M(3,ii)*fac,...
                M(4,ii)*fac, M(5,ii)*fac, M(6,ii)*fac,...
                iexp_all(ii),...
                lon(ii), lat(ii), cmtlabel);
        end
    end
    fclose(fid1);
end

disp('writing psmeca file for GMT plotting...');
disp(sprintf('number of CMT solutions : %i',n));
disp(sprintf('output file (1 of %i): %s',nfile,file1));

if ~isempty(otime)
    % write a list of event IDs (useful to have for other scripts)
    file_eid = [filename '_eid'];
    fid = fopen(file_eid,'w');
    for ii=1:n, fprintf(fid,'%s\n',char(eid{ii})); end
    fclose(fid);
    disp(sprintf('output file: %s',file_eid));
end

%==========================================================================
