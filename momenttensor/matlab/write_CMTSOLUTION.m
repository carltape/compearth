function write_CMTSOLUTION(dir0,ione_file,otime,tshift,hdur,lat,lon,dep,M,  eid,elabel,ftag,nspace_between_entries)
%WRITE_CMTSOLUTION outputs files in CMTSOLUTION format
%
% NOTE: In read_CMTSOLUTION.m we apply the time shift from CMT into the (absolute)
% origin time; thus for most purposes, tshift should be set to zero.
%
% In a typical CMTSOLUTION file, the only important values entered on the
% first line describe the origin time. The other values are based on the
% initial source estimate that Harvard uses to seed their solution.
%
% INPUT - required:
%   dir0        output directory
%   ione_file   =1 to write to a single file (finite source model)
%               =0 to write to individual files
%   otime       origin time
%   tshift      time shift (crucial for finite source models)
%   hdur        half duration (typically scaled from M0)
%   lat,lon     source epicenter
%   dep         source depth in kilometers (positive values point down)
%   M           6 x n set of moment tensors, M = [Mrr Mtt Mpp Mrt Mrp Mtp]
%                  units in N-m
% INPUT - optional:
%   eid         event id
%   elabel      event description
%   ftag        label for single output file (CMTSOLUTION_ftag)
%   nspace_between_entries
%               number of spaces between CMTSOLUTION blocks (0 or 1)
%               for single output file only
%
% Only 9 or 13 arguments are permissible. If 13 arguments are entered, then
% any of the optional arguments can be set to [].
%
% See also read_CMTSOLUTION.m
%
% Carl Tape, 10/2008
%
% calls CMT2m0.m, m02mw.m
% called by test_CMT.m
%

disp('entering write_CMTSOLUTION.m');
enum = length(otime);

% check for empty input arguments
if isempty(tshift)
    tshift = zeros(enum,1);
end
if isempty(hdur)
    hdur = zeros(enum,1);   % we could compute hdur from M0 rather than list zero as default
end

% check dimensions of input variables
if length(unique([ length(otime) length(tshift) length(hdur) length(lat) length(lon) length(dep) ])) ~= 1
    whos otime tshift hdur lat lon dep
    error('input vectors have different lengths');
end
[M,n] = Mdim(M);
if n~=enum, error('M has %i entries, enum = %i',n,enum); end

% convert moment tensor from N-m to dyne-cm (1 N-m = 10^7 dyne-cm)
M = 1e7 * M;

% check OPTIONAL input arguments -- set to defaults if they are missing
if or(nargin==9,nargin==13)
    if nargin==9, eid=[]; elabel=[]; ftag=[]; nspace_between_entries=[]; end
    if isempty(eid)
        %for ii=1:enum, eid{ii} = sprintf('%6.6i',ii); end
        eid = cellstr(num2str([1:enum]','%6.6i'));
        disp('no eid specified: setting to default');
    end
    if isempty(elabel)
        elabel = eid;
        disp('no elabel specified: setting elabel = eid');
    end
    if isempty(ftag)
        ftag = '';
    end
    if isempty(nspace_between_entries)
        nspace_between_entries = 0;
    end
else
    error('write_CMTSOLUTION.m: only 9 or 13 arguments is permissible (some of which can be empty [])');
end

if ione_file == 1       % write to a single file
    
    ofile = [dir0 'CMTSOLUTION_' ftag];
    disp(sprintf('write_CMTSOLUTION.m: writing %s',ofile));
    fid = fopen(ofile,'w');
    for kk = 1:enum
        if mod(kk,100)==0, disp(sprintf('%i/%i',kk,enum)); end
        m0 = CMT2m0(1,M(:,kk));
        mag = m02mw(1,1e-7*m0);  % assumes M is in dyne-cm
       
        % SEE NOTES BELOW
        fprintf(fid,'XXXX %4i %2i %2i %2i %2s %5.2f %8.4f %9.4f %5.1f %3.1f %3.1f %s\n',...
            year(otime(kk)),month(otime(kk)),day(otime(kk)),...
            hour(otime(kk)),datestr(otime(kk),'MM'),second(otime(kk)),...
            lat(kk),lon(kk),dep(kk),mag,mag,elabel{kk});
        
        fprintf(fid,'event name: %11s\n',eid{kk});
        fprintf(fid,'time shift: %11.4f\n',tshift(kk));
        fprintf(fid,'half duration: %8.4f\n',hdur(kk));
        % for SPECFEM3D you can use UTM coordinates, which are in meters
        if or(any(abs(lat) > 90),any(abs(lon) > 360))
            fprintf(fid,'latitude: %13.1f\n',lat(kk));
            fprintf(fid,'longitude: %12.1f\n',lon(kk));
        else
            fprintf(fid,'latitude: %13.4f\n',lat(kk));
            fprintf(fid,'longitude: %12.4f\n',lon(kk));
        end
        fprintf(fid,'depth: %16.4f\n',dep(kk));
        fprintf(fid,'Mrr: %18.6e\n',M(1,kk));
        fprintf(fid,'Mtt: %18.6e\n',M(2,kk));
        fprintf(fid,'Mpp: %18.6e\n',M(3,kk));
        fprintf(fid,'Mrt: %18.6e\n',M(4,kk));
        fprintf(fid,'Mrp: %18.6e\n',M(5,kk));
        fprintf(fid,'Mtp: %18.6e\n',M(6,kk));      % SPACE HERE OR NOT?
        if and(kk ~= enum, nspace_between_entries == 1)
            fprintf(fid,'\n');
        end
    end
    fclose(fid);
    
    % see read_CMTSOLUTION_finite.m
    % some finite source models are huge, so here we also save a matlab
    % file for quick loading
    %save(ofile,'otime','tshift','hdur','lat','lon','dep','M','eid','elabel');
    
else            % write to individual files
    disp(sprintf('write_CMTSOLUTION.m: writing %i CMTSOLUTION files',enum));
    disp(sprintf('  output directory: %s',dir0));
    for kk = 1:enum
        fid = fopen([dir0 'CMTSOLUTION_' eid{kk}],'w');
        m0 = CMT2m0(1,M(:,kk));
        mag = m02mw(1,1e-7*m0);

        %fprintf(fid,'XXXX %4i %2i %2i %2i %2i %5.2f -LA.TUDE -LON.TUDE DPT.H M.B M.S EVENT_NAME\n',...
        %    year(otime(kk)),month(otime(kk)),day(otime(kk)),hour(otime(kk)),minute(otime(kk)),second(otime(kk)));

       % NOTE: (1) cannot use the minute command, since it will round up (see below).
       %       (2) but we need the second command to avoid round-off errors
       %           of milliseconds from using datevec
        fprintf(fid,'XXXX %4i %2i %2i %2i %2s %5.2f %8.4f %9.4f %5.1f %3.1f %3.1f %s\n',...
            year(otime(kk)),month(otime(kk)),day(otime(kk)),...
            hour(otime(kk)),datestr(otime(kk),'MM'),second(otime(kk)),...
            lat(kk),lon(kk),dep(kk),mag,mag,elabel{kk});
        
        fprintf(fid,'event name: %s\n',eid{kk});
        fprintf(fid,'time shift: %11.4f\n',tshift(kk));
        fprintf(fid,'half duration: %8.4f\n',hdur(kk));
        % for SPECFEM3D you can use UTM coordinates, which are in meters
        if or(any(abs(lat) > 90),any(abs(lon) > 360))
            fprintf(fid,'latitude: %13.1f\n',lat(kk));
            fprintf(fid,'longitude: %12.1f\n',lon(kk));
        else
            fprintf(fid,'latitude: %13.4f\n',lat(kk));
            fprintf(fid,'longitude: %12.4f\n',lon(kk));
        end
        fprintf(fid,'depth: %16.4f\n',dep(kk));
        fprintf(fid,'Mrr: %18.6e\n',M(1,kk));
        fprintf(fid,'Mtt: %18.6e\n',M(2,kk));
        fprintf(fid,'Mpp: %18.6e\n',M(3,kk));
        fprintf(fid,'Mrt: %18.6e\n',M(4,kk));
        fprintf(fid,'Mrp: %18.6e\n',M(5,kk));
        fprintf(fid,'Mtp: %18.6e\n',M(6,kk));   % LINE BREAK HERE OR NOT?
        fclose(fid);
    end
    
    % write list of event IDs
    if enum > 1
        disp(sprintf('write_CMTSOLUTION.m: writing %i event IDs to file',enum));
        disp(sprintf('  output directory: %s',dir0));
        fid = fopen(sprintf('%s/eids_%s',dir0,ftag),'w');
        for kk = 1:enum
            fprintf(fid,'%s\n',eid{kk});
        end
        fclose(fid);
    end
end

%-----------------------------------

if 0==1
    % basic example
    M = 3e16*[1 1 1 0 0 0]';
    otime = datenum(3000,1,1);
	write_CMTSOLUTION('./',0,otime,0,0,34.1081,-118.9683,90,M);
    
    % TECHNICAL NOTE
    % The command 'minute' (and perhaps the others)
    % will round up instead of using floor.
    % WHAT ABOUT THE OTHER COMMANDS, like year, month, day, etc?
    n0 = 7.310495416587963e+05;
    hsec = 0.5/3600/24;                 % half second
    nvec = [n0 ; n0+hsec];
    
    for kk=1:2
        n = nvec(kk);
        disp('-------------------');
        fprintf('%s %s %i %s\n',datestr(n,31),datestr(n,'FFF'),minute(n),datestr(n,'MM'));
        [y,m,d,h,mi,s] = datevec(n);
        fprintf('%i %i %i %i %i %.4f\n',y,m,d,h,mi,s);
        fprintf('%i %i %i %i %i %.4f\n',year(n),month(n),day(n),hour(n),minute(n),second(n));
    end
    
end

%==========================================================================
