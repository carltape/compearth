%
% test_euler_rot_tec.m
% Carl Tape, 23-April-2007
%
% This program shows how euler rotations on a sphere work.  Each rotation
% is described in terms of a rotation matrix that can be thought of as a
% pseudo-vector.
%
% The conventions are suited to plate tectonics.
% The programs are based primarily on Cox and Hart (1986), p. 227.
% 
% calls
%   euler_rot_tec.m
%   euler2rotmat.m, rotmat2euler.m
%   euler_mapping.m
%   latlon2xyz.m, xyz2latlon.m
%   globefun3.m, plotfig.m, octpts.m
% called by xxx
%

clc
clear
close all
format short
format compact

% add path to additional matlab scripts (specify bdir)
user_path;

deg = 180/pi;

idisplay = 0;       % boolean: display more information

%-----------------------
% This checks that euler2rotmat.m and rotmat2euler.m
% are indeed inverse operations. 

% GENERATE FINITE ROTATION PSEUDO-VECTORS
% We intentionally take values out of the "standard" range, in order to
% show that this will not crash the program.  For example, if we specify a
% negative rotation angle, then the code will take the antipode vector with
% a positive rotation angle (see euler_mapping.m).
% (NOTE: If we wanted UNIFORMLY distrubuted poles on the sphere, then we
% would use a Fisher distribution.)
n = 1000;
omeg0 = randomvec(-180, 180, n);
lats0 = randomvec(-180, 180, n);
lons0 = randomvec(-360, 360, n);
%omeg0 = zeros(n,1);                % check that it works when omega = 0

for ii = 1:n
    evec0 = [lats0(ii) lons0(ii) omeg0(ii)]';   % input rotation
    [R,evec] =  euler2rotmat(evec0);            % NOTE: evec might be modified
    lats(ii) = evec(1);
    lons(ii) = evec(2);
    omeg(ii) = evec(3);
    evec_recover = rotmat2euler(R);             % recover the pole 
    
    % check the error
    errlat(ii) = norm(evec(1) - evec_recover(1));
    errlon(ii) = norm(evec(2) - evec_recover(2));
    erromeg(ii) = norm(evec(3) - evec_recover(3));
    err(ii) = norm(evec - evec_recover);
    
    if idisplay==1
        disp('----------');
        display([' random euler vector ' num2str(ii) ' out of ' num2str(n)]);
        disp([ evec0 evec evec_recover ]);
    end
end

figure; nr=3; nc=1;
subplot(nr,nc,1); hold on;
plot(lons0,lats0,'.'); plot(lons,lats,'ro');
%axis([-180 180 -90 90]); grid on;
title([num2str(n) ' poles of finite rotations (test_euler_rot_tec.m)'],'interpreter','none');
xlabel(' Longitude'); ylabel(' Latitude');

subplot(nr,nc,2);
plot(err,'.');
title(' error in recovered poles');

subplot(nr,3,7); plot(errlat,'.'); title(' error in latitude');
subplot(nr,3,8); plot(errlon,'.'); title(' error in longitude');
subplot(nr,3,9); plot(erromeg,'.'); title(' error in omega');

orient tall, wysiwyg, fontsize(10);

%-----------------------
% example from Cox and Hart (1986), p. 227

for jj=1:2

    % test points
    if jj==1
        lat = 20; lon = 130;
        initial = latlon2xyz(lat,lon,1);
    else
        %octpts; initial = [A1 A2 A3];
        initial = octpts;
        [lat, lon] = xyz2latlon(initial);
    end

    % euler pole
    elat = -37; elon = 312; omega = 65;
    %evec = [elat elon omega]';
    evec = euler_mapping([elat elon omega]');
    exyz = latlon2xyz(elat,elon,1);

    % apply rotation
    [lat_rot, lon_rot, R] = euler_rot_tec(lat, lon, evec);
    initial_rot = latlon2xyz(lat_rot,lon_rot,1);

    disp('------------------------------------');
    disp('Here is an euler rotation:');
    disp(' input (lat, lon)  : '); disp([lat lon]);
    disp(' rotation matrix   : '); disp(R);
    disp(' norm(R) : '); disp(norm(R));
    disp(' output (lat, lon) : '); disp([lat_rot lon_rot]);

    % check the reverse -- un-rotate the rotated points
    evec = [elat elon -omega];
    [lat_orig, lon_orig, R] = euler_rot_tec(lat_rot, lon_rot, evec);

    disp('  '); disp('Now check the reverse operation:');
    disp(' input (lat, lon)  : '); disp([lat_rot lon_rot]);
    disp(' rotation matrix   : '); disp(R);
    disp(' norm(R) : '); disp(norm(R));
    disp(' output (lat, lon) : '); disp([lat_orig lon_orig]);
    disp('------------------------------------');

    %-----------------------

    ste = [' euler (lat, lon, \Omega)  =  ('  sprintf('%.2f', elat) ', ' ...
                 sprintf('%.2f', elon) ', ' sprintf('%.2f', omega) ')'];

    figure; nr=2; nc=1;
    for ii=1:2
        subplot(nr,nc,ii); hold on; 
        opts = [2 4 1 1 0];
        switch ii
            case 1, plotfig(initial, opts); title('initial points');   
            case 2, plotfig(initial_rot, opts); title({'final points', ['RED DOT: ' ste]}); 
        end
        globefun3(1,elat,elon,1,'r');
    end
    orient tall, wysiwyg, fontsize(10);
end

%==============================================================
