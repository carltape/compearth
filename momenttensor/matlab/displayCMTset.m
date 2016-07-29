function [M0,Mw,M0tot,Mwtot] = displayCMTset(M,tshift,hdur,zele)
%DISPLAYCMTSET display information about a set of moment tensors
%
% It was written with finite-source CMTSOLUTION files in mind.
%
% calls CMT2m0.m, m02mw.m
% called by read_CMTSOLUTION_finite.m
%

% number of subsources
n = length(tshift);

% compute M0 and Mw for each subsource
M0 = CMT2m0(1,M);
Mw = m02mw(1,M0);
M0tot = sum(M0);
Mwtot = m02mw(1,M0tot);
    
disp('----------------------------------');
disp('summary from displayCMTset.m');
%disp(sprintf('file: %s',fname));
disp(sprintf('number of subsources: %i',n));
disp(sprintf('min/median/max Mw of subsources: %.2f / %.2f / %.2f',...
    min(Mw),median(Mw),max(Mw)));
disp(sprintf('total M0: %.3e N-m',M0tot));
disp(sprintf('total Mw: %.3f',Mwtot));
disp(sprintf('min/median/max elevation of subsources:'));
disp(sprintf('    %.3f km / %.3f km / %.3f km',...
    min(zele)/1000,median(zele)/1000,max(zele)/1000));
disp(sprintf('min/median/max tshift of subsources:'));
disp(sprintf('         %.1f s / %.1f s / %.1f s',...
    min(tshift),median(tshift),max(tshift)));
disp(sprintf('      %.2f min / %.2f min / %.2f min',...
    min(tshift)/60,median(tshift)/60,max(tshift)/60));
disp(sprintf('duration of rupture:'));
disp(sprintf('        %5.2f s',max(tshift)-min(tshift)));
disp(sprintf('        %5.2f min',(max(tshift)-min(tshift))/60 ));
disp(sprintf('min/median/max hdur of subsources:'));
disp(sprintf('         %.1f s / %.1f s / %.1f s',min(hdur),median(hdur),max(hdur)));
disp('----------------------------------');

%==========================================================================
