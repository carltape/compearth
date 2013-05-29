function displayCMTshort(M,stfmt)
%DISPLAYCMTSHORT display moment tensor elements

% make sure M is 6 x n
[M,n] = Mdim(M);

% default format
if nargin==1, stfmt = '%16.6e'; end

for kk=1:n
    disp('----------------');
    disp(sprintf(['M11:' stfmt],M(1,kk)));
    disp(sprintf(['M22:' stfmt],M(2,kk)));
    disp(sprintf(['M33:' stfmt],M(3,kk)));
    disp(sprintf(['M12:' stfmt],M(4,kk)));
    disp(sprintf(['M13:' stfmt],M(5,kk)));
    disp(sprintf(['M23:' stfmt],M(6,kk)));
end
disp('====================');