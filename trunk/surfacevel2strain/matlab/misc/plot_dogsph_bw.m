% 
% plot_dogsph_bw.m
% Pablo Muse
%
% This outputs values of a table for our 2009 GPS paper.
% 

close all, clear
clc, format short, format compact

%---------------------

n = 10;          % default 6
nterms = 2^11;   % default 2^9

figure;
A = cell(1,n);
[phi,theta] = sphgrid(2*nterms);

j = 1;
ff = dogsph(1/2^j,theta,phi);
fmat = fst(ff);
plot(abs(fmat(1,:)))
bw = zeros(3,n);
totalsumall = sum(sum(abs(ff).^2))
totalsum = sum(abs(fmat(1,:)).^2)
k = nterms;
cumsum = totalsum;
% while cumsum > 0.999*totalsum
%     cumsum = cumsum - abs(fmat(1,k))^2;
%     k = k-1;
% end
% bw(1,1) = 0.5; bw(2,1) = k;
while cumsum > 0.99*totalsum
    cumsum = cumsum - abs(fmat(1,k))^2;
    k = k-1;
end
bw(3,1) = k;

A{1} = 'a = 1/2';
hold all
for j = 2:n
    ff = dogsph(1/2^j,theta,phi); fmat = fst(ff); plot(abs(fmat(1,:)))
    totalsum = sum(abs(fmat(1,:)).^2)
    A{j} = ['a = 1/2^',num2str(j)];
    k = nterms;
    cumsum = totalsum;
%     while cumsum > 0.999*totalsum
%         cumsum = cumsum - abs(fmat(1,k))^2;
%         k = k-1;
%     end
%     bw(1,j) = 0.5; bw(2,j) = k;
    while cumsum > 0.99*totalsum
        cumsum = cumsum - abs(fmat(1,k))^2;
        k = k-1;
    end
    bw(3,j) = k;
end
legend(A);
bw

%==============================================================
