function lam = Tk2lam(T,k)
%TK2LAM converts Hudson1989 T-k values to eigenvalues sorted lam1 >= lam2 >= lam3
% 
% Carl Tape, 21-Oct-2011
%

T = T(:)';
k = k(:)';
n = length(k);

lam = zeros(3,n);

% Pearce and Rogers (1989), Eq. 1
lam(1,:) = (1 - abs(k)) .* min([2*ones(1,n) ; 2-T]) + 2*k;
lam(2,:) = (1 - abs(k)) .* max([-2*ones(1,n) ; -(2+T)]) + 2*k;
lam(3,:) = (1 - abs(k)) .* T + 2*k;

% eigenvalue convention lam1 >= lam2 >= lam3
lam = sort(lam,'descend');

%==========================================================================
% EXAMPLE

if 0==1
    numT = 13;
    numk = numT;
    Tvec = linspace(-1,1,numT);
    kvec = linspace(-1,1,numk);
    [T,k] = meshgrid(Tvec,kvec);
    T = T(:);
    k = k(:);
    lam = Tk2lam(T,k);
    
    % check
    [Tcheck,kcheck] = lam2Tk(lam);
    % we cannot recover epsilon (T) for isotropic points
    for ii=1:length(k)
       disp(sprintf('%3i k = %6.2f (%6.2f), T = %6.2f (%6.2f ), %6.2f%6.2f%6.2f',...
           ii,k(ii),kcheck(ii),T(ii),Tcheck(ii),lam(:,ii)'))
    end
    norm(Tcheck(:)-T(:))
    norm(kcheck(:)-k(:))
    
    [gamma,delta,M0,mu,lamdev,lamiso] = lam2lune(lam);
    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot(T,k,'.','markersize',12); grid on;
    xlabel('T'); ylabel('k');
    axis equal, axis([-1 1 -1 1]); axis tight
    title('Hudson (1989) T-k points in T-k space');
    
    subplot(nr,nc,2); plot(gamma,delta,'.','markersize',12); grid on;
    axis([-31 31 -91 91]);
    xlabel('gamma'); ylabel('90 - delta');
    title('Hudson (1989) T-k points on the fundamental lune');
end

%==========================================================================
