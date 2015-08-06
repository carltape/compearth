function plotMT_Mij(M,edges,Mlab)
%PLOTMT_MIJ plot histograms of the entries of Mij
%
% M      6 x n moment tensors, M = [M11 M22 M33 M12 M13 M23]
% edges  Ne-dimensional vector for bin edges (or 6 x Ne sets of edges)
% Mlab   optional: labels for each subplot
%

% make sure M is 6 x n
[M,n] = Mdim(M);

if nargin~=3
    Mlab = {'M11','M22','M33','M12','M13','M23'};
end

% edges can be specified for all six entries
% (for now we require the same number of bins)
USE_DIFFERENT_EDGES = false;
if ~any(size(edges)==1)
    USE_DIFFERENT_EDGES = true;
    [n1,n2] = size(edges);
    if n2==6, edges=edges'; end     % default 6 x Ne
end

figure; nr=3; nc=2;
pvec = [1 3 5 2 4 6];
for ii=1:6
   subplot(nr,nc,pvec(ii));
   if USE_DIFFERENT_EDGES
       plot_histo(M(ii,:),edges(ii,:));
   else
       plot_histo(M(ii,:),edges);
   end
   xlabel(Mlab{ii});
end

Mdiag = M(1:3,:);
Moff  = M(4:6,:);

disp('plotMT_Mij.m:');
disp(sprintf('Diagonal entries    : min(Mij) = %f, max(Mij) = %f',min(Mdiag(:)),max(Mdiag(:))));
disp(sprintf('Off-diagonal entries: min(Mij) = %f, max(Mij) = %f',min(Moff(:)),max(Moff(:))));

if min(M(:)) < edges(1),
    disp(sprintf('plotMT_Mij.m WARNING: min(M) < edges(1) (%f)',edges(1)));
end
if max(M(:)) > edges(end),
    disp(sprintf('plotMT_Mij.m WARNING: max(M) < edges(2) (%f)',edges(end)));
end

%==========================================================================
