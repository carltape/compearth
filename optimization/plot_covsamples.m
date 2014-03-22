function plot_covsamples(msamples,rho,tlab,msamples2,rho2,tlab2,mlabs)
% PLOT_COVMSAMPLES generates plots for samples of covariance matrix
%
% INPUT
%    msamples   n x s matrix of vector samples
%    rho        n x n 'analytical' correlation matrix
%    tlab       label for plot
%    msamples2  optional: 2nd set of samples ([] for none)
%    rho2       optional: 2nd 'analytical' correlation matrix ([] for none)
%    tlab2      optional: label for plot
%    mlabs      optional: labels for each variable ([] for default)
%
% EXAMPLE: 
%    plot_covsamples(mpost_samples,rho_post,'mpost',[],[],[],mlabs);
%
% NOTE: We could alternatively estimate the covariance matrix
% (and correlation matrix) directly from the input samples.
%
% called by optimization.m, optim_hw.m
%
% Carl Tape, 01-Jan-2012
%

NMAX = 6;   % max number to make into scatterplot

[n,s] = size(msamples);
disp(sprintf('plot_covsamples.m: n = %i, s = %i',n,s));

if isempty(mlabs)
    %mlabs = strtrim(cellstr(num2str([1:n]')));
    mlabs = repmat(cellstr(''),n,1);
    for ii=1:n, mlabs{ii} = sprintf('i%i',ii); end
end

% whether to plot a second set of samples
if and(~isempty(msamples2),~isempty(rho2)), isecond = 1; else isecond = 0; end

% kk=1: correlation matrices from Cpost
% kk=2: correlation matrices based on inpur SAMPLES
figure; nr=1+isecond; nc=2;
for kk=1:2
    if kk==1
       F1 = rho; 
       if isecond==1, F2 = rho2; end
       stag = '';
    else
       F1 = corrcoef(msamples');
       if isecond==2, F2 = corrcoef(msamples2'); end
       stag = 'sample';
    end

    % note: we could replace imagesc with a non-toolbox function (pcolor)
    pind = kk+isecond*(kk-1);
    subplot(nr,nc,pind); imagesc(F1); caxis([-1 1]), colorbar
    set(gca,'xtick',[1:n],'xticklabel',mlabs,'xaxislocation','top');
    set(gca,'ytick',[1:n],'yticklabel',mlabs);
    title([stag ' correlation matrix for ' tlab]); axis equal, axis tight

    if isecond==1
    subplot(nr,nc,pind+1); imagesc(F2); caxis([-1 1]), colorbar
    set(gca,'xtick',[1:n],'xticklabel',mlabs,'xaxislocation','top');
    set(gca,'ytick',[1:n],'yticklabel',mlabs);
    title([stag ' correlation matrix for ' tlab2]); axis equal, axis tight
    end
end

% scatterplots
if n > NMAX
    disp(sprintf('n = %i is > %i, so no scatterplots made',n,NMAX));
else
    figure; nr=n-1; nc=n-1;
    for ii=1:n-1
       for jj=ii+1:n
           px = msamples(ii,:);
           py = msamples(jj,:);
           iplot = nc*(ii-1) + jj-1;
           %disp([ii jj iplot]);
           subplot(nr,nc,iplot); hold on;
           plot(px,py,'.');
           xlabel(mlabs{ii}); ylabel(mlabs{jj});
           st1 = sprintf('corr(%s) = %.2f (%.2f)',tlab,corr(px(:),py(:)),rho(ii,jj));
           if isecond==1
               px = msamples2(ii,:);
               py = msamples2(jj,:);
               plot(px,py,'r.');
               st2 = sprintf('corr(%s) = %.2f (%.2f)',tlab2,corr(px(:),py(:)),rho2(ii,jj));
               title({st1,st2});
           else
               title(st1);
           end
       end
    end
end

%==========================================================================