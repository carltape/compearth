function plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost)
% PLOT_EPICENTERS plots epicenters within the source-station geometry
%
% called by forward_epicenter.m, optimization.m
%
% Carl Tape, 24-Feb-2010
%

% input options
xrec = opts{1};
yrec = opts{2};
iray = opts{3};
ax0  = opts{4};
[nparm,nsample] = size(mprior_samples);
ndata = length(xrec);

% indices of xs and ys within model vector
ix = 1;
iy = 2;

figure; hold on;
msizer = 16;    % receiver size
msizes = 10;    % source size
rthick = 1;     % receiver edge thickness
sthick = 2;     % source edge thickness
rfsize = 10;

%xmin = ax0(1);
%xmax = ax0(2);
%ymin = ax0(3);
%ymax = ax0(4);

% plot ray paths
if iray==1
    for ii=1:ndata
        plot([minitial(ix) xrec(ii)],[minitial(iy) yrec(ii)],'k','linewidth',1);
    end
end

% plot option depends on if a posterior model is passed
if exist('mpost','var')
    p0 = plot(mprior_samples(ix,:),mprior_samples(iy,:),'.');
    p1 = plot(minitial(ix),minitial(iy),'o','markersize',msizes,'markerfacecolor','k','markeredgecolor','w','linewidth',sthick);
    p2 = plot(mpost(ix),mpost(iy),'o','markersize',msizes,'markerfacecolor','c','markeredgecolor','w','linewidth',sthick);
    pP = plot(mprior(ix),mprior(iy),'o','markersize',msizes,'markerfacecolor','b','markeredgecolor','w','linewidth',sthick);
    pT = plot(mtarget(ix),mtarget(iy),'o','markersize',msizes,'markerfacecolor','r','markeredgecolor','w','linewidth',sthick);
    plot(xrec,yrec,'kV','markersize',msizer,'linewidth',rthick);
    for ii=1:ndata, text(xrec(ii),yrec(ii),num2str(ii),'fontsize',rfsize,'color','r',...
            'horizontalalignment','center','verticalalignment','middle'); end
    legend([p0(1) p1 p2 pP pT],'Cprior sample','minitial','mpost','mprior','mtarget');
    %legend([p0(1) p1 p2 pP pT],'Cprior sample','m00',sprintf('m%2.2i',niter),'mprior','mtarget');
    
else
    if ~isempty(mprior_samples)
        p0 = plot(mprior_samples(ix,:),mprior_samples(iy,:),'.');
    end
    p1 = plot(minitial(ix),minitial(iy),'o','markersize',msizes,'markerfacecolor','k','markeredgecolor','w','linewidth',sthick);
    pP = plot(mprior(ix),mprior(iy),'o','markersize',msizes,'markerfacecolor','b','markeredgecolor','w','linewidth',sthick);
    pT = plot(mtarget(ix),mtarget(iy),'o','markersize',msizes,'markerfacecolor','r','markeredgecolor','w','linewidth',sthick);
    plot(xrec,yrec,'kV','markersize',msizer,'linewidth',rthick);
    for ii=1:ndata, text(xrec(ii),yrec(ii),num2str(ii),'fontsize',rfsize,'color','r',...
            'horizontalalignment','center','verticalalignment','middle'); end
    if ~isempty(mprior_samples)
        legend([p0(1) p1 pP pT],'Cprior sample','minitial','mprior','mtarget');
        x = mprior_samples(ix,:);
        y = mprior_samples(iy,:);
        iin = find( and( and( x >= ax0(1), x <= ax0(2)), and(y >= ax0(3), y <= ax0(4))) );
        title(sprintf('%i/%i (%.2f) prior samples inside the region shown',...
            length(iin),nsample,length(iin)/nsample));
    else
        legend([p1 pP pT],'minitial','mprior','mtarget');
    end
end

axis equal; axis(ax0); %grid on;
%set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
xlabel('X distance (km)'); ylabel('Y distance (km)');

%==========================================================================