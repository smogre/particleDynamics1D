%% script to test reading in protein content
tic
clear;clc
fname = '/data/SimData/particledynamics1D/outfiles/test.out';

[tvals,cellprot,axprot,ntrials] = readprot(fname);
cellprot = cellprot/ntrials;
axprot = axprot/ntrials;
nprot = size(cellprot,2);
protflux = diff(cellprot)/(tvals(2)-tvals(1));
dtvals = tvals(2:end)-tvals(1);
toc
%% plot protein content
cmap = jet(nprot);
hold on
for nc = 1:nprot
	plot(tvals,cellprot(:,nc),'Color',cmap(nc,:),'DisplayName',sprintf('Cumulative cell protein, type %d',nc));
	hold on
	plot(tvals,axprot(:,nc),'--','Color',cmap(nc,:),'DisplayName',sprintf('Axonal protein,type %d',nc));
end
plot_cleanup
xlabel('time')
ylabel('protein content')
% title('type 1: undegraded, type 2: degraded, $k_\mathrm{fuse}$ = 0.1');
title('type 1: undegraded, type 2: degraded, $k_\mathrm{fuse}$ = 0.1');
% ylim([0,70])
legend show
hold off
%% plot protein flux

cmap = jet(nprot);
for nc = 1:nprot
	plot(dtvals,protflux(:,nc),'Color',cmap(nc,:),'DisplayName',sprintf('Protein type %d',nc));
	hold on
end
plot_cleanup
xlabel('time')
ylabel('protein flux')
title(sprintf('$k_\\mathrm{fuse}$ = 0.1/s, $p1_\\mathrm{mean}$ = %0.3f, $p2_\\mathrm{mean}$ = %0.3f',mean(protflux(51:end,:))))
hold off
%% plot total protein flux
cmap = jet(nprot);
plot(dtvals,sum(protflux,2));
plot_cleanup
xlabel('time')
ylabel('total protein flux')
title(sprintf('no fusion, mean steady state value: %0.3f',mean(sum(protflux(floor(length(protflux)/2):end,:),2))))
hold off
%% plot total unprocessed protein
maxconvrate = 5;
uprot = sum(axprot,2)+[0;protflux(:,2)]/maxconvrate;
cmap = jet(nprot);
plot(tvals,uprot);
plot_cleanup
xlabel('time')
ylabel('total unprocessed protein')
title(sprintf('$k_\\mathrm{fuse}$ = 0.1/s, $up_\\mathrm{mean}$ = %0.3f, maxnfuse = 5',mean(uprot(51:end))))
% title(sprintf('no fusion, $up_\\mathrm{mean}$ = %0.3f, maxnfuse = 5',mean(uprot(51:end))))
hold off