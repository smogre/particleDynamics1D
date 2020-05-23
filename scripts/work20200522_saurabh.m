%% steady state protein flux vs krev
% kfuse = 0.1/s
% vl = 4vp

%% test for a single outfile
clear;clc
tic
nk = 100;
nsave = 1000;
kswitch = logspace(-2,4,nk);
cellprot = zeros(nsave,2,nk);
axprot = zeros(nsave,2,nk);
for fc = 1:nk
	disp(fc)
	fname = sprintf('/data/Condor/particledynamics1D/outfiles/scankfuse.%d.out',fc-1);
	[tvals,tmpcell,tmpax,ntrials] = readprot(fname);
	cellprot(:,:,fc) = tmpcell/ntrials;
	axprot(:,:,fc) = tmpax/ntrials;
end
protflux = diff(cellprot)/(tvals(2)-tvals(1));
nprot = size(cellprot,2);
dtvals = tvals(2:end)-tvals(1);
maxconvrate = 5;
uprot = squeeze(sum(axprot,2))+[zeros(1,nk);squeeze(protflux(:,1,:))]/maxconvrate;
toc
%% get steady state values
ssinds = floor(nsave/2):nsave;
uprot_ss = mean(uprot(ssinds,:));
axprot_ss = squeeze(mean(axprot(ssinds,:,:)))';
cellprot_ss = squeeze(mean(cellprot(ssinds,:,:)))';
protflux_ss = squeeze(mean(protflux(ssinds-1,:,:)))';
%% plot protein content/flux
cmap = jet(nprot);
for nc = 1:nprot
	semilogx(kswitch,cellprot_ss(:,nc),'Color',cmap(nc,:),'DisplayName',sprintf('Protein type %d',nc));
	hold on
	semilogx(kswitch,axprot_ss(:,nc),'--','Color',cmap(nc,:),'DisplayName',sprintf('Protein type %d',nc));
end
plot_cleanup
xlabel('time')
ylabel('protein flux')
% title(sprintf('$k_\\mathrm{fuse}$ = 0.1/s, $p1_\\mathrm{mean}$ = %0.3f, $p2_\\mathrm{mean}$ = %0.3f',mean(protflux(51:end,:))))
hold off