%% test effect of fusion rates
kfvals = [0,1,10];
flist = dir('/data/SimData/particledynamics1D/outfiles/*fuse*');
tic
nk = 3;
nsave = 100;
kswitch = logspace(-2,4,nk);
cellprot = zeros(nsave,2,nk);
axprot = zeros(nsave,2,nk);
for fc = 1:nk
	disp(fc)
	fname = [flist(fc).folder,'/',flist(fc).name];
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
	plot(kfvals,cellprot_ss(:,nc),'o-','Color',cmap(nc,:),'DisplayName',...
																									sprintf('Steady state cell protein, type %d',nc));
% 	plot(kfvals,cellprot_ss(:,nc),'o-','Color',cmap(nc,:),'DisplayName',...
% 																									sprintf('Steady state protein flux, type %d',nc));
	hold on
	plot(kfvals,axprot_ss(:,nc),'o--','Color',cmap(nc,:),'DisplayName',...
																								sprintf('Steady state axonal protein, type %d',nc));
end
plot_cleanup
xlabel('fusion rate')
ylabel('steady state protein content')
% title(sprintf('$k_\\mathrm{fuse}$ = 0.1/s, $p1_\\mathrm{mean}$ = %0.3f, $p2_\\mathrm{mean}$ = %0.3f',mean(protflux(51:end,:))))
title('effect of fusion rate')
legend show
hold off
%% plot total cell protein
semilogx(kfvals,uprot_ss,'o-')
plot_cleanup
xlabel('fusion rate')
ylabel('total cell protein')
title('effect of fusion rate')