%% steady state protein flux vs krev
% kfuse = 0.1/s
% vl = 4vp

%% test for a single outfile
clear;clc
tic
nkf = 4;
nks = 100;
nsave = 10000;
kswitch = logspace(-2,4,nks);
kfvals = [0,0.1,1,10];
kfstr = [{'0'},{'1'},{'10'},{'100'}];
cellprot = zeros(nsave,2,nks,nkf);
axprot = zeros(nsave,2,nks,nkf);
for fc = 1:nkf
	for sc = 1:nks
		disp([fc sc])
		% 	fname = sprintf('/data/Condor/particledynamics1D/outfiles/scankfuse.%d.out',fc-1);
		fname = sprintf('/data/Condor/particledynamics1D/scankrev_kf%s.%d.out',kfstr{fc},(sc-1));
		[tvals,tmpcell,tmpax,ntrials] = readprot(fname);
		cellprot(:,:,sc,fc) = tmpcell/ntrials;
		axprot(:,:,sc,fc) = tmpax/ntrials;
	end
end
protflux = diff(cellprot)/(tvals(2)-tvals(1));
nprot = size(cellprot,2);
dtvals = tvals(2:end)-tvals(1);
maxconvrate = 5;
uprot = squeeze(sum(axprot,2))+[zeros(1,nks,nkf);squeeze(protflux(:,1,:,:))]/maxconvrate;
% uprot = squeeze(sum(axprot,2))+[zeros(1,nks,nkf);squeeze(protflux(:,1,:,:))]/maxconvrate;
toc
%% get steady state values
ssinds = floor(nsave/2):nsave;
uprot_ss = squeeze(mean(uprot(ssinds,:,:),1));
axprot_ss = permute(squeeze(mean(axprot(ssinds,:,:,:),1)),[2,1,3]);
cellprot_ss = permute(squeeze(mean(cellprot(ssinds,:,:,:),1)),[2,1,3]);
protflux_ss = permute(squeeze(mean(protflux(ssinds-1,:,:,:),1)),[2,1,3]);
%% plot protein content/flux
% cmap = jet(nprot);
dstr = [{'-'},{'--'},{'-.'}];
cmap = jet(nkf);
plotinds = ssinds;
% plotinds = 5000:5100;
% hold on
ct = 1;
% for kc = 50%[1,50,100]
for fc = 1:nkf
% 	for nc = 1:nprot
% 		semilogx(kswitch,cellprot_ss(:,nc),'Color',cmap(nc,:),...
% 																										'DisplayName',sprintf('Protein type %d',nc))
% 		semilogx(kswitch,squeeze(cellprot_ss(:,nc,fc)),dstr{nc},'Color',cmap(fc,:),...
% 																			'DisplayName',sprintf('Protein type %d, kf: %s',nc,kfstr{fc}))
% 			semilogx(kswitch,squeeze(protflux_ss(:,nc,fc)),dstr{nc},'Color',cmap(fc,:),...
% 																			'DisplayName',sprintf('Protein type %d, kf: %s',nc,kfstr{fc}))
% 			plot(tvals(plotinds),uprot(plotinds,kc,fc),dstr{ct},'Color',cmap(fc,:),...
% 															'DisplayName',sprintf('kf: %s, kswitch: %0.3f',kfstr{fc},kswitch(kc)))
		semilogx(kswitch,uprot_ss(:,fc),'Color',cmap(fc,:),...
																			'DisplayName',sprintf('kf: %s',kfstr{fc}))
		hold on
% 		semilogx(kswitch,uprot_ss,'Color',cmap(nc,:),...
% 																		'DisplayName',sprintf('Protein type %d, kf: %s',nc,kfstr{fc}))
% 	end
end
% ct = ct+1;
% end
plot_cleanup
legend show
xlabel('kswitch')
% xlabel('time')
ylabel('total cell protein')
% title(sprintf('$k_\\mathrm{fuse}$ = 0.1/s, $p1_\\mathrm{mean}$ = %0.3f, $p2_\\mathrm{mean}$ = %0.3f',mean(protflux(51:end,:))))
hold off