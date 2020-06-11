%% test particle distribution
clear;clc
pvals = [0,0.01,0.1,0.5,1];%[0,1];
pstr = [{'0'},{'0pt01'},{'0pt1'},{'0pt5'},{'1'}];%[{'0'},{'1'}];
L = 100;
%%
fusedist = [];
partdist = [];
for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/nfusedist_pf%s.fuse.out',pstr{pc});
	% 	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/testnfuse.fuse.out',pstr{pc});
	[tvals,tmpfusedist,tmppartdist,nreg,regedges,ntrials] = readfusedist(filename);
	fusedist = cat(4,fusedist,tmpfusedist);
	partdist = cat(4,partdist,tmppartdist);
end
nsave = size(fusedist,1);
ntype = size(fusedist,3);
pdensity = squeeze(sum(partdist,2))/L/ntrials;
%% plot particle distribution
ptype = 1;
typestr = [{'lyso'},{'phago'}];
% figure
ssinds = nsave;
cmap = [0,0,1;1,0,0];
fc = 1;
% for tc = 1:20%nsave
	ct = 1;
	for pc = 1:length(pvals)
% 		disp(pc)
	subplot(2,3,ct)
% 		subplot(2,2,ct)
		for ptype = 1:ntype
			
		meanpartdist = smooth(mean(partdist(ssinds,:,ptype,pc),1)./ntrials,11);
% 			meanpartdist = smooth(mean(partdist(tc,:,ptype,pc),1)./ntrials,11);
% 		plot(regedges(1:end-1),meanpartdist,'Color',cmap(ptype,:),'DisplayName',...
			plot(regedges(1:end-1),meanpartdist,'Color',cmap(ptype,:),'DisplayName',...
				sprintf('%s: %0.3f/$\\mu$m',typestr{ptype},mean(pdensity(ssinds,ptype,pc))))
% 			sprintf('pfuse: %0.3f, particle density: %0.3f',pvals(pc),mean(pdensity(ssinds,ptype,pc))))
% 	histogram('BinEdges',regedges,'BinCounts',mean(partdist(ssinds,:,ptype,pc),1)./ntrials)
			plot_cleanup
			hold on
			xlabel('distance')
% 		if(ptype == 1)
% 			ylim([0,0.035])
% 			ylabel('number of lysosomes')
% 		else
% 			ylim([0,0.025])
% 			ylabel('number of autophagosomes')
% 		end
			ylabel('number of organelles')
			ylim([0,0.08])
% 			yticks([0,0.01,0.02])
			title(sprintf('fusion probability: %0.2f',pvals(pc)))
		lgd = legend;
% 		lgd.Location = 'EastOutside';
		lgd.Interpreter = 'latex';
		end
		
% 		pause(0.2)
		hold off
		ct = ct+1;
% 	end
% 	drawnow
% 	M(fc) = getframe(gcf);
% 	fc = fc+1;
end
%%
animation2movie(M,'npartdist_smooth.avi',16)


%%
hold on
for nc = 1:10:length(regedges(1:end-1))
	plot(regedges(nc)*[1,1],ylim,'k--')
end
hold off
%%
meanfusedist(isnan(meanfusedist)) = 0;
cvals = cfunc(meanfusedist,c0,rp,rl);
