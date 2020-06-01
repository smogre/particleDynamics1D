%% test distribution of nfuse
clear;clc
tic
filename = '/data/SimData/particledynamics1D/nfusefiles/testnfuse.fuse.out';
[tvals,fusedist,nreg,regedges,ntrials] = readfusedist(filename);
toc
%%
% for tc = 1:500
% 	histogram('BinEdges',regedges,'BinCounts',mean(fusedist(1:tc,:),1)./ntrials)
	histogram('BinEdges',regedges,'BinCounts',mean(fusedist(5000:end,:),1)./ntrials)
% 	plot(regedges(1:end-1),fusedist(tc,:))
	ylim([0,0.01])
	plot_cleanup
	xlabel('distance')
	ylabel('time- and ensemble-averaged fusion events')
	title(sprintf('time: %0.3f',tvals(tc)))
% 	drawnow
% 	M(tc) = getframe(gcf)
% 	pause(0.01)
% end
%%
animation2movie(M,'nfusedist.avi',24)


