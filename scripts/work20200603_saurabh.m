%% 
trc = 1:ntrials;
ptype = 2;
for tc = length(tvals)
	for ptype = 1:ntype
	plot(regedges(1:end-1),mean(partdist(tc,:,trc,ptype),3))
	ylim([0,0.03])
	hold on
	
	end
	hold off
	plot_cleanup
% 	drawnow
% 	pause(0.01)
end