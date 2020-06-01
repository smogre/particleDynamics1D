%% test reading snapshots, dependence of protein flux on number of fusion events/ reversal rate
clear;clc;
tic
fname = '/data/SimData/particledynamics1D/snapfiles/test.snap.out';
options = struct();
options.progstep = 10;
options.maxpart = 100;
% options.onlyread = 2;
[grouplist,tvals,domlen,ntrials] = readsnapshot(fname,options);
% [grouplist_old,tvals,domlen,ntrials] = readsnapshot_old(fname,options);
toc
%%
nfuse = [];
allnpart = zeros(length(tvals),ntrials);
for gc = 1:length(grouplist)
	for pc = 1:length(grouplist(gc).particles)
		if(grouplist(gc).particles(pc).type==2)
			nfuse = cat(1,nfuse,grouplist(gc).particles(pc).nfuse);
		end
	end
	allnpart(:,gc) = grouplist(gc).npart;
end
%%
histogram(nfuse)
%%
mean(nfuse)
%% 
plot(tvals,mean(allnpart,2))
plot_cleanup