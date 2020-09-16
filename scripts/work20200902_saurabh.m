%% test protein degradation
clear;clc;
protfile = '../testprotdeg.out';
[tvals,cellprot,axprot,ntrials] = readprot(protfile);
%%
npart = 100;
kdeg = 1-tvals;
[x,G] = ode45(@(x,G) degProtdist(x,G),flipud(tvals),0);
%%
clf
plot(tvals,axprot(:,1)./ntrials./npart,'b','DisplayName','undeg')
hold on
plot(tvals,axprot(:,2)./ntrials./npart,'r','DisplayName','deg')
% plot(tvals,sum(axprot,2)./ntrials,'g')
% plot(tvals,exp(-kdeg.*tvals),'ro--')
% plot(tvals,1-exp(-kdeg.*tvals),'bo--')
plot(x,flipud(G),'ro--')
hold off
plot_cleanup
lgd = legend;
%%
function dy = degProtdist(x,y)
	dy = -(1-x)*(1-y);
end

