%% test protein degradation
clear;clc;
protfile = '/data/SimData/particledynamics1D/outfiles/testprotdeg.out';
[tvals,cellprot,axprot,ntrials] = readprot(protfile);
%%
kdeg = 1;
clf
plot(tvals,axprot(:,1)./ntrials,'r')
hold on
plot(tvals,axprot(:,2)./ntrials,'b')
% plot(tvals,sum(axprot,2)./ntrials,'g')
% plot(tvals,exp(-kdeg*tvals),'ro--')
% plot(tvals,1-exp(-kdeg*tvals),'bo--')
hold off
plot_cleanup