%% test reversal of particles
clear;clc;
tic
fname = '/data/SimData/particledynamics1D/snapfiles/testrev.snap.out';
options = struct();
% options.readlastsnap = 1;
[grouplist,tvals,domlen,ntrials] = readsnapshot(fname,options);
toc