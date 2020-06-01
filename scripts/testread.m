clear;clc;
tic
fname = '/data/SimData/particledynamics1D/snapfiles/test.snap.out';
options = struct();
options.progstep = 100;
% options.maxpart = 10;
% options.onlyread = 2;
[grouplist,tvals,domlen,ntrials] = readsnapshot(fname,options);
toc