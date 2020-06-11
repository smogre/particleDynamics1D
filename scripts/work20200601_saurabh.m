%% test distribution of concentration
clear;clc
pvals = [0.01,0.1,0.5];[0,0.01,0.1,0.5,1];%[0,1];
pstr = [{'0pt01'},{'0pt1'},{'0pt5'}];%[{'0'},{'1'}];
cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
%%
c0 = 1;
rp = 400;
rl = 100;
cfunc = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
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
%%
meanfusedist = fusedist./partdist;
meanfusedist(isnan(meanfusedist)) = 0;
cvals = cfunc(meanfusedist,c0,rp,rl);
%%
L = 100;
pdensity = squeeze(sum(partdist,2))/L/ntrials;
%%
% cfunc_th = 
%%
ptype = 1;
typestr = [{'lyso'},{'phago'}];
figure
ssinds = nsave;
ct = 1;
for pc = 1:length(pvals)
	disp(pc)

	meancvals = smooth(mean(cvals(ssinds,:,2,pc),1),11);
% 	semvals = std(cvals(ssinds,:),1)./sqrt(length(ssinds));
% 	curve1 = meancvals+semvals;
% 	curve2 = meancvals-semvals;
% 	x2 = [regedges(1:end-1),fliplr(regedges(1:end-1))];
% 	inBetween = [curve1,fliplr(curve2)];
	plot(regedges(1:end-1),meancvals,'Color',cmap(pc,:),'DisplayName',...
																									sprintf('fusion probability: %0.3f',pvals(pc)))
% 	stdshade(cvals(ssinds,:),0.4,regedges(1:end-1))
% 	hold on
% 	fill(x2,inBetween,cmap(pc,:),'FaceAlpha',0.4,'linestyle','none','DisplayName',...
%																							sprintf('SEM, fusion probability: %0.3f',pvals(pc)))
	
% 	plot(regedges(1:end-1),mean(partdist(ssinds,:,ptype),1)./ntrials,'DisplayName',...
% 																																sprintf('pfuse: %0.3f',pvals(pc)))
	
	hold on
% 	plot(regedges(1:end-1),afunc_th(c0,gamvals(pc),betvals(pc),regedges(1:end-1)),'--',...
% 			'Color',cmap(pc,:));
end
% 	ylim([0,0.001])
hold off
plot_cleanup
xlabel('distance')
% ylabel('time- and ensemble-averaged fusion events')
ylabel('relative autophagosome acidification')
% ylabel('number of autophagosomes')
%  	title(sprintf('time: %0.3f',tvals(tc)))
% 	title(sprintf('pfuse: %0.3f',pfuse))
% 	drawnow
% 	M(tc) = getframe(gcf)
% 	pause(0.01)
% end
legend show
%%
% animation2movie(M,'nfusedist.avi',24)
%%
% domain length in microns
L = 100;
% scaled domain length
Lhat = 1;
% time unit
tau = L/0.75;
% phago velocity
vp = 1;	% 0.75 um/s
% lyso velocity
vy = 2.67;	% 2um/s
% phago radius
rp = 0.4/L; % 400nm
% lyso radius
ry = 0.1/L; % 100nm
% encounter time
T = 2*(rp+ry)/(vp+vy);
% fusion rates
kfvals = -log(1-pvals)/T;
% lyso production rate
kpy = 17.9;
% phago production rate
kpp = 1.6;
% encounter length
ell = (rp+ry);
% lyso acidification
c0 = 1;
% beta
betvals = kpp*ell/vp/vy*kfvals;
% gamma
gamvals = kpy*ell/vp/vy*vy/(vy+vp)*kfvals;
% acidicication
afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));

