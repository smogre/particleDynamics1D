%% test particle distribution
clear;clc
%%
pvals = [0,0.25,0.5,0.75,1];
pstr = [{'0'},{'0pt25'},{'0pt50'},{'0pt75'},{'1'}];
L = 100;
fusedist = [];
partdist = [];
for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/nfusedist_pf%s.fuse.out',pstr{pc});
	[tvals,tmpfusedist,tmppartdist,nreg,regedges,ntrials] = readfusedist(filename);
	fusedist = cat(4,fusedist,tmpfusedist);
	partdist = cat(4,partdist,tmppartdist);
end
dx = regedges(2)-regedges(1);
%%
pvals = 0.25;
pstr = {'0pt25'};
L = 100;
fname = '/data/SimData/particledynamics1D/nfusefiles/testnfuse.fuse.out';
[tvals,fusedist,partdist,nreg,regedges,ntrials] = readfusedist(fname);
dx = regedges(2)-regedges(1);
%%
nsave = size(fusedist,1);
ntype = size(fusedist,3);
pdensity = partdist./ntrials./dx;
meanpdensity = squeeze(mean(pdensity,2));
%% plot particle distribution
typestr = [{'lyso'},{'phago'}];
% figure
ssinds = nsave-100:nsave;
cmap = [0,0,1;1,0,0];
xvals = regedges(1:end-1);

ct = 1;
for pc = 1:length(pvals)
	
	subplot(2,3,ct)
	
	for ptype = 1:ntype
		
		
		pdensity_ss = mean(pdensity(ssinds,:,ptype,pc),1);
		
		plot(xvals,pdensity_ss,'o','Color',cmap(ptype,:),'DisplayName',...
			sprintf('%s: %0.3f/L',typestr{ptype},mean(meanpdensity(ssinds,ptype,pc))))
		
		hold on
		if(ptype==1)
			plot(xvals,rho_y(betvals(pc),xvals),'--','Color',cmap(ptype,:),'DisplayName',...
				sprintf('%s density, an. model',typestr{ptype}))
		else
			plot(xvals,rho_p(xvals),'--','Color',cmap(ptype,:),'DisplayName',...
				sprintf('%s density, an. model',typestr{ptype}))
		end
		plot_cleanup
		xlabel('distance')
		xlim([0,0.2]);
		ylabel('organelle density')
		
		title(sprintf('fusion probability: %0.2f',pvals(pc)))
		% 			lgd = legend;
		% 			lgd.Interpreter = 'latex';
	end
	
	
	hold off
	ct = ct+1;
end
%%
animation2movie(M,'npartdist_smooth.avi',16)
%% analytical model predictions
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
kpy = 17.9*100;
% phago production rate
kpp = 1.6*100;
% encounter length
ell = 2*(rp+ry);
% lyso acidification
c0 = 1;
%%
% beta
betvals = kpp*ell/vp/vy*kfvals;
% gamma
gamvals = kpy*ell/vp/vy*vy/(vy+vp)*kfvals;
% acidicication
afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_y0 = kpy/vy;
rho_y = @(beta,x) rho_y0*exp(-beta.*x);
%%
hold on
for nc = 1:10:length(regedges(1:end-1))
	plot(regedges(nc)*[1,1],ylim,'k--')
end
hold off
%%
meanfusedist(isnan(meanfusedist)) = 0;
cvals = cfunc(meanfusedist,c0,rp,rl);
