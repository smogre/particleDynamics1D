%% plot particle density histograms using fusion distributions
clear;clc;
pvals = [0,0.25,0.5,0.75,1];
pstr = [{'0'},{'0pt25'},{'0pt50'},{'0pt75'},{'1'}];
%%
tic
cmap = [0,0,1;0,0,1];
fusedist = [];
partdist = [];
for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/nfusedist_pf%s.fuse.out',pstr{pc});
% 	filename = '/data/SimData/particledynamics1D/snapfiles/testnfuse.snap.out';
	[tvals,tmpfusedist,tmppartdist,nreg,regedges,ntrials] = readfusedist(filename);
	fusedist = cat(4,fusedist,tmpfusedist);
	partdist = cat(4,partdist,tmppartdist);
end
dx = regedges(2)-regedges(1);
toc
%%
nsave = size(fusedist,1);
ntype = size(fusedist,3);
pdensity = partdist./ntrials./dx;
meanpdensity = squeeze(mean(pdensity,2));
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
% acidification
c0 = 1;
afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
%%
dt = 0.001; % simulation timestep
phfrate = 10; % framerate for phagosomes (binwidth)
lyfrate = 1.01; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
%% plot organelle densities
for pc = 1:length(pvals)
	
	if(length(pvals)>1); subplot(2,3,pc); end
	
	nlyso = histcounts(lysopos{pc},xl);
	nphago = histcounts(phpos{pc},xp);
	
	xlvals = 0.5*(xl(1:end-1)+xl(2:end)); 
% 	xlvals = xl(2:end);
	dxl = xlvals(2)-xlvals(1);
	
	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xlvals,nlyso./ntrials./dxl,'bo','DisplayName','lyso density, sims')
	hold on
	plot(xlvals,rho_y(betvals(pc),xlvals),'b--','DisplayName','lyso density, an. model');
	plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
	plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end
%% plot acidification
% figure
for pc = 1:length(pvals)
	
	if(length(pvals)>1); subplot(2,3,pc); end
	
	[~,~,binind] = histcounts(phpos{pc},xp);
	nfusehist = accumarray(binind,nfusedist{pc});

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xpvals,afunc_sim(nfusehist,c0,rp,ry),'ro','DisplayName','acidification, sims')
	hold on
	plot(xpvals,afunc_th(c0,gamvals(pc),betvals(pc),xpvals),'k--','DisplayName','acidification, an. model');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
% 	ylim([0,10]);
	xlabel('distance')
	ylabel('organelle density')
	
end

