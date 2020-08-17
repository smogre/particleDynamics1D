%% plot distribution of degraded protein
clear;clc;
pvals = [0,0.25,0.5,0.75,1];
pstr = [{'0'},{'0pt25'},{'0pt50'},{'0pt75'},{'1'}];

%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readlastsnap = 0;
lysopos = cell(length(pvals),1);
phpos = cell(length(pvals),1);
nfusedist = cell(length(pvals),1);
cellprot = cell(length(pvals),1);
axprot = cell(length(pvals),1);
protdist = cell(length(pvals),1);

for pc = 1:length(pvals)
	disp(pc)
	snapfile = sprintf('/data/SimData/particledynamics1D/snapfiles/nfusedist_pf%s.snap.out',pstr{pc});
	% 	filename = '/data/SimData/particledynamics1D/snapfiles/testnfuse.snap.out';
	[grouplist,tvals,domlen,ntrials] = readsnapshot(snapfile,options);
	protfile = sprintf('/data/SimData/particledynamics1D/outfiles/nfusedist_pf%s.out',pstr{pc});
	[~,cellprot{pc},axprot{pc},~] = readprot(protfile);
	
	tmplysopos = [];
	tmpphpos = [];
	tmpnfuse = [];
	tmpprot = [];
	
	for tc = 1:length(grouplist)
		pos = vertcat(grouplist(tc).particles(:).pos);
		ptypes = vertcat(grouplist(tc).particles(:).type);
		nfuse = vertcat(grouplist(tc).particles(:).nfuse);
		prot = vertcat(grouplist(tc).particles(:).prot);
		
		tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
		tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
		tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
		tmpprot = cat(1,tmpprot,prot(ptypes==2,:));
	end
	lysopos{pc} = tmplysopos;
	phpos{pc} = tmpphpos;
	nfusedist{pc} = tmpnfuse;
	protdist{pc} = tmpprot;
end
toc
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
gamvals = kpy*ell/vp/vy*ry^3/(ry^3+rp^3)*kfvals;
% acidicication
afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_y0 = kpy/vy;
rho_y = @(beta,x) rho_y0*exp(-beta.*x);
% acidification
c0 = 1;
% afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
afunc_sim = @(n,c0,rp,rl) n*c0.*rl^3/rp^3;
%%
dt = 0.001; % simulation timestep
phfrate = 10; % framerate for phagosomes (binwidth)
lyfrate = 1.01; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
%% alternative binning
nreg = 100;
xp = linspace(0,1,nreg+2);
xl = xp;
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
% 	plot(xlvals,rho_y(betvals(pc),xlvals),'b--','DisplayName','lyso density, an. model');
	plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
% 	plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
% 	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end