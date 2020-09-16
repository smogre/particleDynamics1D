%% predict number of fusions
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

pvals = 0.25:0.25:1;

% fusion rates
kfvals = -log(1-pvals)/T;

% lyso production rate
kpy = 17.9*10;
% phago production rate
kpp = 1.6*10;
% encounter length
ell = 2*(rp+ry);
% beta
betvals = kpp*ell/vp/vy*kfvals;
% gamma
gamvalsN = kpy*ell/vp/vy*kfvals;
% gamma
gamvals = kpy*ell/vp/vy*ry^3/(ry^3+rp^3)*kfvals;

N0 = gamvalsN./betvals.*(1-exp(-betvals*Lhat))


A0 = 1 - exp(gamvals./betvals.*(exp(-betvals*Lhat) - 1))

%% plot particle density histograms using snapshots
clear;clc;
pvals = [0,0.25,0.5,0.75,1];
pstr = [{'0'},{'0pt25'},{'0pt50'},{'0pt75'},{'1'}];
%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readsinglesnap = 1;
lysopos = cell(length(pvals),1);
phpos = cell(length(pvals),1);
nfusedist = cell(length(pvals),1);
for pc = 2:3%1:length(pvals)
	disp(pc)
	filename = sprintf('../test/nfusedist_pf%s.snap.out',pstr{pc});
% 	filename = '/data/SimData/particledynamics1D/snapfiles/testnfuse.snap.out';
	[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
	tmplysopos = [];
	tmpphpos = [];
	tmpnfuse = [];
	for tc = 1:length(grouplist)
		pos = vertcat(grouplist(tc).particles(:).pos);
		ptypes = vertcat(grouplist(tc).particles(:).type);
		nfuse = vertcat(grouplist(tc).particles(:).nfuse);
		tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
		tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
		tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
	end
	lysopos{pc} = tmplysopos;
	phpos{pc} = tmpphpos;
	nfusedist{pc} = tmpnfuse; % number of fusions for each phagosome in all trials
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
%afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
afunc_sim = @(n,c0,rp,rl) n*c0*rl^3./(rp^3+rl^3);
%%
dt = 0.001; % simulation timestep
phfrate = 1.0; % framerate for phagosomes (binwidth)
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
	plot(xlvals,rho_y(betvals(pc),xlvals),'b--','DisplayName','lyso density, an. model');
	plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
	plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
	hold off
	%plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end




%% alternative binning
nreg = 100;
xp = linspace(0,1,nreg+2);
xl = xp;
%% plot organelle densities
for pc = 2:3%1:length(pvals)
	
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
	%plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end
%% plot acidification
% figure
cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
for pc = 2:3%1:length(pvals)
	
	if(length(pvals)>1); subplot(2,3,pc); end
	
	[nphago,~,binind] = histcounts(phpos{pc},xp);
	nfusehist = accumarray(binind,nfusedist{pc})./nphago';

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xpvals,afunc_sim(nfusehist,c0,rp,ry),'o','Color',cmap(pc,:),...
        'DisplayName','acidification, sims')
	hold on
	plot(xpvals,afunc_th(c0,gamvals(pc),betvals(pc),xpvals),'--','Color',cmap(pc,:),...
        'DisplayName','acidification, an. model');
	hold off
	%plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
% 	ylim([0,10]);
	xlabel('distance')
	ylabel('organelle density')
	
end

