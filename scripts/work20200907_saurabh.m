%% plot particle density/acidification/protein content histograms using snapshots
% get acidification directly from sims
clear;clc;
%%
pvals = [0,0.25,0.5,0.75,1];
pstr = ["0","0pt25","0pt50","0pt75","1"];
%%
tic
options = struct();
options.readlastsnap = 1;
lysopos = cell(length(pvals),1);
phpos = lysopos;
phacid = lysopos;
nfusedist = lysopos;
cellprot = lysopos;
axprot = lysopos;
protdist = lysopos;

for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('/data/SimData/particledynamics1D/snapfiles/nfusedist_pf%s.snap.out',pstr(pc));
	[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
	
	protfile = sprintf('/data/SimData/particledynamics1D/outfiles/nfusedist_pf%s.out',pstr(pc));
	[~,cellprot{pc},axprot{pc},~] = readprot(protfile);
	
	tmplysopos = [];
	tmpphpos = [];
	tmpphacid = [];
	tmpnfuse = [];
	tmpprot = [];
	
	for tc = 1:length(grouplist)
		pos = vertcat(grouplist(tc).particles(:).pos);
		acid = vertcat(grouplist(tc).particles(:).acidity);
		ptypes = vertcat(grouplist(tc).particles(:).type);
		nfuse = vertcat(grouplist(tc).particles(:).nfuse);
		prot = vertcat(grouplist(tc).particles(:).prot);
		
		tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
		tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
		tmpphacid = cat(1,tmpphacid,acid(ptypes==2));
		tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
		tmpprot = cat(1,tmpprot,prot(ptypes==2,:));
	end
	
	lysopos{pc} = tmplysopos;
	phpos{pc} = tmpphpos;
	phacid{pc} = tmpphacid;
	nfusedist{pc} = tmpnfuse;
	protdist{pc} = tmpprot;
end

% save('/data/SimData/particledynamics1D/matfiles/partdist_1e2trials_scanpf3vals_directacidity.mat')
toc
%%
% clear;clc
% load('/data/SimData/particledynamics1D/matfiles/partdist_scanpf_scankrev.mat')
load('/data/SimData/particledynamics1D/matfiles/partdist_1e2trials_scanpf5vals_newacidity.mat')
% load('/data/SimData/particledynamics1D/matfiles/partdist_1e2trials_scanpf5vals.mat')
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
% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_y0 = kpy/vy;
rho_y = @(beta,x) rho_y0*exp(-beta.*x);
% acidification from sims
c0 = 1;
afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/ry^3);
% afunc_sim = @(n,c0,rp,rl) n*c0.*rl^3/rp^3;
% number of fusion events:
Nfunc_th = @(beta,x) kpy/kpp*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,x) c0*Nfunc_th(beta,x)./(Nfunc_th(beta,x) + rp^3/ry^3);
% amount of degraded protein near cell body
k0 = 1;
kM = 0.01;
% degprot_th = @(c0,k0,kM,gam,beta) (1-exp(-(k0/vp)*afunc_th(c0,gam,beta,0)/(afunc_th(c0,gam,beta,0)+kM)));

%%
dt = 0.001; % simulation timestep
phfrate = 10; % framerate for phagosomes (binwidth)
lyfrate = 1.001;  % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
%% alternative binning
nreg = 370;
xp = linspace(0,1,nreg+2);
xl = xp;
%% plot organelle densities
% plotinds = [7,12,17];
% figure
plotinds = 1:length(pvals);

ct = 1;
for pc = plotinds
	
	if(length(pvals(plotinds))>1); subplot(1,2,ct); end
	nlyso = histcounts(lysopos{pc},xl);
	nphago = histcounts(phpos{pc},xp);
	
	xlvals = 0.5*(xl(1:end-1)+xl(2:end));
	% 	xlvals = xl(2:end);
	dxl = xlvals(2)-xlvals(1);
	
	xpvals = 0.5*(xp(1:end-1)+xp(2:end));
	% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);
	
	% 	xthvals = linspace(0,1,1e4);
	
	plot(xlvals,nlyso./ntrials./dxl,'bo','DisplayName','lyso density, sims')
	hold on
	plot(xlvals,rho_y(betvals(pc),xlvals),'b--','DisplayName','lyso density, an. model');
	plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
	plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.1]);
	% 	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
	ct = ct+1;
end
%% plot acidification
% figure
plotinds = 1:length(pvals);
cmap = colormapinterp([0,0,1;1,0,0],length(pvals(plotinds)));
ct = 1;
for pc = plotinds
	
	if(length(pvals)>1); subplot(2,3,ct); end
	
	[nphago,~,binind] = histcounts(phpos{pc},xp);
	% 	acvals = afunc_sim(nfusedist{pc,rc},c0,rp,ry);
	acvals = phacid{pc};
	
	% 	nfusehist = accumarray(binind,nfusedist{pc});
	achist = accumarray(binind,acvals);
	
	xpvals = 0.5*(xp(1:end-1)+xp(2:end));
	% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);
	
	plot(xpvals,achist./nphago','o','Color',cmap(ct,:),'DisplayName','acidification, sims')
	% 	plot(xpvals,achist./nphago','o','Color',cmap(ct,:),'DisplayName','acidification, sims')
	% 	plot(xpvals,afunc_sim(nfusehist,c0,rp,ry),'o','Color',cmap(pc,:),...
	% 																																'DisplayName','acidification, sims')
	hold on
	plot(xpvals,afunc_th(betvals(pc),xpvals),'^--','Color',cmap(ct,:),...
		'DisplayName','acidification, an. model');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.1]);
	% 	ylim([0,10]);
	xlabel('distance')
	ylabel('acidification')
	ct = ct+1;
end
%% plot distribution of degraded protein content
% figure
% plotinds = [7,12,17];
plotinds = 1:length(pvals);
% cmap = colormapinterp([0,0,1;1,0,0],length(runlenvals));
cmap = colormapinterp([0,0,1;1,0,0],length(plotinds));
ct = 1;
for pc = 1:length(plotinds)
	pind = plotinds(pc);
	pfval = pvals(pind);
	if(length(pvals)>1); subplot(1,2,ct); end
	% 	if(length(plotinds)>1); subplot(1,3,ct); end
	
	[nphago,~,binind] = histcounts(phpos{pind},xp);
	undegprot = accumarray(binind,protdist{pind}(:,1));
	degprot = accumarray(binind,protdist{pind}(:,2));
	
	xpvals = 0.5*(xp(1:end-1)+xp(2:end));
	
	plot(xpvals,degprot./nphago','o-','Color',cmap(ct,:),'DisplayName',...
		sprintf('pf: %0.3f, sims',pfval))
	hold on
	% 		[xth,degprot_th] = ode45(@(x,y) degProtdist(x,y,k0,kM,gamvals(pind),betvals(pind),vp),...
	% 																																							linspace(1,0,1e3),0);
	[xth,degprot_th] = ode45(@(x,y) degProtdist(x,y,k0,kM,gamvals(pind),betvals(pind),vp),...
		fliplr(xpvals),0);
	plot(xth,degprot_th,'^--','Color',cmap(ct,:),'DisplayName',sprintf('pf: %0.3f, sims',pfval))
	hold off
	plot_cleanup
	title(sprintf('pf: %0.3f',pfval))
	lgd = legend;
	lgd.Interpreter = 'latex';
	xlim([0,0.1]);
	% 	ylim([0,10]);
	xlabel('distance')
	ylabel('degraded protein content')
	ct = ct+1;
end
%% get cell protein and axonal protein amounts
num_ssinds = 100;
undeg_cellprot_ss = zeros(num_ssinds+1,length(pvals));
deg_cellprot_ss = undeg_cellprot_ss;
undeg_axprot_ss = undeg_cellprot_ss;
deg_axprot_ss = undeg_cellprot_ss;
deg_protflux_ss = undeg_cellprot_ss;
undeg_protflux_ss = undeg_cellprot_ss;
for pc = 1:length(pvals)
	undeg_cellprot_ss(:,pc) = cellprot{pc}((end-num_ssinds:end),1)/ntrials;
	deg_cellprot_ss(:,pc) = cellprot{pc}((end-num_ssinds:end),2)/ntrials;
	undeg_axprot_ss(:,pc) = axprot{pc}((end-num_ssinds:end),1)/ntrials;
	deg_axprot_ss(:,pc) = axprot{pc}((end-num_ssinds:end),2)/ntrials;
	deg_protflux_ss(:,pc) = diff(cellprot{pc}((end-num_ssinds-1):end,2))/dt/ntrials;
	undeg_protflux_ss(:,pc) = diff(cellprot{pc}((end-num_ssinds-1):end,1))/dt/ntrials;
end
%%
% figure
maxconvrate = 1;
uprot = mean(deg_axprot_ss+undeg_protflux_ss/maxconvrate,1);
mean_deg_protflux_ss = squeeze(mean(deg_protflux_ss));
std_deg_protflux_ss = squeeze(std(deg_protflux_ss));
%%
% figure
clf
cmap = colormapinterp([0,0,1;1,0,0],length(runlenvals));
for rc = 1:length(runlenvals)
	% 	errorbar(pvals,mean_deg_protflux_ss(:,rc),std_deg_protflux_ss(:,rc),'o-','Color',cmap(rc,:),...
	% 					'LineWidth',2);%,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	errorbar(pvals,mean_deg_protflux_ss,std_deg_protflux_ss,'o-','Color',cmap(rc,:),...
		'LineWidth',2);%,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	hold on
end
hold off
plot_cleanup
xlabel('pf')
ylabel('flux of degraded protein')
lgd = legend;
lgd.Interpreter = 'latex';

%%
mean_deg_cellprot_ss = squeeze(mean(deg_cellprot_ss));
mean_undeg_cellprot_ss = squeeze(mean(undeg_cellprot_ss));
mean_deg_axprot_ss = squeeze(mean(deg_axprot_ss));
mean_undeg_axprot_ss = squeeze(mean(undeg_axprot_ss));
std_deg_cellprot_ss = squeeze(std(deg_cellprot_ss));
std_undeg_cellprot_ss = squeeze(std(undeg_cellprot_ss));
std_deg_axprot_ss = squeeze(std(deg_axprot_ss));
std_undeg_axprot_ss = squeeze(std(undeg_axprot_ss));

for rc = 1:length(runlenvals)
	subplot(1,2,1)
	% 	errorbar(pvals,mean_undeg_cellprot_ss(:,rc),std_undeg_cellprot_ss(:,rc),'o-','Color',cmap(rc,:),...
	% 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',vy/runlenvals(rc)))
	errorbar(pvals,mean_undeg_cellprot_ss,std_undeg_cellprot_ss,'o-','Color',cmap(rc,:),...
		'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	hold on
	plot_cleanup
	xlabel('pf')
	ylabel('protein content')
	title(sprintf('undegraded cell prot\n'))
	lgd = legend;
	lgd.Interpreter = 'latex';
	
	subplot(1,2,2)
	% 	errorbar(pvals,mean_undeg_axprot_ss(:,rc),std_undeg_axprot_ss(:,rc),'o-','Color',cmap(rc,:),...
	% 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',vy/runlenvals(rc)))
	errorbar(pvals,mean_undeg_axprot_ss,std_undeg_axprot_ss,'o-','Color',cmap(rc,:),...
		'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	hold on
	plot_cleanup
	xlabel('pf')
	title(sprintf('undegraded ax prot\n'))
	lgd = legend;
	lgd.Interpreter = 'latex';
	
	% 	subplot(1,2,1)
	% % 	errorbar(pvals,mean_deg_cellprot_ss(:,rc),std_deg_cellprot_ss(:,rc),'o-','Color',cmap(rc,:),...
	% % 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	% 	errorbar(pvals,mean_deg_cellprot_ss,std_deg_cellprot_ss,'o-','Color',cmap(rc,:),...
	% 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	% 	hold on
	% 	plot_cleanup
	% 	xlabel('pf')
	% 	title('degraded cell prot')
	% 	lgd = legend;
	% 	lgd.Interpreter = 'latex';
	%
	% 	subplot(1,2,2)
	% % 	errorbar(pvals,mean_deg_axprot_ss(:,rc),std_deg_axprot_ss(:,rc),'o-','Color',cmap(rc,:),...
	% % 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	% 	errorbar(pvals,mean_deg_axprot_ss,std_deg_axprot_ss,'o-','Color',cmap(rc,:),...
	% 												'LineWidth',2,'DisplayName',sprintf('runlen: %0.3f',runlenvals(rc)))
	% 	hold on
	% 	plot_cleanup
	% 	xlabel('pf')
	% 	title('degraded ax prot')
	% 	lgd = legend;
	% 	lgd.Interpreter = 'latex';
end
%%
[x,y] = ode45(@(x,y) degProtdist(x,y,k0,kM,gamvals(pc),betvals(pc),vp),linspace(1,0,1e3),0);
%%
function dG = degProtdist(x,G,k0,kM,gamvals,betvals,vp)
dG = -degrate(k0,kM,acidification(gamvals,betvals,x))*(1-G)./vp;
end
%%
function kd = degrate(k0,kM,aval)
kd = k0.*aval./(aval+kM);
end
%%
function aval = acidification(gam,bet,x)
aval = 1-exp(gam/bet.*(exp(-bet)-exp(-bet.*x)));
end