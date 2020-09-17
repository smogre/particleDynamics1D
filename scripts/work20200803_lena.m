%% plot particle density histograms using fusion distributions
clear;clc;
pvals = [0.25,0.5,0.75];
pstr = {'0pt25'}%,'0pt5','0pt75'};

%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readlastsnap = 1;
lysopos = cell(length(pvals),1);
phpos = cell(length(pvals),1);
phprots = cell(length(pvals),1);
nfusedist = cell(length(pvals),1);
for pc = 1%:length(pvals)
	disp(pc)
	filename = sprintf('../test/snapfiles/nfusedist_pf%s_kpp16.snap.out',pstr{pc});
	[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
	tmplysopos = [];
	tmpphpos = [];
	tmpnfuse = [];
    tmpprot = [];
	for tc = 1:length(grouplist) % cycle over iterations
		pos = vertcat(grouplist(tc).particles(:).pos);
		ptypes = vertcat(grouplist(tc).particles(:).type);
		nfuse = vertcat(grouplist(tc).particles(:).nfuse);
        prots = vertcat(grouplist(tc).particles(:).prot);
		tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
		tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
		tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
        tmpprot = cat(1,tmpprot,prots(ptypes==2,:));
	end
	lysopos{pc} = tmplysopos;
	phpos{pc} = tmpphpos;
	nfusedist{pc} = tmpnfuse; % number of fusions for each phagosome in all trials
    phprots{pc} = tmpprot; % save proteins in all phagosomes
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
kpy = 17.9*10;
% phago production rate
kpp = 1.6*10;
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
% acidification
c0 = 1;
%afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
% number of fusion events:
Nfunc_th = @(beta,x) kpy/kpp*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,x) c0*Nfunc_th(beta,x)./(Nfunc_th(beta,x) + rp^3/ry^3);
afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);

% acidification with no change in volume
%afunc_th = @(gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
%afunc_sim = @(n,c0,rp,rl) n*c0*rl^3./(rp^3+rl^3);

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
%% plot organelle density
for pc = 1%:length(pvals)
	
	%if(length(pvals)>1); subplot(2,3,pc); end
	
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
	title(sprintf('fusion probability: %0.3f,kpp=16',pvals(pc)))
	xlim([0,0.2]);
	%ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end
%% plot acidification
% figure
%cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
cmap = jet(length(pvals));
for pc = 1%:length(pvals)
	
	%if(length(pvals)>1); subplot(2,3,pc+3); end
	
	[nphago,~,binind] = histcounts(phpos{pc},xp);
	nfusehist = accumarray(binind,nfusedist{pc})./nphago';

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xpvals,afunc_sim(nfusehist,c0,rp,ry),'o','Color',cmap(pc,:),...
        'DisplayName','acidification, sims')
	hold on
	plot(xpvals,afunc_th(betvals(pc),xpvals),'--','Color',cmap(pc,:),...
        'DisplayName','acidification, an. model');
	hold off
	%plot_cleanup
	title(sprintf('fusion probability: %0.3f, kpp=16',pvals(pc)))
	xlim([0,0.2]);
% 	ylim([0,10]);
	xlabel('distance')
	ylabel('acidity')	
end

%% plot analytical results for varying beta

% scaled domain length
Lhat = 1;

% beta values
nb = 10;
betvals = logspace(-2,2,nb);
%betvals = linspace(0.1,100,nb);

kpykpp = 11;

% phago radius
rp = 0.4; % 400nm
% lyso radius
ry = 0.1; % 100nm
% ratio of volumes
Vrat = rp^3/ry^3;

% acidification
c0 = 1;
%afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
% number of fusion events:
Nfunc_th = @(beta,x) kpykpp*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,x) c0*Nfunc_th(beta,x)./(Nfunc_th(beta,x) + rp^3/ry^3);

xvals = linspace(0,Lhat,100);

cmat = jet(nb)*0.95;
for bc = 1:nb
    plot(xvals,afunc_th(betvals(bc),xvals),'Color',cmat(bc,:),'LineWidth',2)
    hold all
    legstrings{bc} = sprintf('$\\beta = %0.2f$',betvals(bc));
end
hold off
set(gca,'FontSize',14,'defaultTextInterpreter','latex','TickLabelInterpreter','latex')
xlabel('position along axon')
legend(legstrings,'Interpreter','latex')
ylabel('acidification')