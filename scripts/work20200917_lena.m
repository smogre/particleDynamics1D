%% plot particle density histograms using fusion distributions
clear;clc;
pvals = [0.25,0.5];
pstr = {'0pt25','0pt5'}%,'0pt75'};

%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readlastsnap = 1;
lysopos = cell(length(pvals),1);
phpos = cell(length(pvals),1);
phprots = cell(length(pvals),1);
nfusedist = cell(length(pvals),1);
for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('../test/snapfiles/nfusedist_pf%s_kpp16.snap.out',pstr{pc});
	[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
	tmplysopos = [];
	tmpphpos = [];
	tmpnfuse = [];
    tmpprot = []; tmpacid = [];
	for tc = 1:length(grouplist) % cycle over iterations
		pos = vertcat(grouplist(tc).particles(:).pos);
		ptypes = vertcat(grouplist(tc).particles(:).type);
		nfuse = vertcat(grouplist(tc).particles(:).nfuse);
        prots = vertcat(grouplist(tc).particles(:).prot);
        acids = vertcat(grouplist(tc).particles(:).acidity);
		tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
		tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
		tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
        tmpprot = cat(1,tmpprot,prots(ptypes==2,:));
        tmpacid = cat(1,tmpacid,acids(ptypes==2,:));
	end
	lysopos{pc} = tmplysopos;
	phpos{pc} = tmpphpos;
	nfusedist{pc} = tmpnfuse; % number of fusions for each phagosome in all trials
    phprots{pc} = tmpprot; % save proteins in all phagosomes
    acidities{pc} = tmpacid;
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
%gamvals = kpy*ell/vp/vy*ry^3/(ry^3+rp^3)*kfvals;

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
phfrate = 10; % framerate for phagosomes (binwidth)
lyfrate = 1.01; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes

%% plot organelle density
for pc = 1:2%:length(pvals)
	
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
for pc = 1:2%:length(pvals)
	
	%if(length(pvals)>1); subplot(2,3,pc+3); end
	
	[nphago,~,binind] = histcounts(phpos{pc},xp);
	nfusehist = accumarray(binind,nfusedist{pc})./nphago';
    acidhist = accumarray(binind,acidities{pc})./nphago';

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xpvals,afunc_sim(nfusehist,c0,rp,ry),'o','Color',cmap(pc,:),...
        'DisplayName','acidification, sims')
    hold on
    plot(xpvals,acidhist,'*','Color',cmap(pc,:),...
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
    hold all
end
hold off

%% analytical solution for degraded proteins
pc = 1;
kdmax = 0.2;
KM = 0.08;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);
A = afunc_th(betvals(pc),xint);
plot(xint,A)
P = kdmax*A./(A+KM)/vp;
H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
% this is int_0^x - int_0^1 = - int_x^1
Hx1 = H - H(end);
G = 1 - exp(Hx1);

plot(xint,G)
%% plot degraded protein
% figure
%cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
cmap = jet(length(pvals));
for pc = 1:2%:length(pvals)
    
    % prediction for degraded prots
    A = afunc_th(betvals(pc),xint);
    P = kdmax*A./(A+KM)/vp;
    H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
    % this is int_0^x - int_0^1 = - int_x^1
    Hx1 = H - H(end);
    G = 1 - exp(Hx1);
    
    
	%if(length(pvals)>1); subplot(2,3,pc+3); end
	
    % where are all the phagosomes located
	[nphago,~,binind] = histcounts(phpos{pc},xp);
    % avg degraded prot per phagosome in each bin
	prothist = accumarray(binind,phprots{pc}(:,2))./nphago';

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
	dxp = xpvals(2)-xpvals(1);

	plot(xpvals,prothist,'o','Color',cmap(pc,:),...
        'DisplayName','acidification, sims')
	hold on
    plot(xint,G,'--','Color',cmap(pc,:))
	%plot(xpvals,afunc_th(betvals(pc),xpvals),'--','Color',cmap(pc,:),...
    %    'DisplayName','acidification, an. model');
	hold off
	%plot_cleanup
	title(sprintf('fusion probability: %0.3f, kpp=16',pvals(pc)))
	%xlim([0,0.2]);
% 	ylim([0,10]);
	xlabel('distance')
	ylabel('degraded protein')	
    hold all
end
hold off

%% effect of pf on flux of degraded protein
pvals = logspace(-3,0);

kdmax = 0.2;
KM = 0.8;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

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

degflux = zeros(size(pvals));
for pc = 1:length(pvals)
    A = afunc_th(betvals(pc),xint);
    P = kdmax*A./(A+KM)/vp;
    H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
    % this is int_0^x - int_0^1 = - int_x^1
    Hx1 = H - H(end);
    G = 1 - exp(Hx1);
    
    degflux(pc) = G(1)*vp;
end

%
semilogx(pvals,degflux)