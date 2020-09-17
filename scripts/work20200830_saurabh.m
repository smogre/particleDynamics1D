%% test numerical ode solver
pvals = 0.5;[1e-4,0.1,0.5,0.75];
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
% amount of degraded protein near cell body
k0 = 1;
kM = 0.01;
degprot_th = @(c0,k0,Km,gam,beta) (1-exp(-(k0/vp)*afunc_th(c0,gam,beta,0)/(afunc_th(c0,gam,beta,0)+Km)));

%%
plotinds = 1:length(pvals);
cmap = colormapinterp([0,0,1;1,0,0],length(plotinds));
G0vals = zeros(length(pvals),1);
fluxvals = zeros(length(pvals),1);
ct = 1;
for pc = plotinds
% 	[x,G] = ode45(@(x,G) degProtdist(x,G,k0,kM,gamvals(pc),betvals(pc),vp),linspace(1,0,1e3),0);
	[x,G] = ode45(@(x,G) degProtdist(x,G,k0,kM,gamvals(pc),betvals(pc),vp),linspace(1,0,1e3),0);
	G0vals(pc) = G(end);
	fluxvals(pc) = -vp*(G(end-1)-G(end))./(x(end-1)-x(end));
% 	plot(x,G,'Color',cmap(ct,:),'DisplayName',sprintf('pf: %0.2g',pvals(pc)))
% 	hold on
% 	plot(x,1-exp(k0.*(x-1)),'--','Color',cmap(ct,:))
	ct = ct+1;
end
% hold off
% plot_cleanup
% lgd = legend;
% lgd.Interpreter = 'latex';
% xlim([0,0.1])
% tmpy = ylim;
% ylim([0,tmpy(2)])
%%
plot(x,G)
% xlim([0,0.1])
plot_cleanup
%%
semilogx(pvals,G0vals)
% semilogx(pvals,fluxvals)
plot_cleanup

%%
function dy = degProtdist(x,y,k0,kM,gamvals,betvals,vp)
	dy = -degrate(k0,kM,acidification(gamvals,betvals,x),vp)*(1-y);
end
%%
function kd = degrate(k0,kM,aval,vp)
% kd = k0.*aval./(aval+kM)./vp;
% kd = k0;
	kd = aval;
end
%%
function aval = acidification(gam,bet,x)
% 	aval = 1-exp(gam/bet.*(exp(-bet)-exp(-bet.*x)));
	aval = 1-x;
end
