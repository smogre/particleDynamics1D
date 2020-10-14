%% degraded protein flux as function of n0
% n0 = kpy*kf*ell/vp/vy
clear;clc;
%% effect of pf on flux of degraded protein for fixed n0
pvals = logspace(-3,0,1e2);

kdmax = 1e0;
KM = 1e-1;
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
% encounter length
ell = 2*(rp+ry);
% phago production rate
kpp = 1e-1;
% beta
betavals = kpp*ell/vp/vy*kfvals;
% n0
n0 = 1;
% lyso production
kpyvals = n0.*vp.*vy./kfvals./ell./Lhat;

% acidification
c0 = 1;
% number of fusion events:
Nfunc_th = @(beta,kpy,x) kpy./kpp.*(exp(-beta.*x) - exp(-beta.*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpy,x) c0*Nfunc_th(beta,kpy,x)./(Nfunc_th(beta,kpy,x) + rp^3/ry^3);

degflux = zeros(size(pvals));
for pc = 1:length(pvals)
	A = afunc_th(betavals(pc),kpyvals(pc),xint);
	P = kdmax*A./(A+KM)/vp;
	H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
	% this is int_0^x - int_0^1 = - int_x^1
	Hx1 = H - H(end);
	G = 1 - exp(Hx1);
	
	degflux(pc) = G(1)*vp;
end

[~,idx] = max(degflux);

semilogx(pvals,degflux)
hold on
plot(pvals(idx),degflux(idx),'ko')
hold off
plot_cleanup
xlabel('fusion probability $p_f$')
ylabel('flux of degraded protein at soma')
ylim([0,1.1*degflux(idx)])
title(sprintf('$n_0 \\equiv k_{p,y}k_f\\ell L/v_pv_y = %0.3g, k_{p,p} = %0.4g$',n0,kpp))
% text(0.1,0.9,sprintf('$p_\\mathrm{max} = %0.2f$',pvals(idx)),'Units','normalized',...
% 																															'FontSize',20,'Interpreter','latex')
% text(pvals(idx)/2,1.05*degflux(idx),sprintf('$p_\\mathrm{max} = %0.4g$',pvals(idx)),'Units','data',...
% 	'FontSize',20,'Interpreter','latex')
%% sims for n0 vs pvals
nvals = 1e2;
pvals = logspace(-6,0,nvals+1);
% pvals = linspace(1/(nvals+1),0.99,nvals+1);

kdmax = 1e0;
KM = 1e-3;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% n0
betavals = logspace(0,1,nvals);

degflux = zeros(length(pvals),nvals);

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
% encounter length
ell = 2*(rp+ry);
% beta
betavals = kpp*ell/vp/vy*kfvals;
% acidification
c0 = 1;
% phago production rate
kpp = 1e-2;
% number of fusion events:
Nfunc_th = @(beta,kpy,x) kpy./kpp.*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpy,x) c0*Nfunc_th(beta,kpy,x)./(Nfunc_th(beta,kpy,x) + rp^3/ry^3);

for nc = 1:nvals
	
	disp(nc)
	
	% lyso production rate
	kpyvals = betavals(nc)*vp*vy./ell./kfvals;
	
	curbetaval = betavals(nc);
	
	parfor pc = 1:length(pvals)

		A = afunc_th(curbetaval,kpyvals(pc),xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(pc,nc) = G(1)*vp;
	end
	
end
% %% get surface plot for n0 vs pvals
pcolor(pvals,betavals,degflux')
shading flat
colormap jet;
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = 'degraded protein flux $J_d$';
cb.Label.Rotation = 270;
cb.Label.Units = 'normalized';
cb.Label.Position = [6,0.5,0];
set(gca,'XScale','log','YScale','log')
% set(gca,'YScale','log')
plot_cleanup
xlabel('fusion probability')
ylabel('$n_0\equiv k_{p,y}k_f\ell/v_pv_y$')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$k_{p,p} = %0.4g$',kdmax,KM,kpp))
%% fix n0 and change kpp
nvals = 1e2;
pvals = logspace(-1,log10(0.99),nvals+1);

kdmax = 1e0;
KM = 1e-1;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% xi
n0 = 1e2;
kppvals = logspace(-2,2,nvals);

degflux = zeros(length(pvals),nvals);

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
% encounter length
ell = 2*(rp+ry);
% beta
betavals = kpp*ell/vp/vy*kfvals;
% acidification
c0 = 1;
% number of fusion events:
Nfunc_th = @(beta,kpy,kpp,x) kpy./kpp.*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpy,kpp,x) c0*Nfunc_th(beta,kpy,kpp,x)./(Nfunc_th(beta,kpy,kpp,x) + rp^3/ry^3);

for nc = 1:nvals
	
	disp(nc)
	
	% lyso production rate
	kpyvals = n0*vp*vy./ell./kfvals;
	
	% phago production rate
	curkpp = kppvals(nc);
	
	curbetaval = curkpp*ell/vp/vy*kfvals;
	
	parfor pc = 1:length(pvals)

		A = afunc_th(curbetaval(pc),kpyvals(pc),curkpp,xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(pc,nc) = G(1)*vp;
	end
	
end
% %% get surface plot for beta
pcolor(pvals,kppvals,degflux')
shading flat
colormap jet;
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = 'degraded protein flux $J_d$';
cb.Label.Rotation = 270;
cb.Label.Units = 'normalized';
cb.Label.Position = [6,0.5,0];
set(gca,'XScale','log','YScale','log')
plot_cleanup
xlabel('fusion probability')
ylabel('$k_{p,p}$')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$\\xi\\equiv k_{p,y}k_f\\ell/v_pv_y=%0.4g$',...
																																									kdmax,KM,n0))
%% flux vs ratio and beta
nvals = 1e2;
rvals = logspace(1,3,nvals+1);
% pvals = linspace(1/(nvals+1),0.99,nvals+1);

kdmax = 1e-4;
KM = 1e-2;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% beta
betavals = logspace(-2,1,nvals);

degflux = zeros(nvals+1,nvals);

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
% encounter length
ell = 2*(rp+ry);
% acidification
c0 = 1;
% phago production rate
kpp = 1e3;
% fusion rates
kfvals = betavals*vp*vy./kpp./ell;
% number of fusion events:
Nfunc_th = @(beta,kpy,x) kpy./kpp.*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpy,x) c0*Nfunc_th(beta,kpy,x)./(Nfunc_th(beta,kpy,x) + rp^3/ry^3);
%%
for nc = 1:nvals
	
	disp(nc)
	
	
	curbetaval = betavals(nc);
	
	for pc = 1:(nvals+1)
		% lyso production rate
		kpyval = rvals(pc)*kpp;
		A = afunc_th(curbetaval,kpyval,xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(pc,nc) = G(1)*vp;
	end
	
end
% %% get surface plot for r and beta
pcolor(rvals,betavals,degflux'./kdmax)
% pcolor(rvals,betavals,degflux'.*KM./kdmax)
shading flat
colormap jet;
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
% cb.Label.String = 'normalized autophagic flux $J_d/k_\mathrm{max}$';
cb.Label.Rotation = 270;
cb.Label.Units = 'normalized';
cb.Label.Position = [6,0.5,0];
set(gca,'XScale','log','YScale','log')
% caxis([0,0.5])
% set(gca,'YScale','log')
plot_cleanup
xlabel('production rate ratio ($r\equiv k_{p,y}/k_{p,p}$)')
ylabel('$\beta\equiv k_{p,p}k_f\ell/v_pv_y$')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$',kdmax,KM))
%%
KM = 0.01;
kdmax = 1;
avals = logspace(-4,0);
kvals = kdmax.*avals./(KM+avals);
semilogx(avals,kvals)
plot_cleanup
xlabel('acidification (A)')
ylabel('degradation rate ($k_d$)')
title('$k_d = k_\mathrm{max}A/(A+K_M),~k_\mathrm{max} = 1,~K_M = 0.01$')