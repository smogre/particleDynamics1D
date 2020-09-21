%% explore parameter space for degrated protein content and flux
clear;clc;
%% effect of pf on flux of degraded protein
pvals = logspace(-3,0,1e3);

kdmax = 1e0;
KM = 1e0;
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
% production rate multiplier
prodmult = 1e3;
% lyso production rate
kpy = 17.9*prodmult;
% phago production rate
kpp = 1.6*prodmult;
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

[~,idx] = max(degflux);

semilogx(pvals,degflux)
hold on
plot(pvals(idx),degflux(idx),'ko')
hold off
plot_cleanup
xlabel('fusion probability $p_f$')
ylabel('flux of degraded protein at soma')
ylim([0,1.1*degflux(idx)])
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$ k_{p,p} = %0.4g,~k_{p,y} = %0.4g$',...
	kdmax,KM,kpp,kpy))
% text(0.1,0.9,sprintf('$p_\\mathrm{max} = %0.2f$',pvals(idx)),'Units','normalized',...
% 																															'FontSize',20,'Interpreter','latex')
text(pvals(idx)/2,1.05*degflux(idx),sprintf('$p_\\mathrm{max} = %0.4g$',pvals(idx)),'Units','data',...
	'FontSize',20,'Interpreter','latex')
%%

pvals = logspace(-3,0,1e3);

kdmax = 1e0;
KM = 1e0;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% production rate multiplier
prodmult = logspace(0,3,1e3);

% ratio of prod rates kpy/kpp
prodratio = logspace(-2,3,1e3);

degflux = zeros(length(pvals),length(prodmult));
pmax =  zeros(length(prodmult),1);
Jmax = pmax;

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
% acidification
c0 = 1;

% for mc = 1:length(prodmult)
for mc = 1:length(prodratio)
	
	disp(mc)
	% phago production rate
% 	kpp = 1.6*prodmult(mc);
	kpp = 10;
	% lyso production rate
% 	kpy = 17.9*prodmult(mc);
	kpy = prodratio(mc)*kpp;
	
	% beta
	betvals = kpp*ell/vp/vy*kfvals;
	% density
	rho_p = @(x) kpp/vp.*ones(size(x));
	rho_y0 = kpy/vy;
	rho_y = @(beta,x) rho_y0*exp(-beta.*x);
	
	% number of fusion events:
	Nfunc_th = @(beta,x) kpy/kpp*(exp(-beta*x) - exp(-beta*Lhat));
	% acidification with increasing volume
	afunc_th = @(beta,x) c0*Nfunc_th(beta,x)./(Nfunc_th(beta,x) + rp^3/ry^3);
	
	
	parfor pc = 1:length(pvals)
		A = afunc_th(betvals(pc),xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(pc,mc) = G(1)*vp;
	end
	
	[Jmax(mc),idx] = max(degflux(:,mc));
	pmax(mc) = pvals(idx);
end
%% plot optimum pf vs production rate multiplier
% kpp = 1.6;
% kpy = 17.9;

semilogx(prodmult,pmax)
xlabel('production rate multiplier ($n$)')
ylabel('$p_f$ with maximum flux')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$ k_{p,p} = %0.4gn,~k_{p,y} = %0.4gn$',...
	kdmax,KM,kpp,kpy))
plot_cleanup
%% get surface plot for production rate multiplier
pcolor(pvals,prodmult,degflux)
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
ylabel('production rate multiplier ($n$)')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$ k_{p,p} = %0.4gn,~k_{p,y} = %0.4gn$',...
	kdmax,KM,kpp,kpy))
%% plot optimum pf vs ratio of production rates
% kpp = 1.6;
% kpy = 17.9;

semilogx(prodratio,pmax)
xlabel('production rate ratio ($r$)')
ylabel('$p_f$ with maximum flux')
xlim([min(prodratio),max(prodratio)])
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$, $ k_{p,p} = %0.4g,~k_{p,y} = rk_{p,p}$',...
	kdmax,KM,1))
plot_cleanup
%% get surface plot for production rate ratio
pcolor(pvals,prodratio,degflux)
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
ylabel('production rate ratio ($r$)')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$, $ k_{p,p} = %0.4g,~k_{p,y} = rk_{p,p}$',...
	kdmax,KM,kpp))
