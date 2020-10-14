%% degraded protein flux as function of beta
clear;clc;
%% effect of pf on flux of degraded protein for fixed beta
pvals = logspace(-3,log10(0.99),1e1);

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
prodmult = 1e0;
% encounter length
ell = 2*(rp+ry);
% set beta value
% betval = kpp*ell/vp/vy*kfvals;
betval = 0.001;
% lyso production rate
kpy = 100;
% kpy = 100;
% phago production rate
% kpp = 10;
kppvals = betval*vp*vy./ell./kfvals;

% acidification
c0 = 1;
% number of fusion events:
Nfunc_th = @(beta,kpp,x) kpy./kpp.*(exp(-beta.*x) - exp(-beta.*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpp,x) c0*Nfunc_th(beta,kpp,x)./(Nfunc_th(beta,kpp,x) + rp^3/ry^3);

degflux = zeros(size(pvals));
for pc = 1:length(pvals)
	A = afunc_th(betval,kppvals(pc),xint);
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
title(sprintf('$\\beta = %0.3g$',betval))
% text(0.1,0.9,sprintf('$p_\\mathrm{max} = %0.2f$',pvals(idx)),'Units','normalized',...
% 																															'FontSize',20,'Interpreter','latex')
text(pvals(idx)/2,1.05*degflux(idx),sprintf('$p_\\mathrm{max} = %0.4g$',pvals(idx)),'Units','data',...
	'FontSize',20,'Interpreter','latex')
%% sims for beta vs pvals
nvals = 1e2;
pvals = logspace(-1,0,nvals+1);

kdmax = 1e0;
KM = 1e-1;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% beta
betvals = logspace(-3,0,nvals);

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
% acidification
c0 = 1;
% lyso production rate
kpy = 1;
% number of fusion events:
Nfunc_th = @(beta,kpp,x) kpy./kpp.*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpp,x) c0*Nfunc_th(beta,kpp,x)./(Nfunc_th(beta,kpp,x) + rp^3/ry^3);

for mc = 1:nvals
	
	disp(mc)
	
	% phago production rate
	kppvals = betvals(mc)*vp*vy./ell./kfvals;
	
	curbeta = betvals(mc);
	
	parfor pc = 1:length(pvals)

		A = afunc_th(curbeta,kppvals(pc),xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(pc,mc) = G(1)*vp;
	end

end
%% get surface plot for beta
pcolor(pvals,betvals,degflux')
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
ylabel('$\beta\equiv k_{p,p}\ell k_f/v_pv_y$')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$k_{p,y} = %0.4g$',kdmax,KM,kpy))