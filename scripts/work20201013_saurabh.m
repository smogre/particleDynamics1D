%% flux vs ratio and beta (vary kpp)
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
% fusion rates
kfvals = 1e4;
% phago production rate
kppvals = betavals*vp*vy./kfvals./ell;
% number of fusion events:
Nfunc_th = @(beta,kpy,kpp,x) kpy./kpp.*(exp(-beta*x) - exp(-beta*Lhat));
% acidification with increasing volume
afunc_th = @(beta,kpy,kpp,x) c0*Nfunc_th(beta,kpy,kpp,x)./(Nfunc_th(beta,kpy,kpp,x) + rp^3/ry^3);

for nc = 1:nvals
	
	disp(nc)
	
	
	curbetaval = betavals(nc);
	curkppval = kppvals(nc);
	
	for pc = 1:(nvals+1)
		% lyso production rate
		kpyval = rvals(pc)*curkppval;
		A = afunc_th(curbetaval,kpyval,curkppval,xint);
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