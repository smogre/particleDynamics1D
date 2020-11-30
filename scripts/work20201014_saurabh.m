%% make plots for organelle distribution, acidity, and degraded protein flux
% use analytical model
clear; clc;
%% set default values
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
%% set remaining params
% phagosome production
kpp = 1.6e0;
% lysosome production
kpy = 17.9e0;
% fusion probability
% pvals = [0.25,0.5,0.75];
pvals = logspace(-2,0,1e3);
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betvals = kpp*ell/vp/vy*kfvals;
% production ratio
prodratio = kpy/kpp;
% max degradation rate
kmax = 1;
% saturating acidification
KM = 1e0;
% x values
xvals = linspace(0,1,5e3); dx = xvals(2)-xvals(1);
% options
options = struct();
options.kpp = kpp;
options.dilute = false;


%% plot lysosome density and acidification vs distance
clf
cmap = colormapinterp([0,0,1;1,0,0],length(betvals));
for bc = 1:length(betvals)
% 	subplot(1,3,bc)
	[~,lysodist,~,~,avals] = getProtFlux(betvals(bc),prodratio,kmax,KM,xvals,options);
	
	yyaxis left
	hold on
	plot(xvals,lysodist./max(lysodist),'-','Color',cmap(bc,:),'DisplayName',...
																										 sprintf('fusion probability: %0.2g',pvals(bc)))
	plot_cleanup
	if(bc==length(betvals)); hold off; end
	ylabel('lysosome density relative to soma')
	
	hold on
	yyaxis right
	plot(xvals,avals,'--','Color',cmap(bc,:))
	plot_cleanup
	if(bc==length(betvals)); hold off; end
	ylabel('relative phagosome acidification')
	xlabel('distance from soma as fraction of axon length')
% 	xlim([0,0.1])
end
hold off
title('solid lines: density, dashed lines: acidification')
%% plot acidification vs lysosome density
clf
cmap = colormapinterp([0,0,1;1,0,0],length(betvals));
for bc = 1:length(betvals)
	[~,lysodist,~,~,avals] = getProtFlux(betvals(bc),prodratio,kmax,KM,xvals,options);
	plot(lysodist./max(lysodist),avals,'Color',cmap(bc,:))
	hold on
	plot_cleanup
	xlabel('lysosome density relative to soma')
	ylabel('relative phagosome acidification')
	xlim([min(lysodist),max(lysodist)])
	xlim([0,1])
end
hold off
%% plot optimum fusion probability
hold off
% hold on
degflux = zeros(length(betvals),1);
for bc = 1:length(betvals)
	degflux(bc) = getProtFlux(betvals(bc),prodratio,kmax,KM,xvals,options);
end
% %%
semilogx(pvals,degflux)
plot_cleanup('FontSize',30)
% ylim([0,0.27])
xlabel('fusion probability')
xticks(logspace(-5,0,6));
ylabel('autophagic flux per phagosome')
[jmax,idx] = max(degflux);
text(pvals(idx)/2,1.02*degflux(idx),sprintf('$k_\\mathrm{max} = %0.4g$/s',kmax),'Units','data',...
	'FontSize',30,'Interpreter','latex')
%% surface plots for autophagic flux vs beta and prodratio
betvals = logspace(-2,1,1e2);
prodratiovals = logspace(1,3,101);
degflux = zeros(length(betvals),length(prodratiovals));

for bc = 1:length(betvals)
	disp(bc)
	for rc = 1:length(prodratiovals)
		degflux(bc,rc) = getProtFlux(betvals(bc),prodratiovals(rc),kmax,KM,xvals,options);
	end
end
%%
pcolor(prodratiovals,betvals,degflux)
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
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$',kmax,KM))