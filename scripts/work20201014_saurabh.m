%% make plots for organelle distribution, acidity, and degraded protein flux
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
kmax = 5;
% saturating acidification
KM = 1e0;
% x values
xvals = linspace(0,1,5e3); dx = xvals(2)-xvals(1);
% options
options = struct();
options.kpp = kpp;


%%
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
%%
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
%%
% hold off
hold on
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