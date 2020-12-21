%% make plots to explore half-saturation point
% the half saturation point is the x value at which the acidity of the
% autophagosome is 50% of that at the soma
% acidification is dependent on vATPase pump activity
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
% vy = 0.4*2.67;	% 2um/s
% phago radius
rp = 0.4/L; % 400nm
% lyso radius
ry = 0.2/L; % 100nm
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
pvals = [0.25,0.5,0.75];
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
options.kmax_a = 1;
options.KM_a = 10;

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
	plot_cleanup('FontSize',30,'LineWidth',4)
	if(bc==length(betvals)); hold off; end
	ylabel('lysosome density relative to soma')
	
	hold on
	yyaxis right
	plot(xvals,avals,'--','Color',cmap(bc,:))
	plot_cleanup('FontSize',30,'LineWidth',4)
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
	plot(lysodist./max(lysodist),avals,'Color',cmap(bc,:),'DisplayName',...
																												sprintf('fusion probability: %g',pvals(bc)))
	hold on
	plot_cleanup('FontSize',30,'LineWidth',4)
	xlabel('lysosome density relative to soma')
	ylabel('relative phagosome acidification')
	xlim([min(lysodist),max(lysodist)])
	xlim([0,1])
end
lgd = legend;
lgd.Interpreter = 'latex';
lgd.Location = 'northwest';
hold off
%% plot optimum fusion probability
pvals = logspace(-2,0,1e3);
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betvals = kpp*ell/vp/vy*kfvals;
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
%% surface plots for 50% acidification vs KM_a and kmax_a
pvals = 0.5;
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betval = kpp*ell/vp/vy*kfvals;

ratiovals = logspace(-1,2,50);
KMvals = logspace(-1,2,51);
% kmaxvals = logspace(-1,2,50);
% kmaxvals = ratiovals.*KMVals;
avals = cell(length(ratiovals),length(KMvals));

halfsat = zeros(length(ratiovals),length(KMvals));
amax = zeros(length(ratiovals),length(KMvals));

options = struct();
options.dilute = false;
for kc2 = 1:length(KMvals)
	disp(kc2);
		options.KM_a = KMvals(kc2);
% for kc1 = 1:length(ratiovals)
	for kc1 = 1:length(ratiovals)
	
% 	options.kmax_a = kmaxvals(kc1);
	options.kmax_a = ratiovals(kc1).*options.KM_a;
	
		[~,~,~,~,avals{kc1,kc2},~] = getProtFlux(betval,prodratio,kmax,KM,xvals,options);
		halfsat(kc1,kc2) = xvals(find(avals{kc1,kc2}<avals{kc1,kc2}(1)/2,1));
		amax(kc1,kc2) = max(avals{kc1,kc2});
	end
end

%% pcolor plots for half saturation/ max acidity
% pcolor(kmaxvals./KMvals',KMvals,halfsat')
pcolor(ratiovals,KMvals,halfsat')
% pcolor(kmaxvals,KMvals,amax')
shading flat
set(gca,'xScale','log','yScale','log')
cb = colorbar;
% cb.Label.String = 'maximum acidification';
cb.Label.String = 'scaled half-saturation length';
cb.Label.Interpreter = 'latex';
cb.FontSize = 30;
cb.TickLabelInterpreter = 'latex';
plot_cleanup('FontSize',30)
xlabel('scaled maximum acidification rate $k_\mathrm{max}^a/K_M^a$')
ylabel('acidification saturation parameter $K_M^a$')
title(sprintf('fusion probability: %g',pvals));
%%
kc2 = 11;
for kc2 = 26:5:51
semilogx(ratiovals,halfsat(:,kc2),'DisplayName',sprintf('KM: %0.2g',KMvals(kc2)))
% semilogx(kmaxvals./KMvals(kc2),amax(:,kc2),'DisplayName',sprintf('KM: %g',KMvals(kc2)))
hold on
end
hold off
plot_cleanup
% title(sprintf('KM: %g, pf: %g',KMvals(kc2),pvals));
title(sprintf('pf: %g',pvals));
xlabel('kmax/KM')
ylabel('half-saturation point')
% ylabel('maximum acidification')
% xlim([min(kmaxvals),max(kmaxvals)]);
%%
kc1 = 20;
kc1inds = 0:10:50;
kc1inds(1) = 1;
for kc1 = kc1inds
semilogx(KMvals,halfsat(kc1,:),'DisplayName',sprintf('kmax: %g',kmaxvals(kc1)));
hold on
end
hold off
plot_cleanup
xlabel('KM')
% ylabel('half-saturation point')
ylabel('maximum acidification')
% title(sprintf('kmax: %g',kmaxvals(kc1)));
xlim([min(KMvals),max(KMvals)]);
%% calculate half saturation/max acidity for varying pf values
pvals = 1e-6;logspace(-6,0,1e3);
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betvals = kpp*ell/vp/vy*kfvals;

halfsat = zeros(length(betvals),1);
amax = zeros(length(betvals),1);

options = struct();
options.dilute = false;
options.kpp = kpp;
options.kmax_a = 10;
options.KM_a = 10;
for bc = 1:length(betvals)
	disp(bc)
	[~,lysodist,~,nfusedist,avals,~] = getProtFlux(betvals(bc),prodratio,kmax,KM,xvals,options);
	hsloc = find(avals<avals(1)/2,1);
	if(~isempty(hsloc))
		halfsat(bc) = xvals(hsloc);
		amax(bc) = max(avals);
	else
		halfsat(bc) = nan;
		amax(bc) = nan;
	end
end
%%
% hold on
semilogx(pvals,halfsat)
plot_cleanup('FontSize',30,'LineWidth',4)
% ylim([0,0.27])
xlabel('fusion probability')
% xticks(logspace(-5,0,6));
ylabel('half-saturation length')
% ylabel('max acidification')
% [hmax,idx] = max(halfsat);
[vmax,idx] = max(amax);
% text(pvals(idx)/2,1.02*halfsat(idx),sprintf('$p_\\mathrm{max} = %0.2g$',pvals(idx)),'Units','data',...
% 'FontSize',30,'Interpreter','latex')
text(pvals(idx)/2,1.02*halfsat(idx),sprintf('$k_\\mathrm{max}^{a} = %0.2g$/s',options.kmax_a),'Units','data',...
	'FontSize',30,'Interpreter','latex')
