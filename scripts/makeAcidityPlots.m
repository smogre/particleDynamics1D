%% make plots to explore half-saturation point and maximum acidification
% the half saturation point is the x value at which the acidity of the
% autophagosome is 50% of that at the soma
% using analytical model
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
ry = 0.2/L; % 100nm/200nm
% encounter time
T = 2*(rp+ry)/(vp+vy);
% encounter length
ell = 2*(rp+ry);
% max rate for protein degradation
kmax_p = 1;
% saturation parameter for protein degradation
KM_p = 1;
%% set other params
% phagosome production rate
kpp = 1.6e0;
% lysosome production rate
kpy = 17.9e0;
% fusion probability
pvals = [0.25,0.5,0.75];
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betvals = kpp*ell/vp/vy*kfvals;
% production ratio
prodratio = kpy/kpp;
% x values
xvals = linspace(0,1,5e3); dx = xvals(2)-xvals(1);
% options
options = struct();
% phagosome production rate
options.kpp = kpp;
% use dilution model to calculate acidification?
options.dilute = false;
% maximum acidification rate
options.kmax_a = 1;
% acidification saturation parameter
options.KM_a = 1;
% lysosome radius
options.ry = ry;
% string to use for titles
dilutestr = "vATPase model";
if(options.dilute); dilutestr = "dilution model"; end

%% plot lysosome density and acidification vs distance
clf

cmap = colormapinterp([0,0,1;1,0,0],length(betvals));
for bc = 1:length(betvals)
	
	[~,lysodist,~,~,avals] = getProtFlux(betvals(bc),prodratio,kmax_p,KM_p,xvals,options);
	
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
end
hold off
if(options.dilute)
	title(sprintf('solid lines: density, dashed lines: acidification\n %s, $r_y: %0.3g\\mu$m',...
		dilutestr,ry*L))
else
	title(sprintf('solid lines: density, dashed lines: acidification\n %s, $k_\\mathrm{max}^a: %0.3g, K_M^a:%0.3g, r_y: %0.3g\\mu$m',...
		dilutestr,options.kmax_a,options.KM_a,ry*L))
end
%% plot acidification vs lysosome density
clf
cmap = colormapinterp([0,0,1;1,0,0],length(betvals));
for bc = 1:length(betvals)
	[~,lysodist,~,~,avals] = getProtFlux(betvals(bc),prodratio,kmax_p,KM_p,xvals,options);
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
if(options.dilute)
	title(sprintf('%s, $r_y: %0.3g\\mu$m',dilutestr,ry*L))
else
	title(sprintf('%s, $k_\\mathrm{max}^a: %0.3g, K_M^a:%0.3g, r_y: %0.3g\\mu$m',...
		dilutestr,options.kmax_a,options.KM_a,ry*L))
end
hold off
%% calculations for half-saturation length vs rate ratio and beta
% for a fixed KM value
nvals = 50;

KM_a = 10;

ratiovals = logspace(-1,2,nvals);
betvals = logspace(-2,1,nvals+1);

avals = cell(nvals,nvals+1);

halfsat = zeros(nvals,nvals+1);
amax = zeros(nvals,nvals+1);

options.KM_a = KM_a;

for rc = 1:nvals
	disp(rc);
	options.kmax_a = ratiovals(rc).*options.KM_a;
	
	for bc = 1:(nvals+1)
		betval = betvals(bc);
		
		[~,~,~,~,avals{rc,bc},~,opt] = getProtFlux(betval,prodratio,kmax_p,KM_p,xvals,options);
		halfsat(rc,bc) = xvals(find(avals{rc,bc}<avals{rc,bc}(1)/2,1));
		amax(rc,bc) = max(avals{rc,bc});
	end
end

%% pcolor plots for half saturation/ max acidity vs kmax and beta
% pcolor(ratiovals,betvals,halfsat')
pcolor(ratiovals,betvals,amax')
shading flat
set(gca,'xScale','log','yScale','log')
cb = colorbar;
% cb.Label.String = 'scaled half-saturation length';
cb.Label.String = 'maximum acidification';
cb.Label.Interpreter = 'latex';
cb.FontSize = 30;
cb.TickLabelInterpreter = 'latex';
plot_cleanup('FontSize',30)
xlabel('scaled maximum acidification rate $k_\mathrm{max}^a/K_M^a$')
ylabel({'lysosome density decay lengthscale','$\beta\equiv k_{p,p}k_f\ell/v_pv_y$'})

title(sprintf('%s, $K_M^a:%0.3g$',dilutestr,options.KM_a))
%% slices from pcolor plots
rc = 11;
for rc = 26:5:51
	semilogx(ratiovals,halfsat(:,rc),'DisplayName',sprintf('KM: %0.2g',KMvals(rc)))
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
%% plot half saturation/max acidity for varying pf values
pvals = logspace(-2,0,1e3);
% fusion rates
kfvals = -log(1-pvals)/T;
% beta
betvals = kpp*ell/vp/vy*kfvals;

halfsat = zeros(length(betvals),1);
amax = zeros(length(betvals),1);

if(options.dilute)
ryvals = [0.1,0.2]/L;
	for rc = 1:length(ryvals)
		options.ry = ryvals(rc);
		for bc = 1:length(betvals)
			disp(bc)
			[~,lysodist,~,nfusedist,avals,~] = getProtFlux(betvals(bc),prodratio,kmax_p,KM_p,xvals,options);
			hsloc = find(avals<avals(1)/2,1);
			if(~isempty(hsloc))
				halfsat(bc) = xvals(hsloc);
				amax(bc) = max(avals);
			else
				halfsat(bc) = nan;
				amax(bc) = nan;
			end
		end
		
% 			semilogx(pvals,halfsat)
		semilogx(pvals,amax)
		
% 			[hmax,idx] = max(halfsat);
% 		text(pvals(idx)/2,1.02*halfsat(idx),sprintf('$r_y = %0.2g\\mu$m',options.ry*L),'Units','data',...
% 			'FontSize',30,'Interpreter','latex')
		
		[vmax,idx] = max(amax);
		text(pvals(idx)/2,1.02*amax(idx),sprintf('$r_y = %0.2g\\mu$m',options.ry*L),'Units','data',...
			'FontSize',30,'Interpreter','latex')
		
		hold on
	end
	hold off
	plot_cleanup('FontSize',30,'LineWidth',4)
	xlabel('fusion probability')
% 	ylabel('half-saturation length')
	ylabel('max acidification')
	
	title(dilutestr)
else
	options.KM_a = 1;
	kmaxvals = [1,5,10];
	
	for kc = 1:length(kmaxvals)
		options.kmax_a = kmaxvals(kc);
		for bc = 1:length(betvals)
			disp(bc)
			[~,lysodist,~,nfusedist,avals,~] = getProtFlux(betvals(bc),prodratio,kmax_p,KM_p,xvals,options);
			hsloc = find(avals<avals(1)/2,1);
			if(~isempty(hsloc))
				halfsat(bc) = xvals(hsloc);
				amax(bc) = max(avals);
			else
				halfsat(bc) = nan;
				amax(bc) = nan;
			end
		end
		
% 		semilogx(pvals,halfsat)
		semilogx(pvals,amax)
		
% 			[hmax,idx] = max(halfsat);
% 			text(pvals(idx)/2,1.02*halfsat(idx),sprintf('$k_\\mathrm{max}^{a} = %0.2g$/s',options.kmax_a),'Units','data',...
% 			'FontSize',30,'Interpreter','latex')
		
		[vmax,idx] = max(amax);
		text(pvals(idx)/2,1.02*amax(idx),sprintf('$k_\\mathrm{max}^a = %0.2g$/s',options.kmax_a),'Units','data',...
			'FontSize',30,'Interpreter','latex')
		
		hold on
	end
	hold off
	plot_cleanup('FontSize',30,'LineWidth',4)
	xlabel('fusion probability')
% 	ylabel('half-saturation length')
	ylabel('max acidification')
	
	title(sprintf('%s, $K_M^a:%0.3g$',dilutestr,options.KM_a))
end