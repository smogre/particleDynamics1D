%% plot particle distribution movie
clear;clc
pvals = 0;
pstr = {'0pt50'};%[{'0'},{'1'}];
L = 100;
fname = './nfusedist_pf0pt50.snap.out';
[tvals,fusedist,partdist,nreg,regedges,ntrials] = readfusedist(fname);
dx = regedges(2)-regedges(1);
dt = 1e-3;
%%
nsave = size(fusedist,1);
ntype = size(fusedist,3);
pdensity = partdist./ntrials./dx;
meanpdensity = squeeze(mean(pdensity,2));
%% plot particle distribution
typestr = [{'lyso'},{'phago'}];
% figure
cmap = [0,0,1;1,0,0];
fc = 1;
xvals = regedges(1:end-1);
for tc = 1:1:nsave
	ct = 1;
	for pc = 1:length(pvals)
		
		for ptype = 1:ntype
			
			plot(xvals,pdensity(tc,:,ptype),'Color',cmap(ptype,:),'DisplayName',...
				sprintf('%s: %0.3f/L',typestr{ptype},mean(meanpdensity(tc,ptype,pc))))
			hold on
			plot(xvals,1e3*ones(size(xvals)),'.k')
			plot_cleanup
			xlabel('distance')
			xlim([0,1]);
			ylim([0,2e3]);
			ylabel('organelle density')

		end
		
		% 		pause(0.2)
		hold off
		ct = ct+1;
	end
	title(sprintf('fusion probability: %0.3f, time: %0.3f',pvals(pc),tc*dt))
	drawnow
% 		M(fc) = getframe(gcf);
% 		fc = fc+1;
end
%%
animation2movie(M,'partdensity_pf0pt25.avi',16)
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
ell = (rp+ry);
% lyso acidification
c0 = 1;
% beta
betvals = kpp*ell/vp/vy*kfvals;
% gamma
gamvals = kpy*ell/vp/vy*vy/(vy+vp)*kfvals;
% acidicication
afunc_th = @(c0,gam,beta,x) c0*(1-exp(gam/beta*(exp(-beta*Lhat)-exp(-beta*x))));
% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_y0 = kpy/vy;
rho_y = @(beta,x) rho_y0*exp(-beta.*x);
%%
hold on
for nc = 1:10:length(regedges(1:end-1))
	plot(regedges(nc)*[1,1],ylim,'k--')
end
hold off
%%
meanfusedist(isnan(meanfusedist)) = 0;
cvals = cfunc(meanfusedist,c0,rp,rl);
