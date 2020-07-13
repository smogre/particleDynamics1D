%% test particle distribution
clear;clc
pvals = [0,0.25,0.5,0.75,1];%[0,1];
pstr = [{'0'},{'0pt25'},{'0pt50'},{'0pt75'},{'1'}];%[{'0'},{'1'}];
L = 100;
%%
fusedist = [];
partdist = [];
for pc = 1:length(pvals)
	disp(pc)
	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/nfusedist_pf%s.fuse.out',pstr{pc});
	% 	filename = sprintf('/data/SimData/particledynamics1D/nfusefiles/testnfuse.fuse.out',pstr{pc});
	[tvals,tmpfusedist,tmppartdist,nreg,regedges,ntrials] = readfusedist(filename);
	fusedist = cat(4,fusedist,tmpfusedist);
	partdist = cat(4,partdist,tmppartdist);
end
nsave = size(fusedist,1);
ntype = size(fusedist,3);
pdensity = partdist./ntrials.*L;
meanpdensity = squeeze(sum(partdist,2))/L/ntrials;
%% plot particle distribution
ptype = 1;
typestr = [{'lyso'},{'phago'}];
% figure
ssinds = nsave;
cmap = [0,0,1;1,0,0];
fc = 1;
xvals = regedges(1:end-1);
% for tc = 1:20%nsave
	ct = 1;
	for pc = 1:length(pvals)
% 		disp(pc)
	subplot(2,3,ct)
% 		subplot(2,2,ct)
		for ptype = 1:ntype
			
% 		meanpartdist = smooth(mean(partdist(ssinds,:,ptype,pc),1)./ntrials,11);
% 		meanpartdist = mean(partdist(ssinds,:,ptype,pc),1)./ntrials;
		pdensity_ss = mean(pdensity(ssinds,:,ptype,pc),1);
% 		pd = fitdist(meanpartdist','kernel');
% 		yvals = pdf(pd,regedges(1:end-1));
% 			meanpartdist = smooth(mean(partdist(tc,:,ptype,pc),1)./ntrials,11);
% 		plot(regedges(1:end-1),meanpartdist,'Color',cmap(ptype,:),'DisplayName',...
% 			plot(xvals,meanpartdist,'Color',cmap(ptype,:),'DisplayName',...
% 				sprintf('%s: %0.3f/$\\mu$m',typestr{ptype},mean(meanpdensity(ssinds,ptype,pc))))
			plot(xvals,pdensity_ss,'Color',cmap(ptype,:),'DisplayName',...
				sprintf('%s: %0.3f/$\\mu$m',typestr{ptype},mean(meanpdensity(ssinds,ptype,pc))))
% 			sprintf('pfuse: %0.3f, particle density: %0.3f',pvals(pc),mean(pdensity(ssinds,ptype,pc))))
% 	histogram('BinEdges',regedges,'BinCounts',mean(partdist(ssinds,:,ptype,pc),1)./ntrials)
% 		plot(regedges(1:end-1),yvals,'Color',cmap(ptype,:),'DisplayName',...
% 				sprintf('%s: %0.3f/$\\mu$m',typestr{ptype},mean(pdensity(ssinds,ptype,pc))))
			hold on
			if(ptype==1)
				plot(xvals,rho_y(betvals(pc),xvals),'--','Color',cmap(ptype,:),'DisplayName',...
																sprintf('%s density, an. model',typestr{ptype}))
			else
				plot(xvals,rho_p(xvals),'--','Color',cmap(ptype,:),'DisplayName',...
																sprintf('%s density, an. model',typestr{ptype}))
			end
			plot_cleanup
			xlabel('distance')
			xlim([0,0.2]);
			ylabel('organelle density')

			title(sprintf('fusion probability: %0.2f',pvals(pc)))
			lgd = legend;
			lgd.Interpreter = 'latex';
		end
		
% 		pause(0.2)
		hold off
		ct = ct+1;
% 	end
% 	drawnow
% 	M(fc) = getframe(gcf);
% 	fc = fc+1;
end
%%
animation2movie(M,'npartdist_smooth.avi',16)
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
