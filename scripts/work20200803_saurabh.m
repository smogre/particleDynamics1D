%% plot organelle densities with reversal
clear;clc;
pvals = [0.25,0.5,0.75];%[0,0.25,0.5,0.75,1];
pstr = [{'0pt25'},{'0pt50'},{'0pt75'}];
ellvals = [0.0158,0.0066,0.0033];%0.1;
ellstr = {'beta'};
%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readsinglesnap = 1;
lysopos = cell(length(ellstr),length(pvals));
phpos = cell(length(ellstr),length(pvals));
nfusedist = cell(length(ellstr),length(pvals));
for lc = 1:length(ellstr)
	for pc = 1:length(pvals)
		disp([lc pc])
		filename = sprintf('/data/SimData/particledynamics1D/snapfiles/nfusedist_pf%s_ell%s.snap.out',...
																																							pstr{pc},ellstr{lc});
		% 	filename = '/data/SimData/particledynamics1D/snapfiles/testnfuse.snap.out';
		[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
		tmplysopos = [];
		tmpphpos = [];
		tmpnfuse = [];
		for tc = 1:length(grouplist)
			pos = vertcat(grouplist(tc).particles(:).pos);
			ptypes = vertcat(grouplist(tc).particles(:).type);
			nfuse = vertcat(grouplist(tc).particles(:).nfuse);
			tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
			tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
			tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
		end
		lysopos{lc,pc} = tmplysopos;
		phpos{lc,pc} = tmpphpos;
		nfusedist{lc,pc} = tmpnfuse;
	end
end
toc
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
delt = 2*(rp+ry);
% switch rates
ksvals = vy./ellvals;
% lam1
lam1 = @(ks,kf) ks./vy+kf.*kpp.*delt./(vp*vy);
% lam2
lam2 = @(ks) ks./vy;
% lam3
lam3 = @(lam1,lam2) sqrt(lam1.^2-lam2.^2);
% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_yp = @(l1,l2,l3,x) kpy/vy.*(l2.*sinh(l3.*(Lhat-x)))/(l3*cosh(l3*Lhat)+l1*sinh(l3*Lhat));
rho_ym = @(l1,l2,l3,x) kpy/vy.*(l3^2.*cosh(l3*(Lhat-x))+l1*l3*sinh(l3*(Lhat-x)))./...
											(l3^2.*cosh(l3*Lhat)+l1*l3*sinh(l3*Lhat));
% acidification
c0 = 1;
% afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
afunc_sim = @(n,c0,rp,rl) n*c0.*rl^3/rp^3;
%%
dt = 0.001; % simulation timestep
phfrate = 10; % framerate for phagosomes (binwidth)
lyfrate = 1.01; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
%% alternative binning
nreg = 100;
xp = linspace(0,1,nreg+2);
xl = xp;
%%
%% plot organelle densities
for pc = 1:length(pvals)
	
	if(length(pvals)>1); subplot(1,3,pc); end
	
	for lc = 1:length(ellstr)
		
		nlyso = histcounts(lysopos{lc,pc},xl);
		nphago = histcounts(phpos{lc,pc},xp);
		
		xlvals = 0.5*(xl(1:end-1)+xl(2:end));
		% 	xlvals = xl(2:end);
		dxl = xlvals(2)-xlvals(1);
		
		xpvals = 0.5*(xp(1:end-1)+xp(2:end));
		% 	xpvals = xp(2:end);
		dxp = xpvals(2)-xpvals(1);
		
		l1val = lam1(ksvals(pc),kfvals(pc));
		l2val = lam2(ksvals(pc));
		l3val = lam3(l1val,l2val);
		
		rho_ypvals = rho_yp(l1val,l2val,l3val,xlvals);
		rho_ymvals = rho_ym(l1val,l2val,l3val,xlvals);
		
		rhovals = rho_ypvals+rho_ymvals;
		
		plot(xlvals,nlyso./ntrials./dxl,'bo','DisplayName','lyso density, sims')
		hold on
		plot(xlvals,rhovals,'b--','DisplayName','lyso density, an. model');
		plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
		plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
	end
	hold off
	plot_cleanup
	title(sprintf('pf: %0.3f, ell: %0.3f',pvals(pc),ellvals(pc)))
	xlim([0,0.2]);
	ylim([0,750]);
	xlabel('distance')
	ylabel('organelle density')
end
