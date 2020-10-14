%% plot organelle densities with reversal
clear;clc;
pvals = [0.25,0.5,0.75];%[0,0.25,0.5,0.75,1];
pstr = [{'0pt25'},{'0pt50'},{'0pt75'}];
runlenvals = [0.01,0.1,inf];%0.1;
runlenstr =[{'0pt01'},{'0pt1'},{'inf'}];
%%
tic
cmap = [0,0,1;0,0,1];
options = struct();
options.readlastsnap = 1;
lysopos = cell(length(runlenstr),length(pvals));
phpos = lysopos;
phacid = lysopos;
nfusedist = lysopos;
cellprot = lysopos;
axprot = lysopos;
protdist = lysopos;

for rc = 1:length(runlenstr)
	for pc = 1:length(pvals)
		disp([rc pc])
		filename = sprintf('/data/SimData/particledynamics1D/snapfiles/nfusedist_pf%s_runlen%s.snap.out',...
																																							pstr{pc},runlenstr{rc});
		[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);
		
		protfile = sprintf('/data/SimData/particledynamics1D/outfiles/nfusedist_pf%s_runlen%s.out',...
			pstr{pc},runlenstr{rc});
% 		protfile = sprintf('../nfusedist_pf%s.out',pstr(pc));
		[~,cellprot{pc,rc},axprot{pc,rc},~] = readprot(protfile);
		
		tmplysopos = [];
		tmpphpos = [];
		tmpphacid = [];
		tmpnfuse = [];
		tmpprot = [];
		
		for tc = 1:length(grouplist)
			pos = vertcat(grouplist(tc).particles(:).pos);
			acid = vertcat(grouplist(tc).particles(:).acidity);
			ptypes = vertcat(grouplist(tc).particles(:).type);
			nfuse = vertcat(grouplist(tc).particles(:).nfuse);
			prot = vertcat(grouplist(tc).particles(:).prot);
			
			tmplysopos = cat(1,tmplysopos,pos(ptypes==1));
			tmpphpos = cat(1,tmpphpos,pos(ptypes==2));
			tmpphacid = cat(1,tmpphacid,acid(ptypes==2));
			tmpnfuse = cat(1,tmpnfuse,nfuse(ptypes==2));
			tmpprot = cat(1,tmpprot,prot(ptypes==2,:));
		end
		
		lysopos{rc,pc} = tmplysopos;
		phpos{rc,pc} = tmpphpos;
		phacid{rc,pc} = tmpphacid;
		nfusedist{rc,pc} = tmpnfuse;
		protdist{rc,pc} = tmpprot;
		
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
ksvals = vy./runlenvals;
% lam1
lam1 = @(ks,kf) ks./vy+kf.*kpp.*delt./(vp*vy);
% lam2
lam2 = @(ks) ks./vy;
% gam
gam = @(lam1,lam2) sqrt(lam1.^2-lam2.^2);
% theta
thet = @(l1,l2) acos(l2/l1);

% density
rho_p = @(x) kpp/vp.*ones(size(x));
rho_yp = @(thet,gam,x) kpy/vy.*(cos(thet).*sinh(gam.*(Lhat-x)))./...
																											 (sin(thet).*cosh(gam.*Lhat)+sinh(gam.*Lhat));
rho_ym = @(thet,gam,x) kpy/vy.*(sinh(gam.*(Lhat-x))+sin(thet).*cosh(gam.*(Lhat-x)))./...
																											 (sin(thet).*cosh(gam.*Lhat)+sinh(gam.*Lhat));
% acidification from sims
c0 = 1;
afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
% number of fusion events
Nfunc_th = @(thet,gam,x) (csc(thet)-cot(thet)).*kpy./kpp.*((1+cos(thet)).*(cosh(gam.*(Lhat-x))-1)+...
						sin(thet).*sinh(gam.*(Lhat-x)))./(sin(thet).*cosh(gam.*Lhat)+sinh(gam.*Lhat));
% acidification with increasing volume
afunc_th = @(thet,gam,x) c0*Nfunc_th(thet,gam,x)./(Nfunc_th(thet,gam,x) + rp^3/ry^3);
% amount of degraded protein near cell body
k0 = 1;
kM = 0.01;					
%%
dt = 0.001; % simulation timestep
phfrate = 8; % framerate for phagosomes (binwidth)
lyfrate = 1.1; % framerate for lysosomes (binwidth)
xp = 0:(phfrate*vp*dt):1; % possible positions for phagosomes
xl = 0:(lyfrate*vy*dt):1; % possible positions for lysosomes
%% alternative binning
nreg = 370;
xp = linspace(0,1,nreg+2);
xl = xp;
%% plot organelle densities
ct = 1;
for pc = 1:length(pvals)
	
	for rc = 1:length(runlenstr)
		if(length(pvals)>1); subplot(3,3,ct); end
		
		nlyso = histcounts(lysopos{rc,pc},xl);
		nphago = histcounts(phpos{rc,pc},xp);
		
		xlvals = 0.5*(xl(1:end-1)+xl(2:end));
		% 	xlvals = xl(2:end);
		dxl = xlvals(2)-xlvals(1);
		
		xpvals = 0.5*(xp(1:end-1)+xp(2:end));
		% 	xpvals = xp(2:end);
		dxp = xpvals(2)-xpvals(1);
		
		l1val = lam1(ksvals(rc),kfvals(pc));
		l2val = lam2(ksvals(rc));
		thetval = thet(l1val,l2val);
		gamval = gam(l1val,l2val);
		
		rho_ypvals = rho_yp(thetval,gamval,xlvals);
		rho_ymvals = rho_ym(thetval,gamval,xlvals);
		
		rhovals = rho_ypvals+rho_ymvals;
		
		plot(xlvals,nlyso./ntrials./dxl,'bo','DisplayName','lyso density, sims')
		hold on
		plot(xlvals,rhovals,'b--','DisplayName','lyso density, an. model');
		plot(xpvals,nphago./ntrials./dxp,'ro','DisplayName','phago density, sims')
		plot(xpvals,rho_p(xpvals),'r--','DisplayName','phago density, an. model');
		hold off
		plot_cleanup
		title(sprintf('pf: %0.3f, runlen: %0.3f',pvals(pc),runlenvals(rc)))
		xlim([0,0.1]);
		ylim([0,750]);
		xlabel('distance')
		ylabel('density')
		ct = ct+1;
	end
end
%% plot acidification
ct = 1;
cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
for pc = 1:length(pvals)
	
	for rc = 1:length(runlenstr)
		if(length(pvals)>1); subplot(3,3,ct); end
		
		[nphago,~,binind] = histcounts(phpos{rc,pc},xp);
		acvals = phacid{rc,pc};
		achist = accumarray(binind,acvals);
		
		xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
		dxp = xpvals(2)-xpvals(1);
				
		l1val = lam1(ksvals(rc),kfvals(pc));
		l2val = lam2(ksvals(rc));
		thetval = thet(l1val,l2val);
		gamval = gam(l1val,l2val);
		
	
		plot(xpvals,achist./nphago','o','Color',cmap(pc,:),'DisplayName','acidification, sims')
		hold on
		plot(xpvals,afunc_th(thetval,gamval,xpvals),'^--','Color',cmap(pc,:),...
																													'DisplayName','acidification, an. model');
		hold off
		plot_cleanup
		title(sprintf('pf: %0.3f, runlen: %0.3f',pvals(pc),runlenvals(rc)))
		xlim([0,0.1]);
% 		ylim([0,750]);
		xlabel('distance')
		ylabel('acidification')
		ct = ct+1;
	end
end
%% distribution of degraded protein
ct = 1;
cmap = colormapinterp([0,0,1;1,0,0],length(pvals));
xint = linspace(0,1,5000); dx = xint(2)-xint(1);
for pc = 1:length(pvals)
	
	for rc = 1:length(runlenstr)
		if(length(pvals)>1); subplot(3,3,ct); end
		
		[nphago,~,binind] = histcounts(phpos{rc,pc},xp);
		undegprot = accumarray(binind,protdist{rc,pc}(:,1));
		degprot = accumarray(binind,protdist{rc,pc}(:,2));
		
		xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
% 	xpvals = xp(2:end);
		dxp = xpvals(2)-xpvals(1);
				
		l1val = lam1(ksvals(rc),kfvals(pc));
		l2val = lam2(ksvals(rc));
		thetval = thet(l1val,l2val);
		gamval = gam(l1val,l2val);
		
		% an. model for deg prot
		A = afunc_th(thetval,gamval,xint);
		P = k0*A./(A+kM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
	
		plot(xpvals,degprot./nphago','o','Color',cmap(pc,:))
		hold on
		plot(xint,G,'--','Color',cmap(pc,:));
		hold off
		plot_cleanup
		title(sprintf('pf: %0.3f, runlen: %0.3f',pvals(pc),runlenvals(rc)))
		xlim([0,0.1]);
% 		ylim([0,750]);
		xlabel('distance')
		ylabel('deg prot')
		ct = ct+1;
	end
end
%% surface plot for beta vs switching rate
nvals = 1e2;
ellvals = [0,logspace(-4,0,nvals)];
ksvals = vy./ellvals;

kdmax = 1e0;
KM = 1e-1;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% beta
betvals = logspace(-2,0,nvals);

degflux = zeros(length(ellvals),nvals);

pvals = 0.99;

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
% lam1
lam1 = @(ks,kf,kpp) ks./vy+kf.*kpp.*delt./(vp*vy);
% lam2
lam2 = @(ks) ks./vy;
% gam
gam = @(lam1,lam2) sqrt(lam1.^2-lam2.^2);
% theta
thet = @(l1,l2) acos(l2/l1);
% number of fusion events
Nfunc_th = @(thet,gam,kpp,x) (csc(thet)-cot(thet)).*kpy./kpp.*((1+cos(thet)).*(cosh(gam.*(Lhat-x))-1)+...
						sin(thet).*sinh(gam.*(Lhat-x)))./(sin(thet).*cosh(gam.*Lhat)+sinh(gam.*Lhat));
% acidification with increasing volume
afunc_th = @(thet,gam,kpp,x) c0*Nfunc_th(thet,gam,kpp,x)./(Nfunc_th(thet,gam,kpp,x) + rp^3/ry^3);

for mc = 1:nvals
	
	disp(mc)
	
	% phago production rate
	kppval = betvals(mc)*vp*vy./delt./kfvals;
	
	parfor sc = 1:length(ksvals)
		
		l1val = lam1(ksvals(sc),kfvals,kppval);
		l2val = lam2(ksvals(sc));
		thetval = thet(l1val,l2val);
		gamval = gam(l1val,l2val);
		
		A = afunc_th(thetval,gamval,kppval,xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(sc,mc) = G(1)*vp;
	end
	
end
%% get surface plot for beta
pcolor(ellvals,betvals,degflux')
shading flat
colormap parula;
cb = colorbar;
cb.TickLabelInterpreter = 'latex';
cb.Label.Interpreter = 'latex';
cb.Label.String = 'degraded protein flux $J_d$';
cb.Label.Rotation = 270;
cb.Label.Units = 'normalized';
cb.Label.Position = [6,0.5,0];
set(gca,'XScale','log','YScale','log')
plot_cleanup
xlabel('run length')
ylabel('$\beta \equiv k_fk_{p,p}\delta/v_pv_y$')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$p_f = %0.4g$, $k_{p,y} = %0.4g$',...
								kdmax,KM,pvals,kpy))

%% surface plot for switching rate vs fusion probability
nvals = 1e2;
ellvals = logspace(-2,2,nvals);
ksvals = vy./ellvals;

kdmax = 1e0;
KM = 1e-1;
xint = linspace(0,1,5000); dx = xint(2)-xint(1);

% fusion probability
pvals = logspace(-3,0,nvals);

degflux = zeros(length(ellvals),nvals);

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
% lam1
lam1 = @(ks,kf) ks./vy+kf.*kpp.*delt./(vp*vy);
% lam2
lam2 = @(ks) ks./vy;
% gam
gam = @(lam1,lam2) sqrt(lam1.^2-lam2.^2);
% theta
thet = @(l1,l2) acos(l2/l1);
% acidification from sims
c0 = 1;
afunc_sim = @(n,c0,rp,rl) n*c0./(n+rp^3/rl^3);
% number of fusion events
Nfunc_th = @(thet,gam,x) (csc(thet)-cot(thet)).*kpy./kpp.*((1+cos(thet)).*(cosh(gam.*(Lhat-x))-1)+...
						sin(thet).*sinh(gam.*(Lhat-x)))./(sin(thet).*cosh(gam.*Lhat)+sinh(gam.*Lhat));
% acidification with increasing volume
afunc_th = @(thet,gam,x) c0*Nfunc_th(thet,gam,x)./(Nfunc_th(thet,gam,x) + rp^3/ry^3);

for mc = 1:nvals
	
	disp(mc)
	
	kfval = kfvals(mc);
	
	parfor sc = 1:length(ksvals)
		
		l1val = lam1(ksvals(sc),kfval);
		l2val = lam2(ksvals(sc));
		thetval = thet(l1val,l2val);
		gamval = gam(l1val,l2val);
		
		A = afunc_th(thetval,gamval,xint);
		P = kdmax*A./(A+KM)/vp;
		H = (cumsum(P) - 0.5*P(1)-0.5*P)*dx;
		% this is int_0^x - int_0^1 = - int_x^1
		Hx1 = H - H(end);
		G = 1 - exp(Hx1);
		
		degflux(sc,mc) = G(1)*vp;
	end
	
end
%% get surface plot for beta
pcolor(pvals,ellvals,degflux)
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
ylabel('run length')
xlabel('fusion probability')
title(sprintf('$k_\\mathrm{max} = %0.4g,~K_M = %0.4g$\n$k_{p,p} = %0.4g$, $k_{p,y} = %0.4g$',...
								kdmax,KM,kpp,kpy))
