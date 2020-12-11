function [degflux,lysodist,phdist,nfusevals,avals,alldegflux] = ...
																									getProtFlux(beta,prodratio,kmax,KM,xvals,options)
% returns the degraded protein flux for a given parameter set
% assume that the phagosome starts at xvals(end) with zero acidification
% and single unit of undegraded protein
% default values for a 100 micron domain are set in options
% x values are automatically scaled such that xvals(end) = 1

opt = struct();
% phagosome production
opt.kpp = 1;
% phagosome velocity
opt.vp = 1;
% phagosome radius
opt.rp = 4e-3;
% lysosome velocity
opt.vy = 2.67;
% lysosome radius
opt.ry = 1e-3;
% use volumetric dilution to update acidification?
% if this option is chosen false, the acidification only changes due to
% vATPase pump activity
opt.dilute = true;

% michaelis-menten rates for vATPase pump activity
opt.kmax_a = 1;
opt.KM_a = 1;


if(exist('options','var')); opt = copyStruct(options,opt); end

vy = opt.vy;
vp = opt.vp;
rp = opt.rp;
ry = opt.ry;
kpp = opt.kpp;
kpy = prodratio*opt.kpp;
kmax_a = opt.kmax_a;
KM_a = opt.KM_a;
A0 = 1;

xvals = xvals./xvals(end);
dx = xvals(2)-xvals(1);
Lhat = xvals(end);

% phagosome density
phdist = kpp/vp.*ones(size(xvals));
% lysosome density
lysodist = kpy/vy.*exp(-beta.*xvals);

% number of fusion events per phagosome
nfuse = @(beta,prodratio,x) prodratio.*(exp(-beta.*x) - exp(-beta*Lhat));
nfusevals = nfuse(beta,prodratio,xvals);

% phagosome acidification
if(opt.dilute)
	% acidification changes due to dilution
	acidity = @(beta,prodratio,x) nfuse(beta,prodratio,x)./(nfuse(beta,prodratio,x)+rp^3/ry^3);
	avals = acidity(beta,prodratio,xvals);
else
	% acidification changes due to vATPase activity
	% acidification rate decreases as acidity of the organelle increases
	% ka(N,a) = kmax_a*N/(N+KM_a)*(A0-A)
	% ka: acidification rate
	% kmax_a, KM_a: Michaelis-Menten constants
	% N: number of fusion events
	% A0: maximal acidity (normalized to 1 in this case)
	% A: current acidity
	% solving for acidification, we get:
	% A(x) = A0*(1-exp(kmax_a./v.*int_L^x(N/(N+KM_a))))
	
	% integral of N/(N+KM_a) 
	nintvals = nfusevals./(nfusevals+KM_a);
	nintx = (cumsum(nintvals)-0.5*(nintvals(1)+nintvals))*dx; %integral from 0 to x
	nint1x = nintx-nintx(end); % integral from L to x
	
	avals = A0*(1-exp(kmax_a./vp.*nint1x));
end

% degradation rate within the domain
degrate = kmax*avals./(avals+KM)./vp;

% integral of degradation rate int_0^x
degintx = (cumsum(degrate)-0.5*degrate(1)-0.5*degrate)*dx;

% int_0^x-int_0^1 = -int_x^1
degintx1 = degintx-degintx(end);

% degraded protein flux
alldegflux = vp*(1-exp(degintx1));
degflux = alldegflux(1);

end