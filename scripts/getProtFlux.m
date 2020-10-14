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


if(exist('options','var')); opt = copyStruct(options,opt); end

vy = opt.vy;
vp = opt.vp;
rp = opt.rp;
ry = opt.ry;
kpp = opt.kpp;
kpy = prodratio*opt.kpp;

xvals = xvals./xvals(end);
dx = xvals(2)-xvals(1);
Lhat = xvals(end);

% phagosome density
phdist = kpp/vp.*ones(size(xvals));
% lysosome density
lysodist = kpy/vy.*exp(-beta.*xvals);

% number of fusion events
nfuse = @(beta,prodratio,x) prodratio.*(exp(-beta.*x) - exp(-beta*Lhat));
% phagosome acidification
acidity = @(beta,prodratio,x) nfuse(beta,prodratio,x)./(nfuse(beta,prodratio,x)+rp^3/ry^3);

% number of fusions within the domain
nfusevals = nfuse(beta,prodratio,xvals);
% acidity within the domain
avals = acidity(beta,prodratio,xvals);

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