%% test readsnapshot, dependence of protein content on degradation rate

% func = @(Jp,k,t) Jp.*(1-exp(-k*t))./k;
func = @(Jp,k,t) Jp.*(1-exp(-k*t));


%%
tvals = linspace(0,1,100);
kvals = logspace(-3,2,101)';
Jp = 1;
fluxout = func(Jp,kvals,tvals);
%%
pcolor(tvals,log10(kvals),log10(fluxout))
% pcolor(tvals,kvals,fluxout)
shading flat
plot_cleanup
xlabel('time')
ylabel('decay rate')
colorbar('TickLabelInterpreter','latex')
%%
semilogx(kvals,fluxout(:,end))
plot_cleanup
xlim([min(kvals),max(kvals)])
%%
xvals = linspace(0,1e1,1e3);
plot(xvals,lognpdf(xvals))
plot_cleanup
xlabel('time since start of quarantine')
ylabel('time spent playing poker/ludo/skribbl')
