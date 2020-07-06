%%
clear;clc;
L = 1;
vg = 1;
vp = 2.67;
kg = 1.6;
kl = 17.9;
func = @(t) exp(-kg*(L/vg+t+L/vp)).*(1-exp(-kl*t));
%%
tvals = linspace(0,1,1e3);
plot(tvals,func(tvals),'DisplayName','survival probability')
hold on
plot([1,1]/kl,ylim,'k--','DisplayName','average phagosome production time')
plot_cleanup
xlabel('time since last lysosome produced')
ylabel('probability of missing a phagosome')
hold off
