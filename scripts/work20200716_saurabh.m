%% check particle binning
nreg = 39;
vp = 0.001;
vy = 0.00267;
xbins = linspace(0,1,nreg+2);
xp = 0:20*vp:1;
xl = 0:3*vy:1;
%%
[nlyso,xlvals] = histcounts(lysopos,xl);
[nphago,xpvals] = histcounts(phpos,xp);
subplot(2,1,1)
histogram(lysopos,xl)
hold on
plot(0.5*(xlvals(1:end-1)+xlvals(2:end)),nlyso)
hold off
plot_cleanup
xlim([0,0.2])

subplot(2,1,2)
histogram(phpos,xp)
hold on
plot(0.5*(xpvals(1:end-1)+xpvals(2:end)),nphago)
hold off
plot_cleanup
xlim([0,0.2])
%%
[xpvals,xpedges] = histcounts(xp,xbins);
[xyvals,xyedges] = histcounts(xy,xbins);
subplot(2,1,1)
plot(xpedges(1:end-1),xpvals,'r');
plot_cleanup
title('phagosomes')
subplot(2,1,2)
plot(xyedges(1:end-1),xyvals,'b');
plot_cleanup
title('lysosomes')
%%
