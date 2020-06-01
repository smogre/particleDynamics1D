%% calculate acidification vs number of fusion events
c0 = 1;
rp = 400;
rl = 25;

nvals = 1:10000;
cfunc = @(n,c0,rp,rl) n*c0./(n+rp.^3/rl^3);
%%
rvals = [1,4,8,12,16];
nvals = logspace(0,4,1e3);
cmap = colormapinterp([0,0,1;1,0,0],length(rvals));
% cmap = lines(length(rvals));
for rc = 1:length(rvals)
	semilogx(nvals,cfunc(nvals,c0,rvals(rc),1),'Color',cmap(rc,:),...
										'DisplayName',sprintf('autophagosome radius = %d x lysosome radius',rvals(rc)))
	hold on
end
plot(nvals,ones(size(nvals)),'k--','DisplayName','lysosome acidification level')
hold off
plot_cleanup
ylim([0,1.1])
xlabel('number of fusion events')
ylabel('relative acidification')
