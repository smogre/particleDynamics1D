figure
for pc = 1:length(pvals)
	subplot(2,3,pc)
	[nphago,~,binind] = histcounts(phpos{pc},xp);
	nfusevals = nfusedist{pc}(binind==1);
	histogram(nfusevals)
	plot_cleanup
	xlabel('number of fusions')
	ylabel('number of phagosomes')
	title(sprintf('pf:%0.2f, mean fusions: %0.2f',pvals(pc),mean(nfusevals)));
end
