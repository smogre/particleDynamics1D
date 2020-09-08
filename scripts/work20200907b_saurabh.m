%%
errval = zeros(5,1);
rc = 1;
ct = 1;
acratio = zeros(length(pvals),length(xp)-1);
acdiff = zeros(length(pvals),length(xp)-1);
for pc = 1:5
	if(length(pvals)>1); subplot(2,3,ct); end
	[nphago,~,binind] = histcounts(phpos{pc,rc},xp);
	acvals_calc = afunc_sim(nfusedist{pc,rc},c0,rp,ry);
	acvals_dir = phacid{pc,rc};
	
	acdist_calc = accumarray(binind,acvals_calc);
	acdist_dir = accumarray(binind,acvals_dir);

	xpvals = 0.5*(xp(1:end-1)+xp(2:end)); 
	dxp = xpvals(2)-xpvals(1);
	
	acdist_th = afunc_th(c0,gamvals(pc),betvals(pc),xpvals);
	
	plot(xpvals,acdist_calc./nphago','o','Color',cmap(ct,:),'DisplayName','acidification, calc')
	hold on
	plot(xpvals,acdist_dir./nphago','--','Color',cmap(ct,:),...
																													'DisplayName','acidification, dir');
	plot(xpvals,acdist_th','^','Color',cmap(ct,:),...
																													'DisplayName','acidification, an');
	hold off
	plot_cleanup
	title(sprintf('fusion probability: %0.3f',pvals(pc)))
	xlim([0,0.2]);
	xlabel('distance')
	ylabel('acidification')
	acdiff(pc,:) = ((acdist_dir./nphago')-acdist_th');
	acratio(pc,:) = (acdist_dir./nphago')./acdist_th';
	ct = ct+1;
end