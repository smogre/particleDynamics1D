%% test algorithm for producing poisson distribution
nvals = 1e6;
lambda = 1.8;
truedist = poissrnd(lambda,nvals,1);
gendist = getPoissonRnd(lambda,nvals);
[truevals,trueedges] = histcounts(truedist,'Normalization','probability');
[genvals,genedges] = histcounts(gendist,'Normalization','probability');
plot(0.5*(trueedges(1:end-1)+trueedges(2:end)),truevals,'ro-')
hold on
plot(0.5*(genedges(1:end-1)+genedges(2:end)),genvals,'b*-')
hold off