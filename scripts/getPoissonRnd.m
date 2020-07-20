function rndvals = getPoissonRnd(lambda,nvals)

rndvals = zeros(nvals,1);

expval = exp(-lambda);
for nc = 1:nvals
	uprod = rand;
	nct = 0;
	while(uprod>expval)
		uprod = uprod*rand;
		nct = nct+1;
	end
	rndvals(nc) = nct;
end

end