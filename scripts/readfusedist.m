function [tvals,fusedist,partdist,nreg,regedges,ntrials] = readfusedist(filename)
% reads out the protein content from filename
% output format: 
% 1st line: NREG NTYPE NTRIALS DOMLEN
%	2nd line: regpos(nreg)
% 3nd line onwards, NTYPExNSAVE blocks: | time nfusedist(1,nreg+1)
%																        | time npartdist(1,nreg+1)
% clear;clc
% filename = '/data/SimData/particledynamics1D/nfusefiles/nfusedist_pf1.fuse.out';
fid = fopen(filename,'r');
siminfo = textscan(fid,'%f %f %f %f',1);
nreg = siminfo{1};
ntype = siminfo{2};
ntrials = siminfo{3};
domlen = siminfo{4};
regstr = strtrim(repmat('%f ',1,nreg));

regpos = textscan(fid,regstr,1);
regedges = [0,cat(2,regpos{:}),domlen];
fmtstr = strtrim(repmat('%f ',1,nreg+2));
data = textscan(fid,fmtstr);
tvals = data{1}(1:2*ntype:end);
%%
fusedist = zeros(size(tvals,1),nreg+1,ntype);
partdist = zeros(size(tvals,1),nreg+1,ntype);
for cc = 1:nreg+1
	for nc = 1:ntype
		fusedist(:,cc,nc) = data{cc+1}(2*nc-1:2*ntype:end);
		partdist(:,cc,nc) = data{cc+1}(2*nc:2*ntype:end);
	end
end
fclose(fid);
end

