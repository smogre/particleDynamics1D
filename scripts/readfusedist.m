function [tvals,fusedist,nreg,regedges,ntrials] = readfusedist(filename)
% reads out the protein content from filename
% output format: 
% 1st line: NREG NTRIALS DOMLEN
%	2nd line: regpos(nreg)
% 3nd line onwards: time nfusedist(1,nreg+1)

fid = fopen(filename,'r');
siminfo = textscan(fid,'%f %f %f',1);
nreg = siminfo{1};
ntrials = siminfo{2};
domlen = siminfo{3};
regstr = strtrim(repmat('%f ',1,nreg));

regpos = textscan(fid,regstr,1);
regedges = [0,cat(2,regpos{:}),domlen];
fmtstr = strtrim(repmat('%f ',1,nreg+2));
data = textscan(fid,fmtstr);
tvals = data{1};
fusedist = cat(2,data{2:end});

fclose(fid);
end

