function [tvals,cellprot,axprot,ntrials] = readprot(filename)
% reads out the protein content from filename
% output format: 
% 1st line: NPROT NTRIALS
% 2nd line onwards: time PROTCONT(1,nprot)

fid = fopen(filename,'r');
siminfo = textscan(fid,'%f %f',1);
nprot = siminfo{1};
ntrials = siminfo{2};
fmtstr = strtrim(repmat('%f ',1,2*nprot+1));

data = textscan(fid,fmtstr);
tvals = data{1};
cellprot = zeros(length(tvals),nprot);
axprot = zeros(length(tvals),nprot);
for pc = 1:nprot
	cellprot(:,pc) = data{pc+1};
	axprot(:,pc) = data{pc+nprot+1};
end
fclose(fid);
end

