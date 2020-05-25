% function [tvals,cellprot,axprot,ntrials] = readsnapshot(filename)
% reads out the protein content from filename
% output format for each timestep: 
% Info line: TIME NTRIALS DOMLEN
%						|Ntrial blocks: TRIALNUMBER NPART
%						|							|Npart blocks: pc, ID, type, nprot
%						|							|							pos, prot

filename = '/data/SimData/particledynamics1D/snapfiles/test.snap.out';
fid = fopen(filename,'r');

tvals = [];
npart = [];
particles = struct('time',[],'pos',[],'id',[],prot,'[]');

while(~feof(fid))

	initdata = textscan(fid,'%f %f %f',1); % time ntrials domlen
	tvals = cat(1,tvals,initdata{1});
	ntrials = initdata{2};
	domlen = initdata{3};

	% start reading all trials for a single time step
	trc = 1;
	lastid = -inf(ntrials,1);
	while(trc<=ntrials)
		trinfo = textscan(fid,'%f %f',1); % trialnumber npart
		npart = cat(1,npart,trinfo{2});

		% read all particles for a given trial
		pc = 1;
		while(pc<=trinfo{2})
			pinfo = textscan(fid,'%f %f %f %f',1); %pc ID type nprot
			pid = pinfo{2};
			fmtstr = strtrim(repmat('%f ',1,pinfo{4}+1));
			pdata = textscan(fid,fmtstr,1); %pos prot
			particles(pid).pos = cat(1,particles(pid).pos,pdata{1});
			if(pinfo{4}~=0)
				particles(pid).prot = cat(1,particles(pid).prot,pdata{2:end});
			end
			pc = pc+1;
		end
		
		trc = trc+1;
	end
end
	
fclose(fid);
% end

