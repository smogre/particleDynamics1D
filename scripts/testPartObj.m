% function [tvals,cellprot,axprot,ntrials] = readsnapshot(filename)
% reads out the protein content from filename
% output format for each timestep: 
% Info line: TIME NTRIALS DOMLEN
%						|Ntrial blocks: TRIALNUMBER NPART
%						|							|Npart blocks: pc, ID, type, nprot, nfuse
%						|							|							pos, prot
clear;clc
filename = '/data/SimData/particledynamics1D/snapfiles/test.snap.out';
fid = fopen(filename,'r');

tvals = [];
npart = [];

trflag = 0;
%%
while(~feof(fid) && ~trflag)
	initdata = textscan(fid,'%f %f %f',1); % time ntrials domlen
	%%
	tvals = cat(1,tvals,initdata{1});
	ntrials = initdata{2};
	domlen = initdata{3};
	if(~exist('grouplist','var')); grouplist(ntrials) = partgroupObj; end
	% start reading all trials for a single time step
	trc = 1;
	lastid = -inf(ntrials,1);
	%%
	while(trc<=ntrials)
		trinfo = textscan(fid,'%f %f',1); % trialnumber npart
		tc = trinfo{1};
		npart = trinfo{2};
		%%
		grouplist(tc).npart = cat(1,grouplist(tc).npart,npart);
		% read all particles for a given trial
		pc = 1;
		%%
		while(pc<=trinfo{2})
			pinfo = textscan(fid,'%f %f %f %f %f',1); %pc ID type nprot nfuse
			pid = pinfo{2};
			%%
			if(pid>lastid(trc)); lastid(trc) = pid; end
			%%
			fmtstr = strtrim(repmat('%f ',1,pinfo{4}+1));
			pdata = textscan(fid,fmtstr,1); %pos prot
			%%
			grouplist(tc).particles(pid).id = pid;
			grouplist(tc).particles(pid).time = [grouplist(tc).particles(pid).time;initdata{1}];
			grouplist(tc).particles(pid).pos = [grouplist(tc).particles(pid).pos;pdata{1}];
			if(pinfo{4}~=0)
				grouplist(tc).particles(pid).prot = [grouplist(tc).particles(pid).prot;pdata{2:end}];
			end
			grouplist(tc).particles(pid).nfuse = pinfo{5};
			pc = pc+1;
		end
%%
		trc = trc+1;
	end
% 	trflag = 1;
% 	tmp = textscan(fid,'%f',1)
end
	%%
fclose(fid);
% end

