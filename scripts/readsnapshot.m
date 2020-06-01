function [grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options)
% reads out the protein content from filename
% output format for each timestep:
% Info line: TIME NTRIALS DOMLEN
%						|Ntrial blocks: TRIALNUMBER NPART
%						|							|Npart blocks: pc, ID, type, nprot, nfuse
%						|							|							pos, prot

opt = struct();
% output progress after how many steps
opt.progstep = 1;
% only read a specific particle type
opt.onlyread = 0;
% maximum number of particles to read
opt.maxpart = inf;

if(exist('options','var')); opt = copyStruct(options,opt); end

fid = fopen(filename,'r');
tvals = [];
stepct = 0;

while(~feof(fid))
	stepct = stepct+1;
	if(~mod(stepct,opt.progstep));  disp(stepct); end
	initdata = textscan(fid,'%f %f %f',1); % time ntrials domlen
	tvals = cat(1,tvals,initdata{1});
	ntrials = initdata{2};
	domlen = initdata{3};
	if(~exist('grouplist','var')); grouplist(ntrials,1) = partgroupObj; end

	% start reading all trials for a single time step
	trc = 1;
	while(trc<=ntrials)
		trinfo = textscan(fid,'%f %f',1); % trialnumber npart
		tc = trinfo{1};
		npart = trinfo{2};
		grouplist(tc).npart = cat(1,grouplist(tc).npart,npart);
		
		% read all particles for a given trial
		pdata = textscan(fid,'%f %f %f %f %f',2*npart); %pc ID type nprot nfuse
																									 %pos prot
		pidvals = pdata{2}(1:2:end);
		typevals = pdata{3}(1:2:end);
		nfusevals = pdata{5}(1:2:end);
		posvals = pdata{1}(2:2:end);
		protvals = [pdata{2}(2:2:end),pdata{3}(2:2:end)];
		
		for pc = 1:length(pidvals)
			pid = pidvals(pc);
			grouplist(tc).particles(pid).id = pid;
			grouplist(tc).particles(pid).type = typevals(pc);
			grouplist(tc).particles(pid).time = [grouplist(tc).particles(pid).time;initdata{1}];
			grouplist(tc).particles(pid).pos = [grouplist(tc).particles(pid).pos;posvals(pc)];
			grouplist(tc).particles(pid).prot = [grouplist(tc).particles(pid).prot;protvals(pc,:)];
			grouplist(tc).particles(pid).nfuse = nfusevals(pc);
		end
		
		trc = trc+1;
	end
end

fclose(fid);
end

