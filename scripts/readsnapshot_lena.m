function [grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options)
% reads out the protein content from filename
% output format for each timestep:
% Info line: TIME NTRIALS DOMLEN
%						|Ntrial blocks: TRIALNUMBER NPART
%						|							|Npart blocks: pc, ID, type, nprot, nfuse
%						|							|							pos, prot, acidity

opt = struct();
% output progress after how many steps
opt.progstep = 0;
% only read a specific particle type
opt.onlyread = 0;
% maximum number of particles to read
opt.maxpart = inf;
% only read last snapshot?
opt.readlastsnap = 0;

if(exist('options','var')); opt = copyStruct(options,opt); end

fid = fopen(filename,'r');
tvals = [];
stepct = 0;

while(~feof(fid))
	initdata = textscan(fid,'%f %f %f',1); % time ntrials domlen
	if(isempty(initdata{1})); return; end
	stepct = stepct+1;
	if(~mod(stepct,opt.progstep));  disp(stepct); end
	tvals = cat(1,tvals,initdata{1});
	if(~exist('ntrials','var')); ntrials = initdata{2}; end
	if(~exist('domlen','var')); domlen = initdata{3}; end
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
																									  %pos prot acidity
		if(npart<=0)
			trc = trc+1;
			continue;
		end
		pidvals = pdata{2}(1:2:end);
		typevals = pdata{3}(1:2:end);
		nfusevals = pdata{5}(1:2:end);
		posvals = pdata{1}(2:2:end);
		protvals = [pdata{2}(2:2:end),pdata{3}(2:2:end)];
		acidityvals = pdata{4}(2:2:end);
		
		if(opt.readlastsnap)
			grouplist(tc).npart = npart;
			grouplist(tc).particles = [];
			for pc = 1:length(pidvals)
				grouplist(tc).particles(pc).id = pidvals(pc);
				grouplist(tc).particles(pc).type = typevals(pc);
				grouplist(tc).particles(pc).time = initdata{1};
				grouplist(tc).particles(pc).pos = posvals(pc);
				grouplist(tc).particles(pc).prot = protvals(pc,:);
				grouplist(tc).particles(pc).nfuse = nfusevals(pc);
				grouplist(tc).particles(pc).acidity = acidityvals(pc);
			end
		else
			for pc = 1:length(pidvals)
				pid = pidvals(pc);
				grouplist(tc).particles(pid).id = pid;
				grouplist(tc).particles(pid).type = typevals(pc);
				grouplist(tc).particles(pid).time = [grouplist(tc).particles(pid).time;initdata{1}];
				grouplist(tc).particles(pid).pos = [grouplist(tc).particles(pid).pos;posvals(pc)];
				grouplist(tc).particles(pid).prot = [grouplist(tc).particles(pid).prot;protvals(pc,:)];
				grouplist(tc).particles(pid).nfuse = nfusevals(pc);
				grouplist(tc).particles(pid).acidity = acidityvals(pc);
			end
		end
		
		trc = trc+1;
	end
end

fclose(fid);
end

