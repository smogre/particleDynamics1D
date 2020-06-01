function [grouplist,tvals,domlen,ntrials] = readsnapshot_old(filename,options)
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
	if(~exist('indlist','var')); indlist = cell(ntrials,1); end

	% start reading all trials for a single time step
	trc = 1;
	while(trc<=ntrials)
		trinfo = textscan(fid,'%f %f',1); % trialnumber npart
		tc = trinfo{1};
		npart = trinfo{2};
		grouplist(tc).npart = cat(1,grouplist(tc).npart,npart);
		
		% read all particles for a given trial
		pc = 1;
		while(pc<=trinfo{2})
			pinfo = textscan(fid,'%f %f %f %f %f',1); %pc ID type nprot nfuse
			pid = pinfo{2};
			
			% read position and protein content
			fmtstr = strtrim(repmat('%f ',1,pinfo{4}+1));
			pdata = textscan(fid,fmtstr,1); %pos prot

			% if particle is not the type required, continue
			if(opt.onlyread~=0 && pinfo{3}~=opt.onlyread)
				pc = pc+1;
				continue
			end
			
			% if number of particles exceeded, continue
			if(length(grouplist(tc).particles)>=opt.maxpart)
				pc = pc+1;
				continue
			end
			
			% find location of particle in the index list for that trial
			[~,partind] = find(indlist{trc}==pid,1);
			if(isempty(partind))
				indlist{trc} = [indlist{trc},pid];
				partind = length(indlist{trc});
			end
			grouplist(tc).particles(partind).id = pid;
			grouplist(tc).particles(partind).type = pinfo{3};
			grouplist(tc).particles(partind).time = [grouplist(tc).particles(partind).time;initdata{1}];
			grouplist(tc).particles(partind).pos = [grouplist(tc).particles(partind).pos;pdata{1}];
			if(pinfo{4}~=0)
				grouplist(tc).particles(partind).prot = [grouplist(tc).particles(partind).prot;pdata{2:end}];
			end
			grouplist(tc).particles(partind).nfuse = pinfo{5};
			pc = pc+1;
		end
		
		trc = trc+1;
	end
end

fclose(fid);
end

