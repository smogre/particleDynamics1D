% read in data
filename = '../results/nfusedist_pf1.fuse.out';
% partdist has dimensions: time, position, particle
% what exactly are regedges??
[tvals,fusedist,partdist,nreg,regedges,ntrials] = readfusedist(filename);
length(tvals)
%%
nbins = nreg+1;
binlen = 1/nbins;
Lum = 100;
% particle distributions per um
for tc = 1:length(tvals)
    plot(regedges(1:end-1)*Lum,partdist(tc,:,2)/ntrials/binlen/Lum)
    hold all
    plot(regedges(1:end-1)*Lum,partdist(tc,:,1)/ntrials/binlen/Lum)
    hold off    
    title(sprintf('Time %f',tvals(tc)))
    %ylim([0,0.8])
    drawnow
    pause(0.1)    
end

%% read snapshot
filename = '../results/test.snap.out';
options = struct('progstep',100);
% grouplist.particles(i) has trajectory of particle i
% tvals = saved timepoints
% grouplist.npart = total number of particles at each timepoint
[grouplist,tvals,domlen,ntrials] = readsnapshot(filename,options);

%%
ptypes = [grouplist.particles.type];
ids = [grouplist.particles.id];

for pc = 1:length(grouplist.particles)
    starttimes(pc) = grouplist.particles(pc).time(1);
    endtimes(pc) = grouplist.particles(pc).time(end);
end


%% visualize movie over time
for fc= 603%:length(tvals)
    curtime = tvals(fc);
    % particles present at this time
    ind = find(starttimes<=curtime & endtimes>=curtime)
    for ic = 1:length(ind)
        pc = ind(ic);
        particle = grouplist.particles(pc);
        cc = find(particle.time == curtime);
        allpos(ic) = particle.pos(cc);
        ptype(ic) = particle.type;
        ids(ic) = particle.id;
    end
    
    ind1 = find(ptype==1); ind2 = find(ptype==2);
    plot(allpos(ind1),zeros(size(ind1)),'bo')
    hold all
    plot(allpos(ind2),zeros(size(ind2)),'ro')
end
hold off