function outputSnapshotForPymol(outfile, microtubules,endosomes,alllinkers,peroxisomes,options)
% given coordinates for microtubues, endosomes, linkers, perox
% output to file for conversion to pdb and viewing in pymol

opt = struct();
% scale all coordinates before outputting to pdb
% default: 1 pdb unit = 0.1 snapshot unit = 100nm
opt.scl = 10; 
opt.domlen = range(microtubules{1}(:,3))/2;

if (exist('options','var'))
    opt= copyStruct(options,opt);
end

%% output to file
OF = fopen(outfile,'w');

% MT beads
for mct = 1:length(microtubules)
    MT = microtubules{mct};
    for b = 1:size(MT,1)
        fprintf(OF,'%s %f %f %f\n', sprintf('MT%d',mct), MT(b,:)*opt.scl);
    end
end

% endosomes and linkers
for ect = 1:length(endosomes)
    linkers = alllinkers{ect};
    endo = endosomes{ect};
    for lc = 1:length(alllinkers)
        linkbeads = linkers{lc};
        for b = 1:size(linkbeads,1)
            fprintf(OF,'%s %f %f %f\n', sprintf('L%d%d',ect,lc), linkbeads(b,:)*opt.scl);
        end
    end
    % endosome
    fprintf(OF,'%s %f %f %f\n', sprintf('EN%d',ect), endo*opt.scl);
end

% peroxisomes
for pc = 1:length(peroxisomes)
    perox = peroxisomes{pc};
    fprintf(OF,'%s %f %f %f\n', 'PX', perox*opt.scl);
end

% hypha cylinder
cylcent = [0,0,0];
cyl1 = [0,0,-opt.domlen];
cyl2 = [0,0,opt.domlen];
fprintf(OF,'%s %f %f %f\n', 'AC', cylcent*opt.scl);
fprintf(OF,'%s %f %f %f\n', 'AC1', cyl1*opt.scl);
fprintf(OF,'%s %f %f %f\n', 'AC2', cyl2*opt.scl);


fclose(OF);
end
