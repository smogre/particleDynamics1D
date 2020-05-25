%% test textscan

fname = 'tmp.txt';
fid = fopen(fname,'r');
particles = struct('pos',[]);
tc = 1;
pc = 1;
while(~feof(fid))
	[pc tc]
% 	tmp = cat(1,tmp,textscan(fid,'%f',1));
 	tmp = textscan(fid,'%f',1);
	if(~isempty(tmp{:}));	particles(pc).pos(tc) = tmp{:}; end
	tc = tc+1;
	if(tc==6)
		pc = pc+1;
		tc = 1;
	end
end
fclose(fid);
%%

