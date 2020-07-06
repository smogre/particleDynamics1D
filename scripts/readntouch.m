% function [ntouchvals,ntrials] = readntouch(touchfile)
% reads out number of touches from the file
% output format: 
% 1st line: NTRIALS
%	2nd line onwards ntouch				        | time npartdist(1,nreg+1)
%%
clear;clc
filename = '../testntouch.touch.out';
fid = fopen(filename,'r');
siminfo = textscan(fid,'%f',1);
ntrials = siminfo{1};
ntouchvals = [];
while(~feof(fid))
	ntouchvals = [ntouchvals;cell2mat(textscan(fid,'%f',1))];
end
fclose(fid);

% end