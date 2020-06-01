classdef partgroupObj < handle
	% object to store properties of groups of particles
	properties
		% number of particle types
		ntype = 2;
		
		% number of particles
		npart = 0;
		
		% list of particles
		particles = [];
		
	end
	methods
		function pgroup = partgroupObj
% 			if(nargin~=0)
% 				pgroup(ntrials,1) = pgroup;
% 				for tc = 1:ntrials
% 					disp tc
					pgroup.npart = [];
					pgroup.ntype = 2;
					pgroup.clearParticles;
% 				end
% 			end
		end
		
		function pgroup = clearParticles(pgroup,varargin)
			pgroup.particles = struct('time',[],'pos',[],'id',[],'prot',[],'nfuse',[],'type',[]);
		end
	end
end