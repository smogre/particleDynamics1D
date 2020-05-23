for i in $(seq 0 10); do
#	for j in $(seq 1 10); do
	#for j in 9; do
	echo param.kfuse_${i}
#		cp param.contactprob_${i}endo_${j}MT_diff param.contactprob_${i}endo_${j}MT_pexrad300nm_diff
#		cp param.contactprob_${i}endo_${j}MT_teth param.contactprob_${i}endo_${j}MT_pexrad300nm_teth
#	cp ../param.test param.kfuse_${1}
  echo $((0.1*$i))
#		sed -i "s/NPART 6/NPART $(($j+1))/g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/NTRACK 5/NTRACK ${i}/g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/BDSTEPS 5D5/BDSTEPS 5D4/g" param.contactprob_${j}endo_${i}MT_teth;
#		sed -i "s/BDPRINTEVERY 1D5/BDPRINTEVERY 1D4/g" param.contactprob_${j}endo_${i}MT_teth;		
#		sed -i "s/WALKINDS 1 2 3/WALKINDS $(echo $(seq 1 $j))/g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/TETHINDS 4/TETHINDS $(($j+1))/g" param.contactprob_${j}endo_${i}MT_teth;
#		sed -i "s/\#TETHINDS/TETHINDS/g" param.contactprob_${j}endo_${i}MT_teth;
#		sed -i "s/NTRIALS 10/NTRIALS 2500/g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/SNAPSHOTS 10 /SNAPSHOTS 10 \/data\/Condor\//g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/\/data\/Condor\/spheretouchBD\/\/data\/Condor\//\/data\/Condor\/spheretouchBD\//g" param.contactprob_${j}endo_${i}MT_*;
#		sed -i "s/STERRAD 0.1D0/STERRAD 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.1D0 0.3D0/g" param.contactprob_${i}endo_${j}MT_pexrad300nm*;
#	done
done
