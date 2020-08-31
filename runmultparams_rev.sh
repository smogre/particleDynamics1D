 #run param files with reversal
runlenstr=(0pt01 0pt1 inf)
for rc in {0..2}; do
  echo param.nfusedist_pf$1_runlen${runlenstr[$rc]}
  ./partdynamics1D.exe nfusedist_pf$1_runlen${runlenstr[$rc]}
done
