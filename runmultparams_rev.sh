 #run param files with reversal
pvals="0 0pt25 0pt50 0pt75 1"
for pf in $pvals; do
  echo param.nfusedist_pf${pf}
  ./partdynamics1D.exe nfusedist_pf${pf}_runlen$1
done