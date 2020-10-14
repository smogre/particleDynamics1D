# make param files with reversals

# pvalstr=(0 0pt05 0pt10 0pt15 0pt20 0pt25 0pt30 0pt35 0pt40 0pt45 0pt50 0pt55 0pt60 0pt65 0pt70 0pt75 0pt80 0pt85 0pt90 0pt95 1)
# pvals=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1)
pvalstr=(0pt25 0pt50 0pt75)
pvals=(0.25 0.50 0.75)
runlenstr=(0pt01 0pt1 inf)
krev=(267 26.7 0)
for rc in {0..2}; do
for pc in {0..2}; do
  echo param.nfusedist_pf${pvalstr[$pc]}_runlen${runlenstr[$rc]}

  cp param.nfusedist_pf0pt25 param.nfusedist_pf${pvalstr[$pc]}_runlen${runlenstr[$rc]}
  
  sed -i "s/PFUSE 1 0D0 0.25D0/PFUSE 1 0D0 ${pvals[$pc]}D0/g" param.nfusedist_pf${pvalstr[$pc]}_runlen${runlenstr[$rc]};
  sed -i "s/PFUSE 2 0.25D0 0D0/PFUSE 2 ${pvals[$pc]}D0 0D0/g" param.nfusedist_pf${pvalstr[$pc]}_runlen${runlenstr[$rc]};
  sed -i "s/KREV 0D0 0D0/KREV ${krev[$rc]}D0 0D0/g" param.nfusedist_pf${pvalstr[$pc]}_runlen${runlenstr[$rc]};
done
done