# make param files with reversals

pvals="0 0pt25 0pt50 0pt75 1"
for pf in $pvals; do
  echo param.nfusedist_pf${pf}
  cp param.nfusedist_pf${pf} param.nfusedist_pf${pf}_krev$1
  sed -i "s/KREV 0D0 0D0/KREV $1 0D0/g" param.nfusedist_pf${pf}_krev$1;
done