# make param files 

pvalstr=(0 0pt25 0pt50 0pt75 1)
pvals=(0 0.25 0.50 0.75 1)

for pc in {0..4}; do
  if [ $pc -gt 0 ]
  then
    cp param.nfusedist_pf0 param.nfusedist_pf${pvalstr[$pc]}
    echo param.nfusedist_pf${pvalstr[$pc]}
  
    sed -i "s/DOFUSION F/DOFUSION T/g" param.nfusedist_pf${pvalstr[$pc]};
    sed -i "s/PFUSE 1 0D0 0D0/PFUSE 1 0D0 ${pvals[$pc]}D0/g" param.nfusedist_pf${pvalstr[$pc]};
    sed -i "s/PFUSE 2 0D0 0D0/PFUSE 2 ${pvals[$pc]}D0 0D0/g" param.nfusedist_pf${pvalstr[$pc]};
    # sed -i "s/DELT 1D-3/DELT 1D-4/g" param.nfusedist_pf${pvalstr[$pc]};
  fi
  
done