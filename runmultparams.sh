# for j in {2..5}; do 
#   if [ $1 -ge ${j} ]
#   then
#     echo param.fixtotlen${j}x_$1mt; 
#     ./transportsimMT.exe fixtotlen${j}x_$1mt;  
#   fi
# done


for i in $(seq 3 2 9); do
  echo param.$1_${i}mt
  ./transportsimMT.exe $1_${i}mt 
done