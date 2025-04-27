# !/bin/bash
case_name="cavity_smac_09"

for i in {0..25}
do
    time=$(( i*10000 ))
    ./postProcess -c ${case_name} -t ${time} >> ${case_name}/energy.dat
done
