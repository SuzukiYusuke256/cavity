# !/bin/bash

# cavity
# start="0"
# end="250000"
# interval="10000"
# case_name="cavity_smac_09"

# couette
start="0"
end="2500000"
interval="100000"
case_name="couette_smac_04"

rm ${case_name}/energy.dat
touch ${case_name}/energy.dat
echo "TimeStep KineticEnergy ViscousDisspition WallWork massBalance" > ${case_name}/energy.dat

# for ( i=start; i<=end; i++ ))
for time in $(seq ${start} ${interval} ${end}) 
do
    ./postProcess -c ${case_name} -t ${time} >> ${case_name}/energy.dat
done
