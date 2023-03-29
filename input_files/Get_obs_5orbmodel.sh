
echo "##{Jhund}  {lambda} {delta}  {L2_}  {S2_}  {J2_}  {n_t2g}  {n_eg}"

for Jhund in 0.25 
do

for lambda in 0.000 0.005 0.010 0.015 0.020 0.030
do

for delta in 0.0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 0.3 0.4 0.42 0.44 0.46 0.48 0.5 1.0
do



line_no=$(awk '/OR STATE NO \(EXACT\) 0/ { print NR; exit}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)

line=$(echo "${line_no} + 2" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n0_up=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 3" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n1_up=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 4" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n2_up=$(printf "%1.8f" ${temp3})


line=$(echo "${line_no} + 5" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n3_up=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 6" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n4_up=$(printf "%1.8f" ${temp3})


line=$(echo "${line_no} + 7" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n0_dn=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 8" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n1_dn=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 9" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n2_dn=$(printf "%1.8f" ${temp3})


line=$(echo "${line_no} + 10" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n3_dn=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 11" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $2}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
n4_dn=$(printf "%1.8f" ${temp3})




line=$(echo "${line_no} + 64" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
LzLz=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 64 + 3" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
LpLm=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 64 + 6" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
LmLp=$(printf "%1.8f" ${temp3})

L2_=$(echo "${LzLz} + 0.5*(${LpLm} + ${LmLp})" | bc -l)


line=$(echo "${line_no} + 73" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
SzSz=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 73 + 3" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
SpSm=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 73 + 6" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
SmSp=$(printf "%1.8f" ${temp3})

S2_=$(echo "${SzSz} + 0.5*(${SpSm} + ${SmSp})" | bc -l)


line=$(echo "${line_no} + 109" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
JzJz=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 109 + 3" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
JpJm=$(printf "%1.8f" ${temp3})

line=$(echo "${line_no} + 109 + 6" | bc -l)
temp1=$(awk -v row=${line} 'NR==row {print $1}' U_1.0/Jh_by_U_${Jhund}/Delta_by_U_${delta}/lambda_by_U_${lambda}/out_run.txt)
temp2=${temp1%,*}
temp3=$(echo "${temp2}" | sed 's/(//')
JmJp=$(printf "%1.8f" ${temp3})

J2_=$(echo "${JzJz} + 0.5*(${JpJm} + ${JmJp})" | bc -l)




n_t2g=$(echo "${n0_up} + ${n1_up} + ${n2_up} + ${n0_dn} + ${n1_dn} + ${n2_dn}" | bc -l)
n_eg=$(echo "${n3_up} + ${n4_up} + ${n3_dn} + ${n4_dn}" | bc -l)
echo "${Jhund}  ${lambda}  ${delta}  ${L2_}  ${S2_}  ${J2_}  ${n_t2g}  ${n_eg}"


done #delta
done #lambda
done #Jhund
