#!/bin/bash -x

for U in 1.0
do
mkdir -p U_${U}
cd U_${U}

for Jh_by_U in 0.25
do
mkdir -p Jh_by_U_${Jh_by_U}
cd Jh_by_U_${Jh_by_U}

for Delta_by_U in 0.42 0.44 0.46 0.48 0.3 0.4 0.5 1.0
#0.0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2
do
mkdir -p Delta_by_U_${Delta_by_U}
cd Delta_by_U_${Delta_by_U}


for lambda_by_U in 0.000 0.005 0.010 0.015 0.020 0.030
do
mkdir -p lambda_by_U_${lambda_by_U}
cd lambda_by_U_${lambda_by_U}


JHUND_TEMP=$(echo "${U}*${Jh_by_U}" | bc -l)
JHUND=$(printf "%1.6f" ${JHUND_TEMP})

UPRIME_TEMP=$(echo "$U - 2.0*${JHUND}" | bc -l) #  Uprime
UPRIME=$(printf "%1.6f" ${UPRIME_TEMP})

DELTA_TEMP=$(echo "${U}*${Delta_by_U}" | bc -l)
DELTA=$(printf "%1.6f" ${DELTA_TEMP})

LAMBDA_TEMP=$(echo "${U}*${lambda_by_U}" | bc -l)
LAMBDA=$(printf "%1.6f" ${LAMBDA_TEMP})


rm input.inp
rm lanczos
cp ../../../../inputMultiorb_temp.inp input.inp
cp ../../../../lanczos .

input="input.inp"

##########################################
sed -i -e "s/DELTA_VAL/${DELTA}/g" $input
sed -i -e "s/U_VAL/${U}/g" $input
sed -i -e "s/JHUND_VAL/${JHUND}/g" $input
sed -i -e "s/UPRIME_VAL/${UPRIME}/g" $input
sed -i -e "s/LAMBDA_VAL/${LAMBDA}/g" $input
#########################################

echo "submitted $U $Jh_by_U $Delta_by_U $lambda_by_U"

./lanczos $input > out_run.txt 


cd ..
done #lambda

cd ..
done #Delta


cd ..
done #Jh

cd ..
done #U

