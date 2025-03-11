#!/bin/bash -x

for U in 50.0
do
mkdir -p U_${U}
cd U_${U}

for Jh_by_U in 0.000 0.020 0.050 0.080 0.100 0.120 0.150 0.180 0.200 0.210 0.215 0.220 0.225 0.230 0.235 0.240 0.242 0.244 0.246 0.248 0.250 0.252 0.254 0.256 0.258 0.260 0.265 0.270 0.275 0.280 0.285 0.290 0.295 0.300 0.305 0.310 0.315 0.320
do
mkdir -p Jh_by_U_${Jh_by_U}
cd Jh_by_U_${Jh_by_U}

JHUND_TEMP=$(echo "${U}*${Jh_by_U}" | bc -l)
JHUND=$(printf "%1.6f" ${JHUND_TEMP})

UPRIME_TEMP=$(echo "$U - 2.0*${JHUND}" | bc -l) #  Uprime
UPRIME=$(printf "%1.6f" ${UPRIME_TEMP})


rm input.inp
rm lanczos
cp ../../input_3orb.inp input.inp
cp ../../lanczos .

input="input.inp"

##########################################
sed -i -e "s/U_VAL/${U}/g" $input
sed -i -e "s/JHUND_VAL/${JHUND}/g" $input
sed -i -e "s/UPRIME_VAL/${UPRIME}/g" $input
#########################################

echo "submitted $U $Jh_by_U"

./lanczos $input > out_run.txt


cd ..
done

cd ..
done
