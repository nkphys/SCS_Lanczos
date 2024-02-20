#!/bin/bash -x

for RAND in 321 
do

for L in 4
do
mkdir -p L_${L}
cd L_${L}

for N in 16
do
mkdir -p N_${N}
cd N_${N}

for U in 0.5 
do
mkdir -p U_${U}
cd U_${U}

JHUND="0.0"   #pair hopping term
UPRIME="0.0"
flot1="1.0"

JHUND=$(echo "scale=6;$U/4" | bc)  # J=U/4
if (( $(echo "$JHUND < $flot1" | bc -l) )); then
    JHUND=$(echo "0$JHUND")
fi

UPRIME=$(echo "scale=6;$U - 2.0*$JHUND" | bc) #  Uprime
if (( $(echo "$UPRIME < $flot1" | bc -l) )); then
    UPRIME=$(echo "0$UPRIME")
fi


#echo "${U}  ${UPRIME}   ${JHUND}"
for lambda in 0.0 0.1 0.2 0.3 0.4 0.6 0.8 0.9 1.0  #0.665 0.675 0.685 0.715 0.725 0.735 0.745 

do
mkdir -p lambda_${lambda}
cd lambda_${lambda}

rm input3orb_L${L}_SOC.inp
rm lanczos
cp ../../../../input3orb_L${L}_SOC.inp input3orb_L${L}_SOC_random${RAND}.inp
cp ../../../../lanczos .

input="input3orb_L${L}_SOC_random${RAND}.inp"

##########################################
sed -i -e "s/N_TOTAL_VALUE/${N}/g" $input
sed -i -e "s/U_VALUE/${U}/g" $input
sed -i -e "s/JHUND_VALUE/${JHUND}/g" $input
sed -i -e "s/UPRIME_VALUE/${UPRIME}/g" $input
sed -i -e "s/LAMBDA_SOC_VALUE/${lambda}/g" $input
sed -i -e "s/VALUE_RAND/${RAND}/g" $input
#########################################

echo "submitted $U $N $lambda"

./lanczos $input > out_${RAND}.txt 


cd ..
done #lambda

cd ..
done #U


cd ..
done #N

cd ..
done #L

done #RAND

