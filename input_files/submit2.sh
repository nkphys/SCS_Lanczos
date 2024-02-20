#!/bin/bash -x

#d_anis=1.0

for JHbyU in 0.25
do
mkdir -p JHbyU_${JHbyU}
cd JHbyU_${JHbyU}


for U in 0.0 0.1 0.2 0.5 0.8 1.0 2.0 3.0 4.0 5.0 6.0 6.4 7.0 8.0 8.5 9.0 9.5 10.0 10.5 11.0 12.0 13.0 14.0 15.0 16.0 18.0 19.0 20.0 22.0 24.0 26.0 28.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 
do
mkdir -p U_${U}
cd U_${U}


JHUND="0.0"   #pair hopping term
UPRIME="0.0"
flot1="1.0"

JHUND=$(echo "scale=6;$U*${JHbyU}" | bc)  # J=U/4
if (( $(echo "$JHUND < $flot1" | bc -l) )); then
    JHUND=$(echo "0$JHUND")
fi

UPRIME=$(echo "scale=6;$U - 2.0*$JHUND" | bc) #  Uprime
if (( $(echo "$UPRIME < $flot1" | bc -l) )); then
    UPRIME=$(echo "0$UPRIME")
fi

echo "${U}  ${UPRIME}   ${JHUND}"



rm input2orb.inp
rm lanczos
cp ../../input_2orb_template.inp input2orb.inp
cp ../../lanczos .
cp -rf ../../Variational_states_L6 .

input="input2orb.inp"

##########################################
sed -i -e "s/U_VALUE/${U}/g" $input
sed -i -e "s/JHUND_VALUE/${JHUND}/g" $input
sed -i -e "s/UPRIME_VALUE/${UPRIME}/g" $input
#########################################

rm "$input-e"
echo "submitted $U"

time ./lanczos $input > out.txt 




cd ..
done #U


cd ..
done #N
