#!/bin/bash -x

#d_anis=1.0

for d_anis in 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.085 0.9 0.95 1.0
do
mkdir -p d_anis${d_anis}
cd d_anis${d_anis}

for L in 6
do
mkdir -p L${L}_OBC
cd L${L}_OBC

#0.0 0.1 0.2 0.5 0.8 1.0 2.0 3.0 4.0 5.0 6.0 6.4 7.0 8.0 8.5 9.0 9.5 10.0 10.5 11.0 12.0 13.0 14.0 15.0 16.0 18.0 19.0 20.0 22.0 24.0 26.0 28.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0
for U in 0.0 0.1 0.2 0.5 0.8 1.0 2.0 3.0 4.0 5.0 6.0 6.4 7.0 8.0 8.5 9.0 9.5 10.0 10.5 11.0 12.0 13.0 14.0 15.0 16.0 18.0 19.0 20.0 22.0 24.0 26.0 28.0 30.0 35.0 40.0 45.0 50.0 55.0 60.0 70.0 80.0 
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

echo "${U}  ${UPRIME}   ${JHUND}"


for Nup in 5
do
mkdir -p Nup_${Nup}
cd Nup_${Nup}

for Ndn in 5  
do
mkdir -p Ndn_${Ndn}
cd Ndn_${Ndn}


rm input2orb_L${L}_Nup${Nup}_Ndn${Ndn}.inp
rm lanczos
cp ../../../../../input2orb_template.inp input2orb_L${L}_Nup${Nup}_Ndn${Ndn}.inp
cp ../../../../../lanczos .

input="input2orb_L${L}_Nup${Nup}_Ndn${Ndn}.inp"

##########################################
sed -i -e "s/NUP_VALUE/${Nup}/g" $input
sed -i -e "s/NDN_VALUE/${Ndn}/g" $input
sed -i -e "s/U_VALUE/${U}/g" $input
sed -i -e "s/JHUND_VALUE/${JHUND}/g" $input
sed -i -e "s/UPRIME_VALUE/${UPRIME}/g" $input
sed -i -e "s/ANISOTROPY_VALUE/${d_anis}/g" $input
sed -i -e "s/VALUE_L/${L}/g" $input
#########################################

rm "$input-e"
echo "submitted $L $U $d_anis $Nup $Ndn"

time ./lanczos $input > out.txt 



cd ..
done #lambda

cd ..
done #U


cd ..
done #N


cd ..
done #L

cd ..

done
