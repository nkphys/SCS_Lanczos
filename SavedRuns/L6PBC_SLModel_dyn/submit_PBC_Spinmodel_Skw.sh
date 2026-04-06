eta=0.0001
Length=12
Lengthbytwo=6
TwoTimesS=2

mkdir -p L${Length}
cd L${Length}

for lambda in 0.42
do
mkdir -p lambda${lambda}
cd lambda${lambda}

for n in 0 #k=2*pi*n/L
do
mkdir -p n${n}_eta${eta}
cd n${n}_eta${eta}


string_opr="${Length} Sz"
for site in {0..5}
do
site_S=$(echo "2*${site}" | bc -l)
k=$(echo "(2.0*3.14159265358*${n})/(${Lengthbytwo})" | bc -l)
real=$(echo "2*(c(${k}*${site}))/(sqrt(${Lengthbytwo}))" | bc -l)
imag=$(echo "2*(s(${k}*${site}))/(sqrt(${Lengthbytwo}))" | bc -l)
real=$(printf "%1.5f" ${real})
imag=$(printf "%1.5f" ${imag})
string_opr="${string_opr} (${real},${imag}) ${site_S}"
done

for site in {0..5}
do
site_L=$(echo "2*${site}+1" | bc -l)
k=$(echo "(2.0*3.14159265358*${n})/(${Lengthbytwo})" | bc -l)
real=$(echo "(c(${k}*${site}))/(sqrt(${Lengthbytwo}))" | bc -l)
imag=$(echo "(s(${k}*${site}))/(sqrt(${Lengthbytwo}))" | bc -l)
real=$(printf "%1.5f" ${real})
imag=$(printf "%1.5f" ${imag})
string_opr="${string_opr} (${real},${imag}) ${site_L}"
done


input="input_run.inp"

rm ${input}
rm lanczos

cp ../../../{LSmodel_params.sh,J1file_L6.dat,J2file_L6.dat,J3file_L6.dat} .
cp ../../../lanczos .
cp ../../../input_SpinOnlyTargetSz.inp ${input}

sed -i -e "s/LAMBDA_VALUE/${lambda}/g" LSmodel_params.sh
bash LSmodel_params.sh

sed -i -e "s/DYN_OPR_STRING_VALUE/${string_opr}/g" ${input}
sed -i -e "s/VALUE_ETA/${eta}/g" ${input}

time ./lanczos ${input} > out_run.txt

cd ..
done

cd ..
done

cd ..
