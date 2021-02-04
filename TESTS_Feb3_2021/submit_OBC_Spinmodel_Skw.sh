Length=4
TwoTimesS=2

mkdir -p L${Length}_OBC
cd L${Length}_OBC

mkdir -p TwoTimesS${TwoTimesS}
cd TwoTimesS${TwoTimesS}

for n in {1..4} #k=pi*n/L+1
do
mkdir -p n${n}_reortho
cd n${n}_reortho


string_opr="${Length} Sz"
for site in {0..7}
do
site_OBC=$(echo "${site}+1" | bc -l)
k=$(echo "(3.14159265358*${n})/(${Length}+1.0)" | bc -l)
real=$(echo "(s(${k}*${site_OBC}))*(sqrt(2.0/(${Length}+1.0)))" | bc -l)
imag=$(echo "0.0" | bc -l)
real=$(printf "%1.5f" ${real})
imag=$(printf "%1.5f" ${imag})

string_opr="${string_opr} (${real},${imag}) ${site}"
done

input="input_run.inp"

rm ${input}
rm lanczos

cp ../../../lanczos .
cp ../../../input_SpinOnlyTargetSz_template.inp ${input}

sed -i -e "s/LENGTH_VALUE/${Length}/g" ${input}
sed -i -e "s/TWOTIMESSPIN_VALUE/${TwoTimesS}/g" ${input}
sed -i -e "s/DYN_OPR_STRING_VALUE/${string_opr}/g" ${input}

time ./lanczos ${input} > out_run.txt

cd ..
done


cd ..
cd ..
