Length=20
TwoTimesS=1

mkdir -p L${Length}_OBC
cd L${Length}_OBC

mkdir -p TwoTimesS${TwoTimesS}
cd TwoTimesS${TwoTimesS}

for n in {0..10} #k=2*pi*n/L
do

mkdir -p n${n}
cd n${n}


string_opr="${Length} Sz"
for site in {0..19}
do
k=$(echo "(2.0*3.14159265358*${n})/(${Length})" | bc -l)
real=$(echo "(c(${k}*${site}))/(sqrt(${Length}))" | bc -l)
imag=$(echo "(s(${k}*${site}))/(sqrt(${Length}))" | bc -l)
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
