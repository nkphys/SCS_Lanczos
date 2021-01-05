#Sz (1.0,0.0) 0

L=8
n=4 

printf "Sz "
for site in {0..7}
do
k=$(echo "(2.0*3.14159265358*${n})/(${L})" | bc -l)
real=$(echo "(c(${k}*${site}))/(sqrt(${L}))" | bc -l)
imag=$(echo "(s(${k}*${site}))/(sqrt(${L}))" | bc -l)
real=$(printf "%1.5f" ${real})
imag=$(printf "%1.5f" ${imag})

printf "(${real},${imag}) ${site} "
done

echo ""
