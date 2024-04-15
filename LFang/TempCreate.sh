input="pp1_input.txt"

rm pp1_input.txt
for i in $(seq 0 1 7)
do
for j in $(seq 0 1 7)
do
printf "$i 0 0 $i 0 1 $j 0 1 $j 0 0\n" >> $input
done
done
