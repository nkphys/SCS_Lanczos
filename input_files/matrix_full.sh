
for c in {1..28}
do

for r in {1..28}
do

if [ "$r" -gt "$c"  ]
then
val=$(awk -v row=${c} -v col=${r} 'NR==row {print $col}' temp_ypx.txt)
val1=$(awk -v row=${c} -v col=${r} 'NR==row {print $col}' temp_x.txt)
val2=$(awk -v row=${c} -v col=${r} 'NR==row {print $col}' temp_y.txt)
else
val=$(awk -v row=${r} -v col=${c} 'NR==row {print $col}' temp_ypx.txt)
val1=$(awk -v row=${r} -v col=${c} 'NR==row {print $col}' temp_x.txt)
val2=$(awk -v row=${r} -v col=${c} 'NR==row {print $col}' temp_y.txt)
fi

val_sum=$(echo "${val1}+${val2}" | bc -l)

printf " ${val1} "

done 
echo
done
