#!/bin/bash
for i in {1..16}
do 
for j in {1..16}
do 
temp_i=$i
temp_j=$j
if [ $i -gt $j ]
then
temp_i=$j
temp_j=$i 
fi
val=$(awk -v row=$temp_i -v col=$temp_j 'NR==row {print $col}' Hoppings_4x4_SD.txt)
printf "(-${val},0.0) "
done
echo ""
done
