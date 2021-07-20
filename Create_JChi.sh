J_Chi=$(echo "2*0.5" | bc -l)
file="ThreePointIntrs.txt"
echo "1" > ${file}
echo "SzSpSm 48" >> ${file}
#site1_array=(0 7 7 8 8 15 15 0 1 6 6 9 9 14 14 1 2 5 5 10 10 13 13 2 3 4 4 11 11 12 12 3)
#site2_array=(1 1 6 6 9 9 14 14 2 2 5 5 10 10 13 13 3 3 4 4 11 11 12 12 0 0 7 7 8 8 15 15)
#site3_array=(7 6 8 9 15 14 0 1 6 5 9 10 14 13 1 2 5 4 10 11 13 12 2 3 4 7 11 8 12 15 3 0)
site1_array=(0 2 2 4 4 6 6 0)
site2_array=(1 1 3 3 5 5 7 7)
site3_array=(2 3 4 5 6 7 0 1)



for i in {0..7}
do
site1_=${site1_array[${i}]}
site2_=${site2_array[${i}]}
site3_=${site3_array[${i}]}

echo "${site1_} ${site2_} ${site3_} (0,${J_Chi})" >> ${file}
echo "${site1_} ${site3_} ${site2_} (0,-${J_Chi})" >> ${file}
echo "${site3_} ${site2_} ${site1_} (0,-${J_Chi})" >> ${file}
echo "${site2_} ${site3_} ${site1_} (0,${J_Chi})" >> ${file}
echo "${site3_} ${site1_} ${site2_} (0,${J_Chi})" >> ${file}
echo "${site2_} ${site1_} ${site3_} (0,-${J_Chi})" >> ${file}


done

file="ThreePointOprs.txt"
echo "1" > ${file}
echo "SzSpSm 48" >> ${file}
site1_array=(0 2 2 4 4 6 6 0)
site2_array=(1 1 3 3 5 5 7 7)
site3_array=(2 3 4 5 6 7 0 1)
for i in {0..7}
do
site1_=${site1_array[${i}]}
site2_=${site2_array[${i}]}
site3_=${site3_array[${i}]}

echo "${site1_} ${site2_} ${site3_}" >> ${file}
echo "${site1_} ${site3_} ${site2_}" >> ${file}
echo "${site3_} ${site2_} ${site1_}" >> ${file}
echo "${site2_} ${site3_} ${site1_}" >> ${file}
echo "${site3_} ${site1_} ${site2_}" >> ${file}
echo "${site2_} ${site1_} ${site3_}" >> ${file}


done
