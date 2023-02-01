
rm file1.txt

for alpha in 0 1
do 
for i in {0..7}
do 
for beta in 0 1
do 
for j in {0..7}
do 
printf "t${i}${j}_${alpha}${beta} " >> file1.txt 
done
done
echo "" >> file1.txt
done
done


#Nearest neighbour [Blue hoppings] intra-orbital
Connections_site_col=(1 2 4 3 5 6 6 7)
Connections_site_row=(0 1 1 2 4 3 5 6)
for beta in 0 1
do
for index in {0..7}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
sed -i "s/t${i}${j}_${beta}${beta}/T${i}${j}_${beta}${beta}/" file1.txt
done
done

#Nearest neighbour [Blue hoppings] inter-orbital [1-->0]
Connections_site_col=(1 2 4 3 5 6 6 7 0 1 1 2 4 3 5 6)
Connections_site_row=(0 1 1 2 4 3 5 6 1 2 4 3 5 6 6 7)

for index in {0..14}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
sed -i "s/t${i}${j}_01/T${i}${j}_01/" file1.txt
done




#Second Nearest neighbour [Green hoppings] intra-orbital
Connections_site_col=(4 5 6 7 2 3 6 7 4 5) 
Connections_site_row=(0 1 2 3 0 1 4 5 2 3)
for beta in 0 1
do
for index in {0..9}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
#sed -i "s/t${i}${j}_${beta}${beta}/T${i}${j}_${beta}${beta}/" file1.txt
done
done

#Second Nearest neighbour [Green hoppings] inter-orbital
Connections_site_col=(4 5 6 7 2 3 6 7 4 5 0 1 2 3 0 1 4 5 2 3)
Connections_site_row=(0 1 2 3 0 1 4 5 2 3 4 5 6 7 2 3 6 7 4 5)
for index in {0..19}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
#sed -i "s/t${i}${j}_01/T${i}${j}_01/" file1.txt
done





#Third Nearest neighbour [Red hoppings] intra-orbital
Connections_site_col=(6 4 5)
Connections_site_row=(1 3 2)
for beta in 0 1
do
for index in {0..2}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
#sed -i "s/t${i}${j}_${beta}${beta}/T${i}${j}_${beta}${beta}/" file1.txt
done
done

#Third Nearest neighbour [Red hoppings] inter-orbital
Connections_site_col=(6 4 5 1 3 2)
Connections_site_row=(1 3 2 6 4 5)
for index in {0..5}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
#sed -i "s/t${i}${j}_01/T${i}${j}_01/" file1.txt
done



cp file1.txt Hopp_file.txt


################################################################
#Replacing with actual values---------------
################################################################


#Bonds 01, 23, 45, 67
#T01_00   ( Tij_ab = (j,b)-->(i,a), (j,b)=j+2bi)
#NOTE: T01_10=T_10_01  for real hoppings
j_array=(1 1 1 1)
b_array=(0 1 1 0)
i_array=(0 0 0 0)
a_array=(0 1 0 1)
val_array=(X X X X)

for index in {0..3}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t0_mat_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}

val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done
#for m>l
#val_l0_m0=${val_array[0]}
#val_l1_m1=${val_array[1]}
#val_l0_m1=${val_array[2]}
#val_m0_l1=${val_array[3]}

Connections_site_col=(1 3 5 7 1 3 5 7 1 3 5 7 0 2 4 6)
Connections_site_row=(0 2 4 6 0 2 4 6 0 2 4 6 1 3 5 7)
Connections_orb_col=(0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1)
Connections_orb_row=(0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0)
Connections_value=(${val_array[0]} ${val_array[0]} ${val_array[0]} ${val_array[0]} ${val_array[1]} ${val_array[1]} ${val_array[1]} ${val_array[1]} ${val_array[2]} ${val_array[2]} ${val_array[2]} ${val_array[2]} ${val_array[3]} ${val_array[3]} ${val_array[3]} ${val_array[3]}) 
for index in {0..15}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done



#Bonds 14, 36
j_array=(0 0 0 0)
b_array=(0 1 1 0)
i_array=(1 1 1 1)
a_array=(0 1 0 1)
val_array=(X X X X)

for index in {0..3}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_minus_a2_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done

Connections_site_col=(4 6 4 6 4 6 1 3)
Connections_site_row=(1 3 1 3 1 3 4 6)
Connections_orb_col=(0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0)
Connections_value=(${val_array[0]} ${val_array[0]} ${val_array[1]} ${val_array[1]} ${val_array[2]} ${val_array[2]} ${val_array[3]} ${val_array[3]})
for index in {0..7}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done



#Bonds 12, 56
j_array=(0 0 0 0)
b_array=(0 1 1 0)
i_array=(1 1 1 1)
a_array=(0 1 0 1)
val_array=(X X X X)

for index in {0..3}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_plus_a1_minus_a2_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}

val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done

Connections_site_col=(2 6 2 6 2 6 1 5)
Connections_site_row=(1 5 1 5 1 5 2 6)
Connections_orb_col=(0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0)
Connections_value=(${val_array[0]} ${val_array[0]} ${val_array[1]} ${val_array[1]} ${val_array[2]} ${val_array[2]} ${val_array[3]} ${val_array[3]})

for index in {0..7}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done



#Bonds 04, 15, 26, 37
j_array=(0 0 0 0 1 1 1 1)
b_array=(0 1 1 0 0 1 1 0)
i_array=(0 0 0 0 1 1 1 1)
a_array=(0 1 0 1 0 1 0 1)
val_array=(X X X X X X X X)

for index in {0..7}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_minus_a2_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done

Connections_site_col=(4 6 4 6 4 6 0 2 5 7 5 7 5 7 1 3)
Connections_site_row=(0 2 0 2 0 2 4 6 1 3 1 3 1 3 5 7)
Connections_orb_col=(0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0)
Connections_value=(${val_array[0]} ${val_array[0]} ${val_array[1]} ${val_array[1]} ${val_array[2]} ${val_array[2]} ${val_array[3]} ${val_array[3]} ${val_array[4]} ${val_array[4]} ${val_array[5]} ${val_array[5]} ${val_array[6]} ${val_array[6]} ${val_array[7]} ${val_array[7]})

for index in {0..15}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done




#Bonds 02=46, 13=57
j_array=(0 0 0 0 1 1 1 1)
b_array=(0 1 1 0 0 1 1 0)
i_array=(0 0 0 0 1 1 1 1)
a_array=(0 1 0 1 0 1 0 1)
val_array=(X X X X X X X X)

for index in {0..7}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_plus_a1_minus_a2_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done



Connections_site_col=(2 6 2 6 2 6 0 4 7 3 7 3 7 3 5 1)
Connections_site_row=(0 4 0 4 0 4 2 6 5 1 5 1 5 1 7 3)
Connections_orb_col=(0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0)
Connections_value=(${val_array[0]} ${val_array[0]} ${val_array[1]} ${val_array[1]} ${val_array[2]} ${val_array[2]} ${val_array[3]} ${val_array[3]} ${val_array[4]} ${val_array[4]} ${val_array[5]} ${val_array[5]} ${val_array[6]} ${val_array[6]} ${val_array[7]} ${val_array[7]})

for index in {0..15}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done





#Bonds 24, 35
j_array=(0 0 0 0 1 1 1 1)
b_array=(0 1 1 0 0 1 1 0)
i_array=(0 0 0 0 1 1 1 1)
a_array=(0 1 0 1 0 1 0 1)
val_array=(X X X X X X X X)

for index in {0..7}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_plus_a1_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done

Connections_site_col=(4 4 4 2 5 5 5 3)
Connections_site_row=(2 2 2 4 3 3 3 5)
Connections_orb_col=(0 1 1 1 0 1 1 1)
Connections_orb_row=(0 1 0 0 0 1 0 0)
# The order is changed because of the direction on hopping in t1_plus_a1_hop_(3rd)NN.txt file
Connections_value=(${val_array[0]} ${val_array[1]} ${val_array[3]} ${val_array[2]} ${val_array[4]} ${val_array[5]} ${val_array[7]} ${val_array[6]})
for index in {0..7}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done






#Bonds 34, 25
j_array=(0 0 0 0 1 1 1 1)
b_array=(0 1 1 0 0 1 1 0)
i_array=(1 1 1 1 0 0 0 0)
a_array=(0 1 0 1 0 1 0 1)
val_array=(X X X X X X X X)

for index in {0..7}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_plus_a1_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done



Connections_site_col=(4 4 4 3 5 5 5 2)
Connections_site_row=(3 3 3 4 2 2 2 5)
Connections_orb_col=(0 1 1 1 0 1 1 1)
Connections_orb_row=(0 1 0 0 0 1 0 0)
# The order is changed because of the direction on hopping in t1_plus_a1_hop_(3rd)NN.txt file
Connections_value=(${val_array[4]} ${val_array[5]} ${val_array[7]} ${val_array[6]} ${val_array[0]} ${val_array[1]} ${val_array[3]} ${val_array[2]})
for index in {0..7}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done



#Bond 16
j_array=(0 0 0 0)
b_array=(0 1 1 0)
i_array=(1 1 1 1)
a_array=(0 1 0 1)
val_array=(X X X X)

for index in {0..3}
do
j=${j_array[${index}]}
b=${b_array[${index}]}
i=${i_array[${index}]}
a=${a_array[${index}]}

row=$(echo "${i}+2*${a}+1" | bc -l)
col=$(echo "${j}+2*${b}+1" | bc -l)
val=$(awk -v c=${col} -v r=${row} 'NR==r {printf $c}' t1_plus_a1_minus_2a2_hop_NN.txt)
val2=${val%,*}
val3=${val2#*\(}
val_array[${index}]=${val3}
#echo "${row}  ${col} ${val_array[${index}]}"
done


Connections_site_col=(6 6 6 1)
Connections_site_row=(1 1 1 6)
Connections_orb_col=(0 1 1 1)
Connections_orb_row=(0 1 0 0)
Connections_value=(${val_array[0]} ${val_array[1]} ${val_array[2]} ${val_array[3]})
for index in {0..3}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
oj=${Connections_orb_col[${index}]}
oi=${Connections_orb_row[${index}]}
value=${Connections_value[${index}]}
sed -i "s/T${i}${j}_${oi}${oj}/${value}/" Hopp_file.txt
done




###################################################################









#Removing all extra hoppings
for alpha in 0 1
do
for i in {0..7}
do
for beta in 0 1
do
for j in {0..7}
do
sed -i "s/t${i}${j}_${alpha}${beta}/0/" Hopp_file.txt
sed -i "s/t${i}${j}_${alpha}${beta}/0/" file1.txt
done
done
done
done
