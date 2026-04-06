
EPS=40.0

rm file1.txt

for alpha in 0 1
do 
for i in {0..7}
do 
for beta in 0 1
do 
for j in {0..7}
do 
printf "u${i}${j}_${alpha}${beta} " >> file1.txt 
done
done
echo "" >> file1.txt
done
done

#Same Unitcell interactions------############
#-----------------------------------------

#Same site U01
val=$(awk 'NR==6 {printf $5}' Onsite_interactions.txt)
val2=${val%,*}
val3=${val2#*\(}
val3=$(echo "${val3}*(1.0/${EPS})" | bc -l)
echo "here 1 : ${val3} "
for site in {0..7}
do
sed -i "s/u${site}${site}_01/${val3}/" file1.txt
done


#Nearest neighbour density-density
Connections_site_col=(1 3 5 7 1 3 5 7 1 3 5 7 0 2 4 6)
Connections_site_row=(0 2 4 6 0 2 4 6 0 2 4 6 1 3 5 7)
Connections_orb_col=(0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1)
Connections_orb_row=(0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0)
Connections_orb_in_col=(2 2 2 2 3 3 3 3 3 3 3 3 2 2 2 2)
Connections_orb_in_row=(0 0 0 0 1 1 1 1 0 0 0 0 1 1 1 1)

for index in {0..15}
do
col=${Connections_site_col[${index}]}
row=${Connections_site_row[${index}]}
orb_col=${Connections_orb_col[${index}]}
orb_row=${Connections_orb_row[${index}]}
b=${Connections_orb_in_col[${index}]}
a=${Connections_orb_in_row[${index}]}

val=$(grep "${a}  ${b}  ${a}  ${b}" Onsite_interactions.txt | awk '{print $5}')
val2=${val%,*}
val3=${val2#*\(}
val3=$(echo "${val3}*(1.0/${EPS})" | bc -l)

sed -i "s/u${row}${col}_${orb_row}${orb_col}/${val3}/" file1.txt

done 


######################################################
#######Same Unit cell done ###########################



########### Neighbour pa1 ma2 ###################
##################################################
Connections_site_col=(2 6 2 6 2 6 1 5 3 7 3 7 3 7 0 4)
Connections_site_row=(1 5 1 5 1 5 2 6 0 4 0 4 0 4 3 7)
Connections_orb_col=(0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0)
Connections_subl_col=(0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1)
Connections_subl_row=(1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0)

for index in {0..15}
do

col=${Connections_site_col[${index}]}
row=${Connections_site_row[${index}]}
orb_col=${Connections_orb_col[${index}]}
orb_row=${Connections_orb_row[${index}]}
subl_col=${Connections_subl_col[${index}]}
subl_row=${Connections_subl_row[${index}]}

val=$(grep "${subl_col} ${subl_col} ${orb_col} ${subl_row} ${subl_row} ${orb_row} ${subl_col} ${subl_col} ${orb_col} ${subl_row} ${subl_row} ${orb_row}" NearestNeighbour_pa1_ma2_interactions.txt | awk '{print $13}')
val2=${val%,*}
val3=${val2#*\(}
val3=$(echo "${val3}*(1.0/${EPS})" | bc -l)

sed -i "s/u${row}${col}_${orb_row}${orb_col}/${val3}/" file1.txt

done

#####################################################
########### Neighbour pa1 ma2 done ################





########### Neighbour ma2 ###################
##################################################
Connections_site_col=(4 6 4 6 4 6 1 3 7 5 7 5 7 5 2 0)
Connections_site_row=(1 3 1 3 1 3 4 6 2 0 2 0 2 0 7 5)
Connections_orb_col=(0 0 1 1 1 1 1 1 0 0 1 1 1 1 1 1)
Connections_orb_row=(0 0 1 1 0 0 0 0 0 0 1 1 0 0 0 0)
Connections_subl_col=(0 0 0 0 0 0 1 1 0 0 0 0 0 0 1 1)
Connections_subl_row=(1 1 1 1 1 1 0 0 1 1 1 1 1 1 0 0)

for index in {0..15}
do

col=${Connections_site_col[${index}]}
row=${Connections_site_row[${index}]}
orb_col=${Connections_orb_col[${index}]}
orb_row=${Connections_orb_row[${index}]}
subl_col=${Connections_subl_col[${index}]}
subl_row=${Connections_subl_row[${index}]}

val=$(grep "${subl_col} ${subl_col} ${orb_col} ${subl_row} ${subl_row} ${orb_row} ${subl_col} ${subl_col} ${orb_col} ${subl_row} ${subl_row} ${orb_row}" NearestNeighbour_ma2_interactions.txt | awk '{print $13}')
val2=${val%,*}
val3=${val2#*\(}
val3=$(echo "${val3}*(1.0/${EPS})" | bc -l)

sed -i "s/u${row}${col}_${orb_row}${orb_col}/${val3}/" file1.txt

done

#####################################################
########### Neighbour ma2 done ################


#Removing all extra hoppings
for alpha in 0 1
do
for i in {0..7}
do
for beta in 0 1
do
for j in {0..7}
do
#sed -i "s/u${i}${j}_${alpha}${beta}/0/" Hopp_file.txt
sed -i "s/u${i}${j}_${alpha}${beta}/0/" file1.txt
done
done
done
done



cp file1.txt DenDenInt_file.txt
