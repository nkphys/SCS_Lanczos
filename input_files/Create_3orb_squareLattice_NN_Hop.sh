
rm file1.txt

for alpha in 0 1 2
do 
for i in {0..3}
do 
for beta in 0 1 2
do 
for j in {0..3}
do 
printf "t${i}${j}_${alpha}${beta} " >> file1.txt 
done
done
echo "" >> file1.txt
done
done


#Nearest neighbour intra-orbital
Connections_site_col=(1 3 1 3)
Connections_site_row=(0 2 0 2)
Connections_orb_col=(0 0 2 2)
Connections_orb_row=(0 0 2 2)


for index in {0..3}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
beta=${Connections_orb_col[${index}]}
alpha=${Connections_orb_row[${index}]}
sed -i "s/t${i}${j}_${alpha}${beta}/-1.0/" file1.txt
done

Connections_site_col=(2 3 2 3)
Connections_site_row=(0 1 0 1)
Connections_orb_col=(1 1 2 2)
Connections_orb_row=(1 1 2 2)
for index in {0..3}
do
j=${Connections_site_col[${index}]}
i=${Connections_site_row[${index}]}
beta=${Connections_orb_col[${index}]}
alpha=${Connections_orb_row[${index}]}
sed -i "s/t${i}${j}_${alpha}${beta}/-1.0/" file1.txt
done



cp file1.txt Hopp_file.txt


#Removing all extra hoppings
for alpha in 0 1 2
do
for i in {0..7}
do
for beta in 0 1 2
do
for j in {0..7}
do
sed -i "s/t${i}${j}_${alpha}${beta}/0/" Hopp_file.txt
sed -i "s/t${i}${j}_${alpha}${beta}/0/" file1.txt
done
done
done
done
