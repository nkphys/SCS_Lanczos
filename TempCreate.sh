rm file1.txt
for alpha in a b
do
for i in 0 1 2 3
do
for beta in a b
do
for j in 0 1 2 3
do
printf "t${i}${j}_${alpha}${beta} " >> file1.txt
done
done
echo "" >> file1.txt
done
done

for str in t01_aa t23_aa t02_aa t13_aa
do
sed -i -e "s/${str}/-0.115/g" file1.txt
echo ""
done

for str in t01_bb t23_bb t02_bb t13_bb
do
sed -i -e "s/${str}/-0.492/g" file1.txt
echo ""
done


for str in t01_ab t23_ab t02_ab t13_ab t10_ab t32_ab t20_ab t31_ab
do
sed -i -e "s/${str}/0.24/g" file1.txt
#echo ""
done



