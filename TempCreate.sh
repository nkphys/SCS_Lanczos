for alpha in 0 1
do
for i in 0 1
do
for beta in 0 1
do
for j in 0 1
do
printf "t${i}${j}_${alpha}${beta} " >> file1.txt
done
done
echo "" >> file1.txt
done
done

