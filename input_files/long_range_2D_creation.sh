Lx=6
Ly=4

for iy in {0..3}
do
for ix in {0..5}
do


for jy in {0..3}
do
for jx in {0..5}
do

val=0.0

dis_y=$(echo "${jy}-${iy}" | bc -l)
dis_x=$(echo "${jx}-${ix}" | bc -l)
if [ ${dis_y} -lt 0 ]
then
dis_y=$(echo "-1*${dis_y}" | bc -l)
fi 
if [ ${dis_x} -lt 0 ]
then
dis_x=$(echo "-1*${dis_x}" | bc -l)
fi

if [ ${dis_y} -eq 1 ] || [ ${dis_y} -eq 3 ]
then
dis_y=1
fi
if [ ${dis_x} -eq 1 ] || [ ${dis_x} -eq 5 ]
then
dis_x=1
fi

dis=$(echo "${dis_y}+${dis_x}" | bc -l)

if [ ${dis} -eq 1 ]
then
val="tval"
i=$(echo "${ix} + ${iy}*${Lx}" | bc -l)
j=$(echo "${jx} + ${jy}*${Lx}" | bc -l)

if [ ${j} -gt ${i} ]
then
echo "${i}  ${j}  ${val}"
fi

fi

i=$(echo "${ix} + ${iy}*${Lx}" | bc -l)
j=$(echo "${jx} + ${jy}*${Lx}" | bc -l)

if [ ${i} -gt ${j} ]
then
val=0.0
fi

#printf "${val} "

done
done


#echo ""
done
done
