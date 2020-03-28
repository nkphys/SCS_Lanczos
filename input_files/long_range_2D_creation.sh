Lx=2
Ly=2
for iy in {0..1}
do
for ix in {0..1}
do


for jy in {0..1}
do
for jx in {0..1}
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
if [ ${dis_x} -eq 1 ] || [ ${dis_x} -eq 3 ]
then
dis_x=1
fi

dis=$(echo "${dis_y}+${dis_x}" | bc -l)

if [ ${dis} -eq 1 ]
then
val="tval"
fi

printf "${val} "
done
done


echo ""
done
done
