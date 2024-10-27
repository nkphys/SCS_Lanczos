Lx=6
Ly=4

for iy in {0..3}
do
for ix in {0..5}
do


ix_new=$(echo "${ix}+1" | bc -l)	
iy_new=$(echo "${iy}" | bc -l)

if [ ${ix_new} -eq ${Lx} ]
then
ix_new=0
fi

i=$(echo "${ix} + ${Lx}*${iy}" | bc -l)
i_new=$(echo "${ix_new} + ${Lx}*${iy_new}" | bc -l)

echo "${i} ${i_new}"

done
done



