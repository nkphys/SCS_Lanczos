val=-0.005
for i in {1..100}
do
val=$(echo "${val}+0.005"| bc -l)
printf "%1.4f " ${val}
done
echo ""
