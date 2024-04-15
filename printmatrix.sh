file="CdCdCC_measurement_LanczosState0.txt"
file1="pp.dat"
row=1
rm $file1
for i in $(seq 0 1 7)
do
for j in $(seq 0 1 7)
do
	((row++))
	awk -v r="$row" 'FNR == r {printf ("%.10f", $13)}' $file >>$file1
printf "\t" >>$file1
done
printf "\n" >>$file1
done
