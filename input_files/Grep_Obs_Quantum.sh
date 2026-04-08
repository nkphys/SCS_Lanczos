file="QAH_quantum.txt"
rm $file
echo "#U total local" > $file
for U in 0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 3.95 4.00 4.05 4.10 4.15 4.20 4.25 4.30 4.35 4.40 4.45 4.50 4.55 4.60 4.65 4.70 4.75 4.80 4.85 4.90 4.95 5.00 5.05 5.10 5.15 5.20 5.25 5.30 5.35 5.40 5.45 5.50 5.55 5.60 5.65 5.70 5.75 5.80 5.85 5.90 5.95 6.00 6.05 6.10 6.15 6.20 6.25 6.30 6.35 6.40 6.45 6.50 6.55 6.60 6.65 6.70 6.75 6.80 6.85 6.90 6.95 7.00 7.05 7.10 7.15 7.20 7.25 7.30 7.35 7.40 7.45 7.50 7.55 7.60 7.65 7.70 7.75 7.80 7.85 7.90 7.95 8.00
	#0.20 0.40 0.60 0.80 1.00 1.20 1.40 1.60 1.80 2.00 2.20 2.40 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20 5.40 5.60 5.80 6.00 6.20 6.40 6.60 6.80 7.00 7.20 7.40 7.60 7.80 8.00
	#0.00 0.05 0.10 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40 1.45 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.85 1.90 1.95 2.00 2.05 2.10 2.15 2.20 2.25 2.30 2.35 2.40 2.45 2.50 2.55 2.60 2.65 2.70 2.75 2.80 2.85 2.90 2.95 3.00 3.05 3.10 3.15 3.20 3.25 3.30 3.35 3.40 3.45 3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 3.95 4.00 4.05 4.10 4.15 4.20 4.25 4.30 4.35 4.40 4.45 4.50 4.55 4.60 4.65 4.70 4.75 4.80 4.85 4.90 4.95 5.00 5.05 5.10 5.15 5.20 5.25 5.30 5.35 5.40 5.45 5.50 5.55 5.60 5.65 5.70 5.75 5.80 5.85 5.90 5.95 6.00 6.05 6.10 6.15 6.20 6.25 6.30 6.35 6.40 6.45 6.50 6.55 6.60 6.65 6.70 6.75 6.80 6.85 6.90 6.95 7.00 7.05 7.10 7.15 7.20 7.25 7.30 7.35 7.40 7.45 7.50 7.55 7.60 7.65 7.70 7.75 7.80 7.85 7.90 7.95 8.00
#	0.00 0.20 0.40 0.60 0.80 1.00 1.20 1.40 1.60 1.80 2.00 2.20 2.40 2.60 2.80 3.00 3.20 3.40 3.60 3.80 4.00 4.20 4.40 4.60 4.80 5.00 5.20 5.40 5.60 5.80 6.00 6.20 6.40 6.60 6.80 7.00 7.20 7.40 7.60 7.80 8.00
do
totalval=$(grep "Total Value Sum Quantum" U_${U}/V1_1.0/Run_out.txt | awk '{print $6}')

nonlocalval0=$(grep "Total quantum for set 0" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval1=$(grep "Total quantum for set 1 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval2=$(grep "Total quantum for set 2 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval3=$(grep "Total quantum for set 3" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval4=$(grep "Total quantum for set 4" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval5=$(grep "Total quantum for set 5" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval6=$(grep "Total quantum for set 6" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval7=$(grep "Total quantum for set 7" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval8=$(grep "Total quantum for set 8" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval9=$(grep "Total quantum for set 9" U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval10=$(grep "Total quantum for set 10 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval11=$(grep "Total quantum for set 11 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval12=$(grep "Total quantum for set 12 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval13=$(grep "Total quantum for set 13 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval14=$(grep "Total quantum for set 14 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')
nonlocalval15=$(grep "Total quantum for set 15 " U_${U}/V1_1.0/Run_out.txt | awk '{print $7}')

echo "${U} ${totalval} ${nonlocalval0} ${nonlocalval1}  ${nonlocalval2} ${nonlocalval3} ${nonlocalval4}  ${nonlocalval5}  ${nonlocalval6}  ${nonlocalval7}  ${nonlocalval8} ${nonlocalval9}  ${nonlocalval10}  ${nonlocalval11}  ${nonlocalval12}  ${nonlocalval13}  ${nonlocalval14}  ${nonlocalval15}" >> $file
done

sed -i -e "s/(//g" $file
sed -i -e "s/)//g" $file
sed -i -e "s/,/  /g" $file
