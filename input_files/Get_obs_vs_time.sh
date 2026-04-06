grep "Energy" out2.txt  > Energy_Lanczos.txt
grep "Psi(t)|H|Psi(t)" out2.txt > Energy.txt
grep "Splus local : 0" out2.txt > Splus.txt
grep "Sz local : 0" out2.txt > Sz.txt
grep "Sz2_Total" out2.txt > Sz2Total.txt

sed -i -e "s/ (/  /g" Energy.txt
sed -i -e "s/)/  /g" Energy.txt
sed -i -e "s/,/  /g" Energy.txt

sed -i -e "s/ (/  /g" Sz.txt
sed -i -e "s/)/  /g" Sz.txt
sed -i -e "s/,/  /g" Sz.txt

sed -i -e "s/ (/  /g" Splus.txt
sed -i -e "s/)/  /g" Splus.txt
sed -i -e "s/,/  /g" Splus.txt

sed -i -e "s/ (/  /g" Sz2Total.txt
sed -i -e "s/)/  /g" Sz2Total.txt
sed -i -e "s/,/  /g" Sz2Total.txt
