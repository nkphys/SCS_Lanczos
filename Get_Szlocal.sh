
for site in 0 1 2 3
do

grep "Sz local \[Psi(t)\] : ${site}" out.txt > Szlocal${site}.txt
grep "Splus local \[Psi(t)\] : ${site}" out.txt > Spluslocal${site}.txt


sed -i -e "s/,/ /g" Szlocal${site}.txt
sed -i -e "s/(/ /g" Szlocal${site}.txt
sed -i -e "s/)/ /g" Szlocal${site}.txt

sed -i -e "s/,/ /g" Spluslocal${site}.txt
sed -i -e "s/(/ /g" Spluslocal${site}.txt
sed -i -e "s/)/ /g" Spluslocal${site}.txt

done


grep "Entropy \[H(t)\]" out.txt > Entropy_Ht.txt
grep "Entropy \[Psi(t)\]" out.txt > Entropy_Psit.txt

grep " |Psi(t)> :" out.txt > WF_Psit.txt
sed -i -e "s/,/ /g" WF_Psit.txt
sed -i -e "s/(/ /g" WF_Psit.txt
sed -i -e "s/)/ /g" WF_Psit.txt


grep "CHSH_std_01" out.txt  > CHSH_std_01.txt
sed -i -e "s/,/ /g" CHSH_std_01.txt
sed -i -e "s/(/ /g" CHSH_std_01.txt
sed -i -e "s/)/ /g" CHSH_std_01.txt

grep "CHSH_std_12" out.txt  > CHSH_std_12.txt
sed -i -e "s/,/ /g" CHSH_std_12.txt
sed -i -e "s/(/ /g" CHSH_std_12.txt
sed -i -e "s/)/ /g" CHSH_std_12.txt


