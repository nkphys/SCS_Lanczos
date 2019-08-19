#!/bin/bash -x
Length=$1
Jxzz_="Jxzz_"
Jyzz_="Jyzz_"
Jxpm_="Jxpm_"
Jypm_="Jypm_"
tx_="tx_"
ty_="ty_"
Ux_="Ux_"
Uy_="Uy_"



file_hopping=hopping_longrange_$Length.txt
file_Jzz=Jzz_longrange_$Length.txt
file_Jpm=Jpm_longrange_$Length.txt
file_denden=denden_longrange_$Length.txt

printf "Input params beaing created for L=$Length system \n"
rm $file_hopping;
rm $file_Jzz;
rm $file_Jpm;
rm $file_denden;

term1=('0.0');
term2=('0.0');
term3=('0.0');
term4=('0.0');
counter=1
while [ $counter -lt $Length ];
do
term1=("${term1[@]}" "0.0")
term2=("${term2[@]}" "0.0")
term3=("${term3[@]}" "0.0")
term4=("${term4[@]}" "0.0")


let counter=counter+1
done
#printf ${#J_zz[@]}
#--------------RIGHT HERE THE CONNECTIONS VALUES like A[d];where d=|i-j|------
term1[1]="$ty_"
term1[2]="$tx_"
term2[1]="$Jyzz_"
term2[2]="$Jxzz_"
term3[1]="$Jypm_"
term3[2]="$Jxpm_"
term4[1]="$Uy_"
term4[2]="$Ux_"
#----------------CONNECTIONS VALUES WRITTEN-----------------------------------

#Jzz
site1=0
while [ $site1 -lt $Length ];
do
site2=0
while [ $site2 -lt $Length ];
do
if [ $site1 -lt $site2 ]; then
dis=$(echo "$site2 -$site1" | bc)
fi
if [ $site2 -lt $site1 ]; then
dis=$(echo "$site1 -$site2" | bc)
fi
if [ "$site2" -eq "$site1" ]; then
dis=0
fi

if [ $(echo "$site2 -$site1" | bc) -eq 1 ]; then

        if [ $((site1%2)) -eq 0 ];then
        printf  "(${term1[$dis]},0) "  >> $file_hopping
	printf "(${term2[$dis]},0) ">> $file_Jzz
	printf "(${term3[$dis]},0) ">> $file_Jpm
	printf "(${term4[$dis]},0) ">> $file_denden
        else
        printf  "(0,0) "  >> $file_hopping
        printf "(0,0) ">> $file_Jzz
        printf "(0,0) ">> $file_Jpm
        printf "(0,0) ">> $file_denden
        fi

elif [ $(echo "$site1 -$site2" | bc) -eq 1 ]; then

       if [ $((site1%2)) -eq 0 ];then
        printf  "(0,0) "  >> $file_hopping
        printf "(0,0) ">> $file_Jzz
        printf "(0,0) ">> $file_Jpm
        printf "(0,0) ">> $file_denden
        else
        printf  "(${term1[$dis]},0) "  >> $file_hopping
        printf "(${term2[$dis]},0) ">> $file_Jzz
        printf "(${term3[$dis]},0) ">> $file_Jpm
        printf "(${term4[$dis]},0) ">> $file_denden
       fi

else
        printf "(${term1[$dis]},0)  ">> $file_hopping
	printf "(${term2[$dis]},0) ">> $file_Jzz
	printf "(${term3[$dis]},0) ">> $file_Jpm
	printf "(${term4[$dis]},0) ">> $file_denden
fi

#printf "(${term1[$dis]},0)  ">> $file_hopping
#printf "(${term2[$dis]},0) ">> $file_Jzz
#printf "(${term3[$dis]},0) ">> $file_Jpm
#printf "(${term4[$dis]},0) ">> $file_denden

#printf "$dis "
let site2=site2+1
done
printf "\n" >> $file_hopping
printf "\n ">> $file_Jzz
printf "\n ">> $file_Jpm
printf "\n ">> $file_denden
let site1=site1+1
done
