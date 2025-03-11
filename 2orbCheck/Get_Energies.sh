for U in 50.0
do

for Jh_by_U in 0.000 0.020 0.050 0.080 0.100 0.120 0.150 0.180 0.200 0.210 0.215 0.220 0.225 0.230 0.235 0.240 0.242 0.244 0.246 0.248 0.250 0.252 0.254 0.256 0.258 0.260 0.265 0.270 0.275 0.280 0.285 0.290 0.295 0.300 0.305 0.310 0.315 0.320 
do


E0=$(grep "^0 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E1=$(grep "^1 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E2=$(grep "^2 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E3=$(grep "^3 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E4=$(grep "^4 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E5=$(grep "^5 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E6=$(grep "^6 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E7=$(grep "^7 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')
E8=$(grep "^8 " U_${U}/Jh_by_U_${Jh_by_U}/out_run.txt | awk '{print $2}')



echo "${U}   ${Jh_by_U}  ${E0}  ${E1}  ${E2}  ${E3}  ${E4}  ${E5}  ${E6}  ${E7}  ${E8}"
done

done
