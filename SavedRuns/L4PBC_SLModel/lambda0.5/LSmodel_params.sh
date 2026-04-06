Length=4
t=1.0
U=40.0
J=$(echo "${U}*0.25" | bc -l)
lambda=0.5

Term_offset=$(echo "(-8.0/3.0)*((${t}*${t})/(${U}-3.0*${J})) +  ((-4.0/3.0)*((${t}*${t})/${U})) + 0.0*((${lambda}*${lambda})/(20.0*${J}))" | bc -l)
Term_LL=$( echo "( (4.0/3.0)*( (${t}*${t})/(${U}-3*${J}) ) )  +  (  (-1.0/3.0)*( (${t}*${t})/(${U}) )  )"  | bc -l )
Term_LL2=$( echo "( (4.0/3.0)*( (${t}*${t})/(${U}-3*${J}) ) )  +  (  (1.0/6.0)*( (${t}*${t})/(${U}) )  )   +   ( (-1.0*${t}*${t})/(2.0*(${U}+2.0*${J}))  )"  | bc -l )
Term_SSLL=$( echo "( (2.0/3.0)*( (${t}*${t})/(${U}-3*${J}) ) )  +  (  (1.0/3.0)*( (${t}*${t})/(${U}) )  )"  | bc -l )
Term_SSLL2=$( echo "( (2.0/3.0)*( (${t}*${t})/(${U}-3*${J}) ) )  +  (  (-1.0/6.0)*( (${t}*${t})/(${U}) )  )  +  ( (1.0/2.0)*( (${t}*${t})/(${U}+2*${J}) ) )"  | bc -l )
Term_SS=$( echo "( (-4.0/3.0)*( (${t}*${t})/(${U}-3*${J}) ) )  +  (  (4.0/3.0)*( (${t}*${t})/(${U}) )  )"  | bc -l )

Term_SL=$( echo "(  (${lambda}/2.0)  - 0.0*((${lambda}*${lambda})/(8.0*${J}))   )" | bc -l  )


Term_offset=$(printf "%1.10f" ${Term_offset})
Term_LL=$(printf "%1.10f" ${Term_LL})
Term_LL2=$(printf "%1.10f" ${Term_LL2})
Term_SSLL=$(printf "%1.10f" ${Term_SSLL})
Term_SSLL2=$(printf "%1.10f" ${Term_SSLL2})
Term_SS=$(printf "%1.10f" ${Term_SS})
Term_SL=$(printf "%1.10f" ${Term_SL})

#echo "Term_offset=${Term_offset}"
#echo "Term_LL2=${Term_LL2}"
#echo "Term_SSLL=${Term_SSLL}"
#echo "Term_SSLL2=${Term_SSLL2}"
#echo "Term_SS=${Term_SS}"
#echo "Term_SL=${Term_SL}"
#echo "Term_LL=${Term_LL}"

sed -i "s/Term_SS/${Term_SS}/g" J1file_L${Length}.dat
sed -i "s/Term_LL/${Term_LL}/g" J1file_L${Length}.dat
sed -i "s/Term_SL/${Term_SL}/g" J1file_L${Length}.dat
sed -i "s/Term_LL2/${Term_LL2}/g" J2file_L${Length}.dat
sed -i "s/Term_SSLL/${Term_SSLL}/g" J2file_L${Length}.dat
sed -i "s/Term_SSLL2/${Term_SSLL2}/g" J3file_L${Length}.dat
