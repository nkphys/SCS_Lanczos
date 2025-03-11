t=1.0
U=VALUE_U
fac_=-1.0

JHbyU=VALUE_JH
#0.00 0.02 0.05 0.08 0.10 0.12 0.15 0.17 0.20 0.22 0.25 0.27 0.30 0.32


JH_temp=$(echo "${U}*${JHbyU}" | bc -l)
JH=$(printf "%1.5f" ${JH_temp})

file="SS_LL_X.txt"
cp SS_LL_X_template.txt ${file}

SQ1_=$(echo "((${t}*${t})/(${U}))*(  ( (${U}+${JH})/(${U}+2.0*${JH})  )  + ( (${JH})/(${U}-3.0*${JH}) )  ) "  | bc -l)
SQ1_by2_=$(echo "0.5*${SQ1_}"  | bc -l)

SQ2_=$(echo "-1.0*((${t}*${t})/(${U}))*(  ( (${JH})/(${U}-3.0*${JH}) )  ) "  | bc -l)
SQ2_by2_=$(echo "0.5*${SQ2_}"  | bc -l)

_mSQ2_by2_=$(echo "-1.0*${SQ2_by2_}"  | bc -l)
_mSQ2_=$(echo "-1.0*${SQ2_}"  | bc -l)

# extra minus sign because QyzI=iota X Qyz
SQ3_=$(echo "((-0.5*(${t}*${t}))/(${U}))*(  ( (${U}-${JH})/(${U}-3.0*${JH})  )  - ${fac_}*( (${JH})/(${U}+2.0*${JH}) )  ) "  | bc -l)
SQ3_by2_=$(echo "0.5*${SQ3_}"  | bc -l)


SL_=$(echo "((0.5*${t}*${t})/(${U}))*(  ( (${U}-${JH})/(${U}-3.0*${JH})  )  + ${fac_}*( (${JH})/(${U}+2.0*${JH}) )  ) "  | bc -l)
SL_by2_=$(echo "0.5*${SL_}"  | bc -l)


QQ1_=$(echo "((-1.0*${t}*${t})/(${U}))*(  ( (3.0*${U}-${JH})/(${U}-3.0*${JH})  )   )"  | bc -l)


QQ2_=$(echo "((-1.0*${t}*${t})/(${U}))*(  ( (${U}+${JH})/(${U}+2.0*${JH})  )  + ( (${JH})/(${U}-3.0*${JH}) )  )"  | bc -l)

QQ3_=$(echo "((${t}*${t})/(${U}))*(  ( (2.0*${U}-${JH})/(${U}-3.0*${JH})  )  ) "  | bc -l)

QQ4_=$(echo "(-1.0*(${t}*${t})/(${U}))*(  ( (2.0*${U}-${JH})/(${U}-3.0*${JH})  )  )"  | bc -l)

# extra minus sign because QyzI=iota X Qyz
QQ5_=$(echo "((-0.5*(${t}*${t}))/(${U}))*(  ( (${U}+${JH})/(${U}-3.0*${JH})  )  + ${fac_}*( (${JH})/(${U}+2.0*${JH}) )  ) "  | bc -l)

LL_=$(echo "((0.5*${t}*${t})/(${U}))*(  ( (${U}+${JH})/(${U}-3.0*${JH})  )  - ${fac_}*( (${JH})/(${U}+2.0*${JH}) )  ) "  | bc -l)



sed -i -e "s/_mSQ2_by2_/${_mSQ2_by2_}/g" ${file}
sed -i -e "s/_mSQ2_/${_mSQ2_}/g" ${file}

sed -i -e "s/SQ1_by2_/${SQ1_by2_}/g" ${file}
sed -i -e "s/SQ2_by2_/${SQ2_by2_}/g" ${file}
sed -i -e "s/SQ3_by2_/${SQ3_by2_}/g" ${file} 
sed -i -e "s/SL_by2_/${SL_by2_}/g" ${file}

sed -i -e "s/SQ1_/${SQ1_}/g" ${file}
sed -i -e "s/SQ2_/${SQ2_}/g" ${file}
sed -i -e "s/SQ3_/${SQ3_}/g" ${file}
sed -i -e "s/SL_/${SL_}/g" ${file}
sed -i -e "s/QQ1_/${QQ1_}/g" ${file}
sed -i -e "s/QQ2_/${QQ2_}/g" ${file}
sed -i -e "s/QQ3_/${QQ3_}/g" ${file}
sed -i -e "s/QQ4_/${QQ4_}/g" ${file}
sed -i -e "s/QQ5_/${QQ5_}/g" ${file}
sed -i -e "s/LL_/${LL_}/g" ${file}




echo "${U}  ${JH} ${SQ1_} ${SQ2_} ${SQ3_} ${SL_} ${QQ1_} ${QQ2_} ${QQ3_} ${QQ4_} ${LL_}"

