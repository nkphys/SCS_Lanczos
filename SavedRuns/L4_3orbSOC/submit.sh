for lambda in 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.4 0.41 0.42 0.43 0.44 0.45 0.46 0.48 0.49 0.5 0.52 0.54 0.55 0.6 0.7
do

mkdir -p lambda${lambda}
cd lambda${lambda}

cp ../{input3orb_SOC.inp,lanczos} .
sed -i -e "s/LAMBDA_VALUE/${lambda}/g" input3orb_SOC.inp

echo "lambda=${lambda} running"
time ./lanczos input3orb_SOC.inp > out_run.txt

cd ..
done
