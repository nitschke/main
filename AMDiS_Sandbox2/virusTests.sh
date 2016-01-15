#Cs=({10..70..5})
Cs=({10..35..5})
Cs+=({45..70..5})
echo ${Cs[*]}

for C in ${Cs[*]}
do
  sleep 5
  echo $C
  sed -i -e '/macro\/virus\/virus/s/virus[[:digit:]]*kS/virus'$C'kS/' init/octic.dat.3d
  sed -i -e '/output\/virus/s/virus[[:digit:]]*/virus'$C'/' init/octic.dat.3d
  nohup ./directorFieldOctic init/octic.dat.3d &> nohups/nohupVirus${C}kS.out &
  while [ `ps a | grep -c octic` -gt 6 ]
  do
    #echo '...waiting...'
    sleep 100
  done
done
