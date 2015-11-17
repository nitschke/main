Cs=(0.980 1.000 0.960)

for C in ${Cs[*]}
do
  sleep 5
  sed -ie '/nonic99p_C/s/[[:digit:]][[:punct:]][[:digit:]]*/'$C'/' init/octic.dat.3d
  sed -ie 's/sphere[[:digit:]][[:punct:]][[:digit:]]*/sphere'$C'/' init/octic.dat.3d
  for SEED in {1..25}
  do
    sleep 5
    echo 'start simulation (C = '$C') with seed: '$SEED
    sed -ie '/seed/s/ [[:digit:]]*/ '$SEED'/' init/octic.dat.3d
    nohup ./directorFieldOctic init/octic.dat.3d &> nohups/nohup_${C}_${SEED}.out &
    while [ `ps a | grep -c octic` -gt 6 ]
    do
      #echo '...waiting...'
      sleep 100
    done
  done
done
