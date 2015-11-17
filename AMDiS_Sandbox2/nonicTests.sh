# sed -ie '/nonic99p_C/s/[[:digit:]][[:punct:]][[:digit:]]*/0.9943/' init/octic.dat.3d 

Cs=(0.{5..9}00 1.{0..5}00 1.0{1..9}0 1.09{1..9})

for C in ${Cs[*]}
do
  sleep 5
  echo $C
  sed -i -e '/nonic95p_C/s/[[:digit:]][[:punct:]][[:digit:]]*/'$C'/' init/octic.dat.3d
  sed -i -e 's/sphere[[:digit:]][[:punct:]][[:digit:]]*/sphere'$C'/' init/octic.dat.3d
  #sleep `shuf -i10-30 -n1`&
  nohup ./directorFieldOctic init/octic.dat.3d &> nohups/nohup${C}.out &
  while [ `ps a | grep -c octic` -gt 6 ]
  do
    #echo '...waiting...'
    sleep 100
  done
done
