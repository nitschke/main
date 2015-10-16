# sed -ie '/nonic99p_C/s/[[:digit:]][[:punct:]][[:digit:]]*/0.9943/' init/octic.dat.3d 

Cs=(0.97997 0.97998 0.97999 0.982 0.983 0.984 0.985 0.986 0.987 0.988 0.989 0.990 0.960 1.000)

for C in ${Cs[*]}
do
  sleep 5
  echo $C
  sed -ie '/nonic99p_C/s/[[:digit:]][[:punct:]][[:digit:]]*/'$C'/' init/octic.dat.3d
  sed -ie 's/sphere[[:digit:]][[:punct:]][[:digit:]]*/sphere'$C'/' init/octic.dat.3d
  #sleep `shuf -i10-30 -n1`&
  nohup ./directorFieldOctic init/octic.dat.3d &> nohup${C}.out &
  while [ `ps a | grep -c octic` -gt 7 ]
  do
    #echo '...waiting...'
    sleep 100
  done
done
