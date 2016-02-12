IDS=({41..49})

IF=init/nonicPressed.dat.3d

for ID in ${IDS[*]}
do
  sleep 5
  echo $ID
  initSet.py -f $IF -k 'macro file name' -v ./macro/nonicsPressed/`nonicName.py --id $ID`_64k.3d

  initSet.py -f $IF -k 'output->filename' -v outNonicPressed/`nonicName.py --id $ID`_2D
  initSet.py -f $IF -k 'initField' -v rotated_ey
  nohup ./directorField init/nonicPressed.dat.3d &> nohups/nohup_${ID}_2D.out &
  while [ `ps x | grep -c directorField` -gt 5 ]
  do
    sleep 100
  done

  sleep 5
  initSet.py -f $IF -k 'output->filename' -v outNonicPressed/`nonicName.py --id $ID`_4D
  initSet.py -f $IF -k 'initField' -v ex
  nohup ./directorField init/nonicPressed.dat.3d &> nohups/nohup_${ID}_4D.out &
  while [ `ps x | grep -c directorField` -gt 5 ]
  do
    sleep 100
  done
done
