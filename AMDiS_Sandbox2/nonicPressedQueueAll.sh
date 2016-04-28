IDS=({2..9})
IDS+=({41..49})
IDS+=({71..79})
IDS+=({431..439})
IDS+=({441..449})
IDS+=({4331..4339})

IF=init/nonicPressed.dat.3d

for ID in ${IDS[*]}
do
  sleep 5
  initSet.py -f $IF -k 'macro file name' -v ./macro/nonicsPressed/`nonicName.py --id $ID`_64k.3d
  initSet.py -f $IF -k 'userParameter->B' -v `nonicName.py --id $ID -gb`
  initSet.py -f $IF -k 'userParameter->C' -v `nonicName.py --id $ID -gc`

  echo Start $ID With 2 Defects: `date`
  initSet.py -f $IF -k 'output->filename' -v outNonicPressed/`nonicName.py --id $ID`_2D.
  initSet.py -f $IF -k 'initField' -v rotated_ey
  nohup ./directorField init/nonicPressed.dat.3d &> nohups/nohup_${ID}_2D.out &
  while [ `ps x | grep -c directorField` -gt 6 ]
  do
    sleep 100
  done

  sleep 5
  echo Start $ID With 4 Defects: `date`
  initSet.py -f $IF -k 'output->filename' -v outNonicPressed/`nonicName.py --id $ID`_4D.
  initSet.py -f $IF -k 'initField' -v ex
  nohup ./directorField init/nonicPressed.dat.3d &> nohups/nohup_${ID}_4D.out &
  while [ `ps x | grep -c directorField` -gt 6 ]
  do
    sleep 100
  done
done
