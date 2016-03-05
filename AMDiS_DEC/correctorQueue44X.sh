IDS=({441..449})

for ID in ${IDS[*]}
do
  sleep 2
  echo Start nonic mesh $ID: `date`
  typeset -i IDOLD
  IDOLD=$ID-1
  initSet.py -f init/meshCorrector2.dat.3d -k 'macro file name' -v ../AMDiS_Sandbox2/macro/nonicsPressed/`nonicName.py --id $IDOLD`_64k.3d
  initSet.py -f init/meshCorrector2.dat.3d -k 'outName' -v `nonicName.py --id $ID`
  initSet.py -f init/meshCorrector2.dat.3d -k 'octic->c' -v `nonicName.py --id $ID -gc`
  initSet.py -f init/meshCorrector2.dat.3d -k 'nonic->press' -v `nonicName.py --id $ID -gb`
  initSet.py -f init/meshCorrector2.dat.3d -k 'nonic->old->c' -v `nonicName.py --id $IDOLD -gc`
  initSet.py -f init/meshCorrector2.dat.3d -k 'nonic->old->press' -v `nonicName.py --id $IDOLD -gb`
  nohup ./meshCorrector init/meshCorrector2.dat.3d &> nohups/nohup_${ID}.out
  sleep 1
  cp output/meshOut`nonicName.py --id $ID`_250.3d ../AMDiS_Sandbox2/macro/nonicsPressed/`nonicName.py --id $ID`_64k.3d
done
