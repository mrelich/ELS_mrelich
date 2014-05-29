#!/bin/bash -f
########!/bin/tcsh -f

  NEV=$1
  SETF=$2
  NP=$3

  GENF=./dat-set/init-generator.dat
  #RO="-e -q"
  RO="-st"
  

  # Geant4 Simulation exec
  echo "Geant4 Simulation is started"
  #echo elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF -o $OF -g $GF -k $KF -f $FF -e $EF -m $MF 
  nohup elsbeam-geant4 -v 0 -ne $NEV -np $NP -i $GENF -s $SETF $RO >& /dev/null
  #elsbeam-geant4 -v 0 -ne $NEV -np $NP -i $GENF -s $SETF $RO > log
  #echo elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF $RO
  #elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF -o $OF -g $GF -k $KF -f $FF -e $EF -m $MF 

  echo "Geant4 Simulation is finished"
 
  #-------------------------
  # Faraday Cup ChargeCount
  #-------------------------
  #FARADAYCUP=./faradaycup_capture.dat
  #./c++/faradaycup_chargecount $FF $FARADAYCUP $NEV

exit

