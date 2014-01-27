#!/bin/bash -f
########!/bin/tcsh -f

  NEV=$1
  SETF=$2

  GENF=./dat-set/init-generator.dat
  OF=./detector-plane.dat  
  GF=./generator.dat
  KF=./dedx.dat
  FF=./faradaycup.dat
  EF=./energydeposit.dat
  MF=./mattOutput.dat
  

  # Geant4 Simulation exec
  echo "Geant4 Simulation is started"
  #echo elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF -o $OF -g $GF -k $KF -f $FF -e $EF -m $MF 
  elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF -o $OF -g $GF -k $KF -f $FF -e $EF -m $MF >& ./log #/dev/null
  #elsbeam-geant4 -v 0 -d $NEV -i $GENF -s $SETF -o $OF -g $GF -k $KF -f $FF -e $EF -m $MF 

  echo "Geant4 Simulation is finished"
 
  #-------------------------
  # Faraday Cup ChargeCount
  #-------------------------
  FARADAYCUP=./faradaycup_capture.dat
  ./c++/faradaycup_chargecount $FF $FARADAYCUP $NEV

exit

