#!/bin/bash -f
#########!/bin/tcsh -f


  #######################################################
  # Shower Generator for ELS RUN 0
  #######################################################

  NEV=1 # Number of Incident electron
  #NEV=10 # Number of Incident electron
  SETF=./dat-set/setupfile-elssim.dat

  date 

  ./elssim.sh $NEV $SETF

     rm ./detector-plane.dat
     #rm ./dedx.dat
     rm ./generator.dat
     #rm ./energydeposit.dat
     rm ./screenmonitor3.dat
     rm ./initialbeamenergy.dat
     #rm ./cerenkov.dat
     #rm ./debug-output.dat
     #rm ./geant4-air-condition.dat
     rm ./faradaycup_capture.dat  
     rm ./faradaycup.dat

     #rm ./eventnumber.dat
     

  date 

exit

