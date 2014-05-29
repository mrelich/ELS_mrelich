#!/bin/bash -f



  #######################################################
  # Shower Generator for ELS RUN 0
  #######################################################

  NEV=2 # Number of events to generate
  #NEV=300 # Number of events to generate
  SETF=./dat-set/setupfile-elssim.dat
  NP=1 # Number of incident electron

  date 

  ./elssim.sh $NEV $SETF $NP

  date 

exit

