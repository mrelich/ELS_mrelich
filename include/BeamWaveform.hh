//===================================================================
//    Electron Beam Initial Information :: BeamWaveform.hh
//      Author    : T.Shibata
//      Creation  : 2011.12.13
//  Last Update   : 2011.12.13
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

#include "Randomize.hh"

#include "ReadFile.hh"

#ifndef BeamWaveform_h
#define BeamWaveform_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class Plist;
class Waveform {

public:
  /* Waveform Array*/
  int Ntbin;
  int NewNtbin;
  double dt;
  double Newdt;
  double InputWaveform[10000];
  double NewWaveform[10000];
  double NormalizedNewWaveform[10000];
  double NormalizationFactor;

  double rn(void);   

  int ReadWaveformFile( Plist *plist );

  double Output_time(void);

  //Contstructione
  Waveform();
  //Deconstruction
  ~Waveform();

private:

};
#endif
