//===================================================================
//    Electron Beam Initial Information :: BeamWaveform.cc
//      Author    : T.Shibata
//      Creation  : 2011.12.13
//  Last Update   : 2011.12.13
//===================================================================
#include "BeamWaveform.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
Waveform::Waveform(){

  time_t  t;
  srand( (unsigned) time(NULL) );

  Ntbin = 10000;
  NewNtbin = 10000;
  dt=1.; // unit=1ns
  Newdt=1.;
  NormalizationFactor=1.;

  for( int i=0; i<10000; i++ ){
    InputWaveform[i]=NewWaveform[i]=NormalizedNewWaveform[i]=0.0;
  }

}
//------------------------------------------------
Waveform::~Waveform(){}
//------------------------------------------------
double Waveform::rn(void)
{
  return double(rand())/double(RAND_MAX);
}
//------------------------------------------------
int Waveform::ReadWaveformFile( Plist *plist )
{
  FILE *fp_waveform;

  if(!( fp_waveform = fopen( plist->ElectronBeam_Waveform_File() ,"r"))){
    return -1;
  }

  int c(0);
  float tmp_t;
  float tmp_v;
  int t;
  int i(0);
  while(1){
    c=fscanf( fp_waveform, "%f %f\n",&tmp_t,&tmp_v );
    if(c==-1) break;
    t=int(tmp_t);   
    InputWaveform[i]=double(tmp_v);  
    i++;
  }
  fclose(fp_waveform);
  
  Ntbin = i;   
  dt = plist->ElectronBeam_Waveform_dt();  // unit = nsec

  NewNtbin = plist->ElectronBeam_Waveform_new_ntbin();
  Newdt = dt*double(Ntbin)/double(NewNtbin);

  for( int i=0; i<Ntbin; i++ ){
       NewWaveform[i/int(Ntbin/NewNtbin)]+=InputWaveform[i];
  }

  double max(0.);
  for( int i=0; i<NewNtbin; i++ ){
    if( NewWaveform[i] > max ) max=NewWaveform[i]; 
  }
  NormalizationFactor=1.2*max;
  for( int i=0; i<NewNtbin; i++ ){
    NormalizedNewWaveform[i]=NewWaveform[i]/NormalizationFactor;
  }       

  return 0;
}
//------------------------------------------------
double Waveform::Output_time(void)
{

  double rn_t;
  double y_rn;
  while(1){
    rn_t= rn()*double(NewNtbin);
    y_rn=rn()*NormalizationFactor;
    if( y_rn < NormalizedNewWaveform[int(rn_t)] ) break;       
  }

  return rn_t*Newdt;
}
//------------------------------------------------
