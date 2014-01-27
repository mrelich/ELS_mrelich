#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>
#include <string.h>
#include <signal.h>

//============================================
using namespace std;
//============================================
int main( int argc, char** argv )
{

  // Read File ==============================================
  char FaradayCupFileName[256]; 
  sprintf( FaradayCupFileName, argv[1] );

  FILE *fp_faradaycup;
  if(!(fp_faradaycup=fopen(FaradayCupFileName,"r"))){
    printf("error:Can not Open File %s\n", FaradayCupFileName );
    exit(1);
  }

  // Write File ==============================================
  char FaradayCupFileName_new[256];
  sprintf( FaradayCupFileName_new, argv[2] );

  FILE*fp_faradaycup_new;
  if(!(fp_faradaycup_new=fopen(FaradayCupFileName_new,"w"))){
    printf("error:Can not Open File %s\n", FaradayCupFileName_new  );
    exit(1);
  }

  //===========================================================   
    int Number_of_Beamelectron = atoi(argv[3]);
  //===========================================================   

  int c(0);
  int t1, t2, t3;
  int PID, TID, ChargeDeposit ;  

  int TotalCharge(0);
  int PrimaryCharge(0);
  int SecondaryCharge(0);

  float Capture_efficiency_total(0.);  
  float Capture_efficiency_primary(0.);  
  float Capture_efficiency_secondary(0.);  

  int dummy_d;
  float dummy_f;

  while(1){
    c=fscanf(fp_faradaycup, "%d %d %d %d %d %f %f %f\n", &dummy_d, &dummy_d, 
                                                         &t1, &t2, &t3, 
                                   	                     &dummy_f, &dummy_f, &dummy_f );
    if(c==-1) break;

    PID=t1;
    TID=t2;
    ChargeDeposit=t3;    
    TotalCharge+=ChargeDeposit;
    if( TID==1 ) PrimaryCharge+=ChargeDeposit;
    if( TID!=1 ) SecondaryCharge+=ChargeDeposit;
  }

  Capture_efficiency_total     = float((-1)*TotalCharge)/float(Number_of_Beamelectron)*100.0;
  Capture_efficiency_primary   = float((-1)*PrimaryCharge)/float(Number_of_Beamelectron)*100.0;
  Capture_efficiency_secondary = float((-1)*SecondaryCharge)/float(Number_of_Beamelectron)*100.0;

  fprintf( fp_faradaycup_new, "%f %f %f\n",
	   Capture_efficiency_total,
	   Capture_efficiency_primary,
	   Capture_efficiency_secondary );

  fclose(fp_faradaycup);
  fclose(fp_faradaycup_new);

  return 0;
}
