//===================================================================
//    Electron Beam Initial Information :: PrimaryEbeam.hh
//      Author    : T.Shibata
//      Creation  : 2006.02.13
//  Last Update   : 2010.02.11
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

#include "G4ThreeVector.hh"

#include "Randomize.hh"

#ifndef PrimaryEbeam_h
#define PrimaryEbeam_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class Ebeam {

public:
      double gaussian( double mean, double sigma ); 
      double double_gaussian( double fluction, 
                              double mean1, double sigma1,
                              double mean2, double sigma2 ); 
      double rn_radian(void);

      double pulse( double xs, double xe ); 
      int poisson( double mean, double xmax, double xmin );
      double probability_poisson( int n, double m );

      double bunch_time( int Nevent );

      double crystal_ball( double mean,
	   	           double sigma,
		           double alpha,
		           double n,
		           double xmin, double xmax );

      double probability_crystal_ball( double x,
	  	 	 	       double mean,
				       double sigma,
				       double alpha,
				       double n );
        
      double emittance_position( double emmitance[], double mean );
      double emittance_direction( double emmitance[], double posi );

      void special_beam_direction_20101102( double beam_position[],
                                            double &beam_direction_x,
					    double &beam_direction_y,
                                            double &beam_direction_z );
                                     
      //Contstructione
      Ebeam();

      //Deconstruction
      ~Ebeam();

private:
      double pi;
      double rn(void);

       int Nppb;
       double bunchtime;

}; 
#endif
