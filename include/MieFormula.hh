//===================================================================
//     Mie Scattering Formula
//      Author    : T.Shibata
//
//      Creation  : 2013.09.09
//      Modified  : 2013.09.09
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

#ifndef MieFormula_h
#define MieFormula_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class Mie{

public:

  //Contstruction
   Mie();
  //Deconstruction
  ~Mie();

  void attenuation_length( const int Num, 
                           double in_lambda_min,
                           double in_lambda_max,
                           double in_temperature,
                           double in_pressure,
                           double in_pressure_water,
                           double in_rhomass,
                           vector<double> & Mlambda,
                           vector<double> & Attenuation );

  double Forward_g(void);
  double Backward_g(void);
  double ForwardBackward_r(void);

private:

  double lambda_min;
  double lambda_max;
  double dlambda;
  double temperature;  
  double pressure;
  double pressure_water;
  double rhomass;

  double AttenuationCoefficientTAJAVA( double lambda );

};
#endif
