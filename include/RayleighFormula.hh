//===================================================================
//     Rayleigh Scattering Formula
//      Author    : T.Shibata
//
//      Creation  : 2006.10.26
//      Modified  : 2012.09.24
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

#ifndef RayleighFormula_h
#define RayleighFormula_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class Rayleigh{

public:

  //Contstruction
   Rayleigh();
  //Deconstruction
  ~Rayleigh();

  void attenuation_length( const int Num, 
                           double in_lambda_min,
                           double in_lambda_max,
                           double in_temperature,
                           double in_pressure,
                           double in_pressure_water,
                           double in_rhomass,
                           vector<double> & Mlambda,
                           vector<double> & RefractionIndex,
                           vector<double> & Attenuation );

private:

  double lambda_min;
  double lambda_max;
  double dlambda;
  double temperature;  
  double pressure;
  double pressure_water;
  double rhomass;

  double Refraction(void);
  double RefractionCorrection( double lambda, double temp, double pres, double pres_water );
  double AttenuationCoefficient( double temp, double refraction, double lambda );
  double AttenuationCoefficientTAJAVA( double rho, double lambda );

};
#endif
