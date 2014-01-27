//===================================================================
//    Atmosferic Formula
//      Author    : T.Shibata
//      Creation  : 2010.11.15
//  Last Update   : 2010.11.15
//  Last Update   : 2011.12.13
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>

#include "G4ThreeVector.hh"

#ifndef AtmosfericFormula_h
#define AtmosfericFormula_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class AtmForm{

public:

  AtmForm();
  ~AtmForm();
  
  double teme_Kelvin( double temerature_deg );
  double teme_DegC( double temerature_k );

  double Saturated_vapor_pressure_tetens( double temerature_deg );
  double Saturated_vapor_pressure_weler_hyland( double temerature_deg );        

  double Vapor_pressure( double relative_humidity, 
                         double saturated_vapor_pressure );
  double Vapor_pressure_rate( double air_pressure,
                              double vapor_pressure );  
  
  double Air_mass_density( double temerature_deg,
			   double air_pressure,
                           double vapor_pressure );

  double N2Air_mass_density( double temerature_deg,
			     double air_pressure );
 
private:



protected:



};
#endif

