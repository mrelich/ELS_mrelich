//===================================================================                                         
//    Atmosferic Formula                                                                                      
//      Author    : T.Shibata                                                                                 
//      Creation  : 2010.11.15                                                                                
//  Last Update   : 2010.11.15
//  Last Update   : 2011.12.13 
//===================================================================      
#include "PhysicsParameters.hh"
#include "AtmosfericFormula.hh"
//-----------------------------------------------                                                      
using namespace std;
//------------------------------------------------
AtmForm::AtmForm(){}
AtmForm::~AtmForm(){}
//------------------------------------------------                  
double AtmForm::teme_Kelvin( double temerature_deg )
{
  return temerature_deg-Kelvin;
}
//------------------------------------------------
double AtmForm::teme_DegC( double temerature_k )
{
  return temerature_k+Kelvin;
}
//-----------------------------------------------------------------------
// Cacluration of saturated_vapor_pressure
// Model Formula : Tetens formula(1930) 
// input  = deg-C
// output = hPa                     
double AtmForm::Saturated_vapor_pressure_tetens( double temerature_deg )
{
  return 6.11*pow( 10, 7.5*temerature_deg/(temerature_deg+237.3) );
}
//-----------------------------------------------------------------------
// Cacluration of saturated_vapor_pressure
// Model Formula : Weler-Hyland (1983)
// input  = deg-C/deg-K
// output = hPa                     
double AtmForm::Saturated_vapor_pressure_weler_hyland( double temerature_deg )
{
  double temerature_k=teme_Kelvin(temerature_deg);

  double a[2][6];

  a[0][0] = -0.56745359;
  a[0][1] =  0.63925247;
  a[0][2] = -0.96778430;
  a[0][3] =  0.62215701;
  a[0][4] =  0.20747825;
  a[0][5] =  0.41635019;

  a[1][0] = -0.58002206;
  a[1][1] =  0.13914993;
  a[1][2] = -0.48640239;
  a[1][3] =  0.41764768;
  a[1][4] = -0.14452093;
  a[1][5] =  0.65459673;

  double b[2][6];
  
  b[0][0] = 1.0E4;
  b[0][1] = 1.0E1;
  b[0][2] = 1.0E-2; 
  b[0][3] = 1.0E-6; 
  b[0][4] = 1.0E-12; 
  b[0][5] = 1.0E1;

  b[1][0] = 1.0E4;
  b[1][1] = 1.0E1;
  b[1][2] = 1.0E-1;
  b[1][3] = 1.0E-4; 
  b[1][4] = 1.0E-7; 
  b[1][5] = 1.0E1;

  int i(0);
  if ( temerature_deg > 0.01 ){
    i=1;
  }else if( temerature_deg <= 0.01 ){
    i=0;
  }

  return exp( a[i][0]*b[i][0]/temerature_k + 
              a[i][1]*b[i][1] +
              a[i][2]*b[i][2]*temerature_k +
              a[i][3]*b[i][3]*temerature_k*temerature_k +
              a[i][4]*b[i][4]*temerature_k*temerature_k*temerature_k +
              a[i][5]*b[i][5]*log(temerature_k) )/100.0;       // Pa -> hPa
}
//-----------------------------------------------------------------------
// Cacluration of vapor pressure (hPa)                                                                                 
// input  = %                                                                                                         
// output = hPa        
double AtmForm::Vapor_pressure( double relative_humidity, 
                                double saturated_vapor_pressure )
{
  return relative_humidity*saturated_vapor_pressure*0.01;
}
//-----------------------------------------------------------------------
// Air composition ratio
// output = % 
double AtmForm::Vapor_pressure_rate( double air_pressure,
				     double vapor_pressure )
{
  return (air_pressure-vapor_pressure)/air_pressure;
}
//-----------------------------------------------------------------------    
// Air mass density 
// Rika-Nenpyou. 2010.Phys.26(376)
// output = g/cm3
double AtmForm::Air_mass_density( double temerature_deg,
				  double air_pressure,
				  double vapor_pressure )
{
 double air_mass_density 
   = 0.001293/(1.+0.00367*temerature_deg)*air_pressure/1013.2472;
 air_mass_density *= (1.0-0.378*vapor_pressure/air_pressure);

 return air_mass_density;
}
//-----------------------------------------------------------------------    
// N2 Air mass density
// output = g/cm3 
double AtmForm::N2Air_mass_density( double temerature_deg,
                                    double n2air_pressure )
{
  double N2air_mass_density = n2air_pressure*100.0/(296.92*teme_Kelvin(temerature_deg)); // unit=kg/m^3
  N2air_mass_density*=1.0/1000.0; // unit = g/cm^3
        
  return N2air_mass_density;
}

