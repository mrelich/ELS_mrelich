//===================================================================
//     Rayleigh Scattering Formula
//      Author    : T.Shibata
//
//      Creation  : 2006.10.26
//      Modified  : 2012.09.24
//
//===================================================================
#include "RayleighFormula.hh" 
#include "PhysicsParameters.hh"
#include <iomanip>
//-----------------------------------------------
using namespace std;
//------------------------------------------------
Rayleigh::Rayleigh(){}
//------------------------------------------------
Rayleigh::~Rayleigh(){}
//------------------------------------------------
void Rayleigh::attenuation_length( const int in_Num,
				   double in_lambda_min,
				   double in_lambda_max,
				   double in_temperature,
				   double in_pressure,
				   double in_pressure_water,
				   double in_rhomass,
                   vector<double> & Mlambda,
				   vector<double> & RefractionIndex,
  		    	   vector<double> & Attenuation )
{

   const int Num = in_Num;

   lambda_min   = in_lambda_min;
   lambda_max   = in_lambda_max;
   temperature  = in_temperature;
   pressure     = in_pressure; 
   pressure_water = in_pressure_water;
      
   rhomass      = in_rhomass;

   if( in_lambda_max < in_lambda_min ){
     double tmp = in_lambda_max;
     in_lambda_max = in_lambda_min;
     in_lambda_min = tmp;
   }

   dlambda = (in_lambda_max - in_lambda_min )/Num;  

   for( int i=0; i<Num; i++ ){

     double ilambda = in_lambda_min + (double(i)-1)*dlambda;
     //double refraction = Refraction( ilambda, pressure, temperature );
     //double iattenuationcoeff = AttenuationCoefficient( temperature, refraction, ilambda );
     double irefraction = RefractionCorrection( ilambda, temperature, pressure, pressure_water );
	 double iattenuationcoeff= AttenuationCoefficientTAJAVA( in_rhomass, ilambda );
     //cout << " ray= " << ilambda << " " << setprecision(15) << irefraction << endl;

     Mlambda.push_back(ilambda); 
     RefractionIndex.push_back(irefraction);
     Attenuation.push_back(1.0/iattenuationcoeff);

   }   

}
//------------------------------------------------
/*

*/
double Rayleigh::Refraction( void )
{
  return 1.000292; 
}
//------------------------------------------------
/*
   Owens formula :Applied Optics, 6, 51 (1967) 
   
   CO2 = 0.03% contanimation
*/
double Rayleigh::RefractionCorrection( double lambda,
                                       double temp,
                                       double pres,
                                       double pres_water )
{
  double p_a = pres;
  double p_w = pres_water;
  double T   = temp;

  double wn=1.0/(lambda/1000.0); // um^-1  
  double Ds
	=(p_a/T)*
	( 1.0 + p_a*( 57.90E-8 - 9.3250E-4/T + 0.25844/T/T ) );
  double Dw
	=(p_w/T)*
	( 1.0 + p_w*( 1.0 + 3.7E-4*p_w )*( -2.37321E-3+2.23366/T-710.792/T/T+7.75141E4/T/T/T ));
  double dn
	= (2371.34+683939.7/(130.0-wn*wn)+4547.3/(38.9-wn*wn))*Ds
	+ (6487.31+58.058*wn*wn-0.71150*pow(wn,4)+0.08851*pow(wn,6))*Dw;

  return 1.0+dn/1.E8;
}
//------------------------------------------------
/*
 Einstein-Smoluchowski Formula ( 1910 )
 
*/
double Rayleigh::AttenuationCoefficient( double temp, 
                                         double refraction, 
                                         double lambda )
{

  const double  kT    ( 9.87E-7 );  // cm^2/dyn = cm * s^2 * g^-1

  lambda=lambda*1E-9;   // nm --> m

  double f1 = 8.0*PI*PI*PI/27.0;
  double f2 = k_boltzman*temp*kT;
  double f3 =  ( refraction*refraction - 1 )*( refraction*refraction - 1 )
              *( refraction*refraction + 2 )*( refraction*refraction + 2 );
  double f4 = 1.0/(lambda*lambda*lambda*lambda);
  double f5 = 10.0;    // <-- Dimension Correction
 
  double rn = f1*f2*f3*f4*f5;

  return rn;  // output is "1/m" 
}
//------------------------------------------------
/*
 TA-JAVA
 /sources/telescopeArray/atmosphere/MeasuredAtmosphere.java  
*/
double Rayleigh::AttenuationCoefficientTAJAVA( double rho,  // g/cm^3
                                               double lambda ) // nm 
{
  double alpha=pow(400.0/lambda, 4 )*rho/2974.0;
  return alpha*1.E2; // output is "1/m"
}
