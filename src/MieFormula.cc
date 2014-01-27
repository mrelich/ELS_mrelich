//===================================================================
//     Mie Scattering Formula
//      Author    : T.Shibata
//
//      Creation  : 2013.09.09  
//      Modified  : 2013.09.09  
//
//===================================================================
#include "MieFormula.hh" 
#include "PhysicsParameters.hh"
#include <iomanip>
//-----------------------------------------------
using namespace std;
//------------------------------------------------
Mie::Mie(){}
//------------------------------------------------
Mie::~Mie(){}
//------------------------------------------------
void Mie::attenuation_length( const int in_Num,
				              double in_lambda_min,
				              double in_lambda_max,
				              double in_temperature,
				              double in_pressure,
				              double in_pressure_water,
				              double in_rhomass,
                              vector<double> & Mlambda,
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
	 double iattenuationcoeff= AttenuationCoefficientTAJAVA(ilambda);

     //cout << " refraction = " << ilambda << " " << setprecision(15) << 1./iattenuationcoeff << endl;

     Mlambda.push_back(ilambda); 
     Attenuation.push_back(1.0/iattenuationcoeff);

   }   

}
//------------------------------------------------
double Mie::Forward_g(void)
{
  return 0.970;
}
double Mie::Backward_g(void)
{
  return 0.750;
}
double Mie::ForwardBackward_r(void)
{
  return 0.9990;
}
//------------------------------------------------
double Mie::AttenuationCoefficientTAJAVA(  double lambda ) // nm 
{
  double alpha=(355.0/lambda)/29400.0;  // m
  return alpha; // output is "1/m"
}
//------------------------------------------------

