//==============================================================================
//    Photon Generator :: Geo-Magnetic Field
//      Author    : T.Shibata
//      Creation  : 2007.06.07
//    Last Update : 2007.06.07
//    Last Update : 2012.01.20
//==============================================================================
#include "GeoMagneticField.hh"
#include "PhysicsParameters.hh"

//-----------------------------------------------
using namespace std;
//------------------------------------------------

//-------------------------------------------------------------------------
GeoMagnetic::GeoMagnetic()
{
  OnOff       = 0;

  BiasAngle   = 0.;

  Declination = 0.;
  BH = 0.0*tesla;
  BV = 0.0*tesla;

  Bx = 0.0*tesla;
  By = 0.0*tesla;
  Bz = 0.0*tesla;
}
//-------------------------------------------------------------------------
GeoMagnetic::~GeoMagnetic(){;}
//-------------------------------------------------------------------------
void GeoMagnetic::Geomagneticfield_onoff( int onoff )
{
  OnOff = onoff;
}
//-------------------------------------------------------------------------
void GeoMagnetic::InputGeomagneticInfomation( double Theta, 
                                              double D, double Bh, double Bv )
{
  BiasAngle   = Theta*RADI;
  
  Declination = D*RADI;
  BH          = Bh;
  BV          = Bv;
}
//-------------------------------------------------------------------------
void GeoMagnetic::GetFieldValue(const double Point[3],double *Bfield) const
{
  Bfield[0] = 0;
  Bfield[1] = 0;
  Bfield[2] = 0;

  if( OnOff ){
      Bfield[0] =       BH*cos( BiasAngle + Declination ) * tesla ;
      Bfield[1] = (-1.)*BH*sin( BiasAngle + Declination ) * tesla ; 
      Bfield[2] = BV * tesla ;
  }
 
}
//-------------------------------------------------------------------------
