//==============================================================================
//    Photon Generator :: Geo-Magnetic Field
//      Author    : T.Shibata
//      Creation  : 2007.06.07
//    Last Update : 2007.06.07
//    Last Update : 2012.01.20
//==============================================================================
#ifndef GeoMagneticField_H
#define GeoMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

class GeoMagnetic : public G4MagneticField
{
  public:
  GeoMagnetic();
  ~GeoMagnetic();

  void Geomagneticfield_onoff( int onoff );
  void InputGeomagneticInfomation( double Theta, double D,
                                   double Bh,    double Bv );

  void GetFieldValue( const  double Point[3], double *Bfield ) const;

  private:
  G4int    OnOff; 

  G4double BiasAngle;         // unit = dgree

  G4double Declination;       // unit = dgree
  G4double BH;
  G4double BV;

  G4double Bx,By,Bz;

};

#endif
