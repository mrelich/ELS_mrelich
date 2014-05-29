//==============================================================================
//    ELS Magnetic Field
//      Author    : T.Shibata
//      Creation  : 2009.12.08
//    Last Update : 2010.11.01
//==============================================================================
#ifndef ELSMagneticField_H
#define ELSMagneticField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"

#include "ELSParameters.hh"
#include "G4SystemOfUnits.hh"


class ELSParameters;
class ELSMagnetic : public G4MagneticField
{
  public:
  ELSMagnetic();
  ~ELSMagnetic();

  void elsmagneticfield_initialize( ELSParameters *elsparam );

  void elsmagneticfield_onoff( int QM_onoff,
                               double QM1_field, double QM2_field,
                               int BM_onoff,
                               double BM_filed );

  void GetFieldValue( const  double Point[3],
		      double *Bfield ) const;

  private:

  double c0;
  double me;
  double T0;
  double p0;
  double r0;

  G4double Bx,By,Bz;

  G4double BeamLinePoistion[3];

  // QM1,2
  G4int    QM_OnOff; 
  G4double QM1_Typical_Field;
  G4double qm1_boundary_x[2];
  G4double qm1_boundary_y[2];
  G4double qm1_boundary_z[2];

  G4double QM2_Typical_Field;
  G4double qm2_boundary_x[2];
  G4double qm2_boundary_y[2];
  G4double qm2_boundary_z[2];

  // BM 
  G4int    BM_OnOff; 
  G4double BM_Typical_Field;
  G4double bm_boundary_x[2];
  G4double bm_boundary_y[2];
  G4double bm_boundary_z[2];

};

#endif
