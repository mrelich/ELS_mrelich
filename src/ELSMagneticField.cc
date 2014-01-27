//==============================================================================
//    ELS Magnetic Field
//      Author    : T.Shibata
//      Creation  : 2009.12.08
//    Last Update : 2010.11.11
//==============================================================================
#include "ELSMagneticField.hh"
#include "PhysicsParameters.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
//-------------------------------------------------------------------------
ELSMagnetic::~ELSMagnetic(){}
//-------------------------------------------------------------------------
ELSMagnetic::ELSMagnetic(){}
//-------------------------------------------------------------------------
void ELSMagnetic::elsmagneticfield_initialize( ELSParameters *elsparam )
{

  QM_OnOff       = 0;
  BM_OnOff       = 0;

  Bx = 0.0*tesla;
  By = 0.0*tesla;
  Bz = 0.0*tesla;
  
  c0=2.99792458;
  me=0.51099907;
  T0=40.0;
  p0=sqrt((T0+me)*(T0+me)-me*me);
  r0=20.0; //unit=cm

  QM1_Typical_Field = QM2_Typical_Field = 0.0;

  BM_Typical_Field = p0/(c0*r0)*tesla;

  // Quadrupole Magnetic Filed 1 Region
  if( elsparam->QuadrupoleManget1FiledResion(0,0) <
      elsparam->QuadrupoleManget1FiledResion(0,1) ){
     qm1_boundary_x[0] = elsparam->QuadrupoleManget1FiledResion(0,0);
     qm1_boundary_x[1] = elsparam->QuadrupoleManget1FiledResion(0,1);
  }else{
     qm1_boundary_x[0] = elsparam->QuadrupoleManget1FiledResion(0,1);
     qm1_boundary_x[1] = elsparam->QuadrupoleManget1FiledResion(0,0);
  }

  if( elsparam->QuadrupoleManget1FiledResion(1,0) <
      elsparam->QuadrupoleManget1FiledResion(1,1) ){
     qm1_boundary_y[0] = elsparam->QuadrupoleManget1FiledResion(1,0);
     qm1_boundary_y[1] = elsparam->QuadrupoleManget1FiledResion(1,1);
  }else{
     qm1_boundary_y[0] = elsparam->QuadrupoleManget1FiledResion(1,1);
     qm1_boundary_y[1] = elsparam->QuadrupoleManget1FiledResion(1,0);
  }

  if( elsparam->QuadrupoleManget1FiledResion(2,0) <
      elsparam->QuadrupoleManget1FiledResion(2,1) ){
     qm1_boundary_z[0] = elsparam->QuadrupoleManget1FiledResion(2,0);
     qm1_boundary_z[1] = elsparam->QuadrupoleManget1FiledResion(2,1);
  }else{
     qm1_boundary_z[0] = elsparam->QuadrupoleManget1FiledResion(2,1);
     qm1_boundary_z[1] = elsparam->QuadrupoleManget1FiledResion(2,0);
  }

  // Quadrupole Magnetic Filed 1 Region
  if( elsparam->QuadrupoleManget2FiledResion(0,0) <
      elsparam->QuadrupoleManget2FiledResion(0,1) ){
     qm2_boundary_x[0] = elsparam->QuadrupoleManget2FiledResion(0,0);
     qm2_boundary_x[1] = elsparam->QuadrupoleManget2FiledResion(0,1);
  }else{
     qm2_boundary_x[0] = elsparam->QuadrupoleManget2FiledResion(0,1);
     qm2_boundary_x[1] = elsparam->QuadrupoleManget2FiledResion(0,0);
  }

  if( elsparam->QuadrupoleManget2FiledResion(1,0) <
      elsparam->QuadrupoleManget2FiledResion(1,1) ){
     qm2_boundary_y[0] = elsparam->QuadrupoleManget2FiledResion(1,0);
     qm2_boundary_y[1] = elsparam->QuadrupoleManget2FiledResion(1,1);
  }else{
     qm2_boundary_y[0] = elsparam->QuadrupoleManget2FiledResion(1,1);
     qm2_boundary_y[1] = elsparam->QuadrupoleManget2FiledResion(1,0);
  }

  if( elsparam->QuadrupoleManget2FiledResion(2,0) <
      elsparam->QuadrupoleManget2FiledResion(2,1) ){
     qm2_boundary_z[0] = elsparam->QuadrupoleManget2FiledResion(2,0);
     qm2_boundary_z[1] = elsparam->QuadrupoleManget2FiledResion(2,1);
  }else{
     qm2_boundary_z[0] = elsparam->QuadrupoleManget2FiledResion(2,1);
     qm2_boundary_z[1] = elsparam->QuadrupoleManget2FiledResion(2,0);
  }

  // Bending Magnetic Filed Region
  if( elsparam->BendingMangetFiledResion(0,0) <
      elsparam->BendingMangetFiledResion(0,1) ){
     bm_boundary_x[0] = elsparam->BendingMangetFiledResion(0,0);
     bm_boundary_x[1] = elsparam->BendingMangetFiledResion(0,1);
  }else{
     bm_boundary_x[0] = elsparam->BendingMangetFiledResion(0,1);
     bm_boundary_x[1] = elsparam->BendingMangetFiledResion(0,0);
  }

  if( elsparam->BendingMangetFiledResion(1,0) <
      elsparam->BendingMangetFiledResion(1,1) ){
     bm_boundary_y[0] = elsparam->BendingMangetFiledResion(1,0);
     bm_boundary_y[1] = elsparam->BendingMangetFiledResion(1,1);
  }else{
     bm_boundary_y[0] = elsparam->BendingMangetFiledResion(1,1);
     bm_boundary_y[1] = elsparam->BendingMangetFiledResion(1,0);
  }

  if( elsparam->BendingMangetFiledResion(2,0) <
      elsparam->BendingMangetFiledResion(2,1) ){
     bm_boundary_z[0] = elsparam->BendingMangetFiledResion(2,0);
     bm_boundary_z[1] = elsparam->BendingMangetFiledResion(2,1);
  }else{
     bm_boundary_z[0] = elsparam->BendingMangetFiledResion(2,1);
     bm_boundary_z[1] = elsparam->BendingMangetFiledResion(2,0);
  }

  BeamLinePoistion[0] = elsparam->BeamLinePoistion(0);
  BeamLinePoistion[1] = elsparam->BeamLinePoistion(1);
  BeamLinePoistion[2] = elsparam->BeamLinePoistion(2);
  
}
//-------------------------------------------------------------------------
void ELSMagnetic::elsmagneticfield_onoff( int QM_onoff,
					  double QM1_field, double QM2_field,
                                          int BM_onoff, 
                                          double BM_filed )
{
  QM_OnOff = QM_onoff;
  QM1_Typical_Field = QM1_field;
  QM2_Typical_Field = QM2_field;
   
  BM_OnOff = BM_onoff;

  /* BM_filed ... Kinematic Unit = MeV */
  T0=BM_filed;
  p0=sqrt((T0+me)*(T0+me)-me*me);
  BM_Typical_Field=p0/(c0*r0)*tesla;

}
//-------------------------------------------------------------------------
void ELSMagnetic::GetFieldValue(const double Point[3],double *Bfield) const
{
  Bfield[0] = 0;
  Bfield[1] = 0;
  Bfield[2] = 0;

  if( QM_OnOff ){
      // QM 1
      if( Point[0] > qm1_boundary_x[0] && Point[0] < qm1_boundary_x[1] ){
	if( Point[1] > qm1_boundary_y[0] && Point[1] < qm1_boundary_y[1] ){
	  if( Point[2] > qm1_boundary_z[0] && Point[2] < qm1_boundary_z[1] ){
	    Bfield[0] = 0.0;
	    Bfield[1] = QM1_Typical_Field*(Point[2]-BeamLinePoistion[2])*0.001;
	    Bfield[2] = QM1_Typical_Field*(Point[1]-BeamLinePoistion[1])*0.001;
          }
	}
      }
      // QM 2
      if( Point[0] > qm2_boundary_x[0] && Point[0] < qm2_boundary_x[1] ){
	if( Point[1] > qm2_boundary_y[0] && Point[1] < qm2_boundary_y[1] ){
	  if( Point[2] > qm2_boundary_z[0] && Point[2] < qm2_boundary_z[1] ){
	    Bfield[0] = 0.0;
	    Bfield[1] = QM2_Typical_Field*(Point[2]-BeamLinePoistion[2])*0.001;
	    Bfield[2] = QM2_Typical_Field*(Point[1]-BeamLinePoistion[1])*0.001;
          }
	}
      }
  }

  if( BM_OnOff ){
      // Bending Magnet Typical Magnetic Field is 0.6 Tesla 
      if( Point[0] > bm_boundary_x[0] && Point[0] < bm_boundary_x[1] ){
	if( Point[1] > bm_boundary_y[0] && Point[1] < bm_boundary_y[1] ){
	  if( Point[2] > bm_boundary_z[0] && Point[2] < bm_boundary_z[1] ){
	        Bfield[0] = 0.0;
    	        Bfield[1] = (-1.0)*BM_Typical_Field;
		Bfield[2] = 0.0;
          }
	}
      }

  }

}
//-------------------------------------------------------------------------
