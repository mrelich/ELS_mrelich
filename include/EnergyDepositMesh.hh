//==========================================
// Energy Deposit Mesh 
//   Define the class of dE_Mesh
//  Created 2010.11.09 by T.Shibata 
//  Last Updated 2010.11.09 by T.Shibata
//  Last Updated 2011.12.12 by T.Shibata
//==========================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/Rotation.h>

#include "G4ThreeVector.hh"

#ifndef EnergyDepositMesh_h
#define EnergyDepositMesh_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class EnergyDepositMeth {

public: 
  double  energy_deposit;
  CLHEP::Hep3Vector vertex_elscoord;
  CLHEP::Hep3Vector vertex_brmcoord;
  CLHEP::Hep3Vector vertex_eathcoord;
  double global_time;
  int id_mesh;
    
  EnergyDepositMeth(){
    energy_deposit=0;
    vertex_elscoord.setX(0.);
    vertex_elscoord.setY(0.);
    vertex_elscoord.setZ(0.);

    vertex_brmcoord.setX(0.);
    vertex_brmcoord.setY(0.);
    vertex_brmcoord.setZ(0.);

    vertex_eathcoord.setX(0.);
    vertex_eathcoord.setY(0.);
    vertex_eathcoord.setZ(0.);

    global_time=0.0;

    id_mesh=-1;
  }   

};

class EnergyDepositMeth;
class dE_Mesh {

public:

   void set_delta_PathLength( double delta_pathlength );

   void set_position_of_beaminjection( CLHEP::Hep3Vector position_of_beaminjection ); 

   void set_position_of_mirror_6_BRM_coordinate( CLHEP::Hep3Vector position_of_mirror_6_BRM_coordinate ); 
   void set_angle_ELS_mirror_6( double angle_ELS_mirror_6 ); 
   void set_distance_ELS_from_mirror6( double distance_ELS_from_mirror6 );

   void set_ELS_origin_BRM_coordinate(void);

   void set_angle_North_Yaxis_BRM( double angle_north_yaxis_brm );

   void set_mesh_dV( double mesh_dv[] );
   void set_mesh_Vregion( double mesh_vregion[] );
   void set_mesh_number(void);

   double dump_delta_PathLength(void){ return delta_PathLength;  }
   double dump_segment_Pathlength(void){ return segment_Pathlength; }
   double dump_remaining_Pathlength(void){ return remaining_Pathlength; } 
   int    dump_number_of_segment(void){ return number_of_segment; }

   double dump_segment_EnergyDeposit(void){ return segment_EnergyDeposit; }
   double dumP_remaining_EnergyDeposit(void){ return remaining_EnergyDeposit; }

   typedef vector<EnergyDepositMeth> CreatedEnergyDepositMeth;

   int energy_deposit_mesh_generator( int process_id,
                                      double energy_deposit,
 				      CLHEP::Hep3Vector particle_preposition,
				      CLHEP::Hep3Vector particle_postposition,
                                      CreatedEnergyDepositMeth &created_energy_deposit_meth );

  // added for Cerenkov Study
  CLHEP::Hep3Vector PositionTransportatinoELStoEARTCH( CLHEP::Hep3Vector inputposition );
    
  // Consttructor    
  dE_Mesh();

  // DeConstructor
  ~dE_Mesh();

protected:

private:

  double delta_PathLength;
  
  CLHEP::Hep3Vector Position_of_beaminjection;
  CLHEP::Hep3Vector Position_of_mirror_6_BRM_coordinate;
  double Angle_ELS_mirror_6;
  double Distance_ELS_from_mirror6; 
  double Angle_North_Yaxis_BRM;

  double Mesh_dV[3];
  double Mesh_Vregion[3];

  int    Number_of_Mesh[3];
  int    ID_Mesh;

  double EnergyDeposit;
  double PathLength;
  double segment_Pathlength;
  double remaining_Pathlength;
    
  int    number_of_segment;

  double segment_EnergyDeposit;
  double remaining_EnergyDeposit;

  CLHEP::Hep3Vector Particle_PrePosition;
  CLHEP::Hep3Vector Particle_PostPosition;
  CLHEP::Hep3Vector Particle_Direction;
  CLHEP::Hep3Vector Particle_SegmentPosition;

  CLHEP::Hep3Vector Particle_SegmentPosition_BRM_coordinate;
  CLHEP::Hep3Vector Particle_SegmentPosition_EATH_coordinate;

  CLHEP::Hep3Vector  Position_of_ELS_injection_BRM_coordinate;
  CLHEP::Hep3Vector  Position_of_ELS_BRM_coordinate;
  
  CLHEP::Hep3Vector Coordinate_Transportation_ELS2BRM( CLHEP::Hep3Vector input_positioin );
  CLHEP::Hep3Vector Coordinate_Transportation_BRM2EARTH( CLHEP::Hep3Vector input_positioin );
  int Mesh_number( CLHEP::Hep3Vector input_positioin );

};
#endif


