//==========================================
// Energy Deposit Mesh                      
//   Define the class of dE_Mesh
//  Created 2010.11.09 by T.Shibata
//  Modified 2010.12.06 by T.Shibata
//  Modified(bug fixed) 2011.01.26 by T.Shibata
//  Modified 2011.12.12 by T.Shibata
//  Modified 2011.12.24 by T.Shibata
//==========================================      
#include "PhysicsParameters.hh"
#include "EnergyDepositMesh.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
dE_Mesh::dE_Mesh()
{

  delta_PathLength     = 0.2;

  Position_of_beaminjection.setX(11.4862);
  Position_of_beaminjection.setY(-1.7182);
  Position_of_beaminjection.setZ(0.0);
  
  Position_of_mirror_6_BRM_coordinate.setX(-2.166254023);
  Position_of_mirror_6_BRM_coordinate.setY(12.654034);
  Position_of_mirror_6_BRM_coordinate.setZ(0.0);

  Angle_ELS_mirror_6 = 100.0;
  Distance_ELS_from_mirror6 = 9.0;

  Angle_North_Yaxis_BRM=69.39492;

  Mesh_dV[0] = 0.1;
  Mesh_dV[1] = 0.1;
  Mesh_dV[2] = 0.1;

  Mesh_Vregion[0] = 200.0;
  Mesh_Vregion[1] = 200.0;
  Mesh_Vregion[2] = 200.0;
 
  Number_of_Mesh[0] = Mesh_Vregion[0]/Mesh_dV[0];
  Number_of_Mesh[1] = Mesh_Vregion[1]/Mesh_dV[1];
  Number_of_Mesh[2] = Mesh_Vregion[2]/Mesh_dV[2];

  ID_Mesh = -1;

  EnergyDeposit        = 0.;
  PathLength           = 0.;
  segment_Pathlength   = 0.;
  remaining_Pathlength = 0.;
   
  number_of_segment = 0;

  segment_EnergyDeposit =0.;
  remaining_EnergyDeposit =0.; 

  Particle_PrePosition.setX(0);
    Particle_PrePosition.setY(0);
      Particle_PrePosition.setZ(0);

  Particle_PostPosition.setX(0);
    Particle_PostPosition.setY(0);
      Particle_PostPosition.setZ(0);

  Particle_Direction.setX(0);
    Particle_Direction.setY(0);
      Particle_Direction.setZ(0);

  Particle_SegmentPosition.setX(0);
    Particle_SegmentPosition.setY(0);
      Particle_SegmentPosition.setZ(0);

  Particle_SegmentPosition_BRM_coordinate.setX(0);
    Particle_SegmentPosition_BRM_coordinate.setY(0);
      Particle_SegmentPosition_BRM_coordinate.setZ(0);

  Position_of_ELS_BRM_coordinate.setX(0);
    Position_of_ELS_BRM_coordinate.setY(0);
     Position_of_ELS_BRM_coordinate.setZ(0);

}
//-------------------------------------------------------------------------------------------
dE_Mesh::~dE_Mesh(){ 

}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_delta_PathLength( double delta_pathlength )
{
  delta_PathLength = delta_pathlength;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_position_of_beaminjection( CLHEP::Hep3Vector position_of_beaminjection )
{
  Position_of_beaminjection = position_of_beaminjection;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_position_of_mirror_6_BRM_coordinate( CLHEP::Hep3Vector position_of_mirror_6_BRM_coordinate ){
  Position_of_mirror_6_BRM_coordinate = position_of_mirror_6_BRM_coordinate;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_angle_ELS_mirror_6( double angle_ELS_mirror_6 )
{
  Angle_ELS_mirror_6 = angle_ELS_mirror_6;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_distance_ELS_from_mirror6( double distance_ELS_from_mirror6 )
{
  Distance_ELS_from_mirror6 = distance_ELS_from_mirror6;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_angle_North_Yaxis_BRM( double angle_north_yaxis_brm )
{
  Angle_North_Yaxis_BRM = angle_north_yaxis_brm;
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_mesh_dV( double mesh_dv[] )
{
  for( int i=0; i<3; i++){ Mesh_dV[i]=mesh_dv[i]; }
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_mesh_Vregion( double mesh_vregion[] )
{
  for( int i=0; i<3; i++){ Mesh_Vregion[i] = mesh_vregion[i]; }
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_mesh_number( void )
{
  for( int i=0; i<3; i++){ Number_of_Mesh[i] = Mesh_Vregion[i]/Mesh_dV[i]; }
}
//-------------------------------------------------------------------------------------------
void dE_Mesh::set_ELS_origin_BRM_coordinate(void)
{
  
  double phi = Position_of_beaminjection.phi();
  if( phi < 0 ) phi*=(-1);
  double r   = Position_of_beaminjection.perp();

  double x_els = Position_of_mirror_6_BRM_coordinate.x() 
                  + Distance_ELS_from_mirror6*sin(Angle_ELS_mirror_6*RADI);
  double y_els = Position_of_mirror_6_BRM_coordinate.y() 
                  + Distance_ELS_from_mirror6*cos(Angle_ELS_mirror_6*RADI);
  double z_els = Position_of_mirror_6_BRM_coordinate.z();

  Position_of_ELS_injection_BRM_coordinate.setX(x_els);
  Position_of_ELS_injection_BRM_coordinate.setY(y_els);
  Position_of_ELS_injection_BRM_coordinate.setZ(Position_of_beaminjection.z());

  double delta_x = r * sin( phi + Angle_ELS_mirror_6*RADI );
  double delta_y = r * cos( phi + Angle_ELS_mirror_6*RADI );
  double delta_z = 0.; 

  Position_of_ELS_BRM_coordinate.setX( x_els - delta_x );
  Position_of_ELS_BRM_coordinate.setY( y_els - delta_y );
  Position_of_ELS_BRM_coordinate.setZ( z_els - delta_z );

}
//-------------------------------------------------------------------------------------------
int dE_Mesh::energy_deposit_mesh_generator(  int process_id,
                                             double energy_deposit,  // unit=MeV
		              		                 CLHEP::Hep3Vector particle_preposition,  // unit=m
				                             CLHEP::Hep3Vector particle_postposition, // unit=m
				                             CreatedEnergyDepositMeth &created_energy_deposit_meth )
{  
  // Input the Energy Deposit 
  EnergyDeposit = energy_deposit;
  if( EnergyDeposit < 0 ){ 
    fprintf(stderr, "Energy deposit is less than 0... Error and exit!!\n");
    //return -1;
    return 0;
  } 
  
  // Input the PrePosition and PostPosition
  // and calculate Particle direction.  This direction is not momentum
  Particle_PrePosition  = particle_preposition;
  Particle_PostPosition = particle_postposition;
 
  Particle_Direction = (Particle_PostPosition-Particle_PrePosition).unit(); // Direction Unit Vector 
  
  // Calculate the Pahtlength 
  PathLength = (Particle_PostPosition-Particle_PrePosition).mag();
  if( PathLength < 0){
      fprintf(stderr, "pathlengtht is less than 0... Error and exit!!\n");
      //return -1;
      return 0;
  }

  // Calculate the number of segment energy deposit region.
  if( delta_PathLength <= 0 ){
    fprintf(stderr, "delta pathlengtht is less than 0... Error and exit!!\n");
    //return -1;
    return 0;
  }

  number_of_segment = int(PathLength/delta_PathLength);   
  remaining_Pathlength = PathLength - double(number_of_segment)*delta_PathLength;

  segment_EnergyDeposit   = EnergyDeposit/(PathLength/delta_PathLength);
  remaining_EnergyDeposit = EnergyDeposit-double(number_of_segment)*segment_EnergyDeposit;  

  if( process_id != 1 ){
    // process_id = 1 is Ionization by primary and secondary electrons
    number_of_segment = 0;
    remaining_Pathlength = 0.;
   
    segment_EnergyDeposit   = EnergyDeposit;
    remaining_EnergyDeposit = 0.;
  }

  // EnergyDepositMeth is created.
  EnergyDepositMeth h;

  // Calculate the Segment Energy Deposit
  for( int i=0; i<number_of_segment+1; i++ ){

    if( process_id == 1 ){

        if( i < number_of_segment ){
            Particle_SegmentPosition = ((double(i)+0.5)*delta_PathLength)*Particle_Direction 
                                       + Particle_PrePosition;
        }else if( i == number_of_segment ){
            segment_EnergyDeposit = remaining_EnergyDeposit;      
   	        Particle_SegmentPosition = (number_of_segment*delta_PathLength+remaining_Pathlength*0.5)*Particle_Direction
   	                                  + Particle_PrePosition;
        }
      //Particle_SegmentPosition = Particle_PostPosition; // test in 2012.06.07 

    }else{
          Particle_SegmentPosition = Particle_PostPosition;
    }
  
        Particle_SegmentPosition_BRM_coordinate = Coordinate_Transportation_ELS2BRM( Particle_SegmentPosition );

        Particle_SegmentPosition_EATH_coordinate = 
                                   Coordinate_Transportation_BRM2EARTH(Particle_SegmentPosition_BRM_coordinate); 

        ID_Mesh = Mesh_number(Particle_SegmentPosition_BRM_coordinate); 

        h.energy_deposit   = segment_EnergyDeposit;
        h.vertex_elscoord  = Particle_SegmentPosition;
        h.vertex_brmcoord  = Particle_SegmentPosition_BRM_coordinate;
        h.vertex_eathcoord = Particle_SegmentPosition_EATH_coordinate;
        h.id_mesh          = ID_Mesh;
	
        created_energy_deposit_meth.push_back(h);
  } 
  return 1;
   
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector dE_Mesh::Coordinate_Transportation_ELS2BRM( CLHEP::Hep3Vector input_positioin )
{

  CLHEP::Hep3Vector output_positioin; 

  CLHEP::Hep3Vector position_els_coordinate;   // unit = m
  position_els_coordinate.setX( input_positioin.x() - Position_of_beaminjection.x() );
  position_els_coordinate.setY( input_positioin.y() - Position_of_beaminjection.y() );
  position_els_coordinate.setZ( input_positioin.z() );

  double ANGLE=RADI*(90.0-Angle_ELS_mirror_6);
  CLHEP::Hep3Vector position_els_coordinate_2;
  position_els_coordinate_2.setX( position_els_coordinate.x()*cos( ANGLE ) 
				                  - position_els_coordinate.y()*sin( ANGLE ) );
  position_els_coordinate_2.setY( position_els_coordinate.x()*sin( ANGLE ) 
				                  + position_els_coordinate.y()*cos( ANGLE ) );
  position_els_coordinate_2.setZ(  position_els_coordinate.z() );

  ANGLE=RADI*Angle_ELS_mirror_6;

  output_positioin.setX( position_els_coordinate_2.x() 
                          + Position_of_mirror_6_BRM_coordinate.x() 
			              + Distance_ELS_from_mirror6*sin( ANGLE ) );
  output_positioin.setY( position_els_coordinate_2.y() 
                          + Position_of_mirror_6_BRM_coordinate.y() 
		   	              + Distance_ELS_from_mirror6*cos( ANGLE ) );
  output_positioin.setZ(  position_els_coordinate_2.z() );

  return output_positioin;
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector dE_Mesh::Coordinate_Transportation_BRM2EARTH( CLHEP::Hep3Vector input_positioin )
{
  /*
    Earth Coordinate 
    X = east
    Y = north
  */

  CLHEP::Hep3Vector output_positioin;

  double ANGLE=RADI*Angle_North_Yaxis_BRM;
  
  output_positioin.setX( input_positioin.x()*cos(ANGLE) - input_positioin.y()*sin(ANGLE) );
  output_positioin.setY( input_positioin.x()*sin(ANGLE) + input_positioin.y()*cos(ANGLE) );
  output_positioin.setZ(input_positioin.z());

  return output_positioin;
}
//-------------------------------------------------------------------------------------------
int dE_Mesh::Mesh_number( CLHEP::Hep3Vector input_positioin )
{

  int ix = int( (input_positioin.x()-Position_of_ELS_injection_BRM_coordinate.x()
	  	 -(-1)*Mesh_Vregion[0])/Mesh_dV[0] );
  int iy = int( (input_positioin.y()-Position_of_ELS_injection_BRM_coordinate.y()
	  	 -(-1)*Mesh_Vregion[1])/Mesh_dV[1] );
  int iz = int( (input_positioin.z()-Position_of_ELS_injection_BRM_coordinate.z()
                 -(-1)*Mesh_Vregion[2])/Mesh_dV[2] ); 

  return ix + iy*Number_of_Mesh[0] + iz*Number_of_Mesh[0]*Number_of_Mesh[1];
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector dE_Mesh::PositionTransportatinoELStoEARTCH( CLHEP::Hep3Vector inputposition ) // = m 
{  
  CLHEP::Hep3Vector brm_coord_position=Coordinate_Transportation_ELS2BRM(inputposition);
  CLHEP::Hep3Vector earth_coord_position=Coordinate_Transportation_BRM2EARTH(brm_coord_position); 
  return earth_coord_position;   
}
//==================================================================================
