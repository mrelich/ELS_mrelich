//===================================================================
//    SetUp Parameters 
//      Author    : T.Shibata
//
//      Creation  : 2009.04.01
//    Last Update : 2010.12.03
//    Last Update : 2011.02.03
//    Last Update : 2011.05.15
//    Last Update : 2011.12.12
//    Last Update : 2011.12.13
//    Last Update : 2011.12.14
//    Last Update : 2012.01.07
//    Last Update : 2012.01.20  # added to Geomagnetic Filed
//    Last Update : 2012.01.25
//    Last Update : 2012.01.31
//    Last Update : 2012.09.24
//    Last Update : 2012.10.10
//===================================================================
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <math.h>
#include <string.h>
#include <signal.h>

#ifndef ReadFile_h
#define ReadFile_h
//----------------------------------------------------------------------
using namespace std;
//----------------------------------------------------------------------
class Plist{

public:
  //Contstruction
  Plist();
  //Deconstruction
  ~Plist();

  void readfile(char *setupfile); 

  // Public Detector Parameter //
 
  // World Geometory 
  double SizeOfWorld( int i ){ return size_of_world[i]; }

  double ELSSiteGroundHeight( int i ){ return elssite_ground_height[i]; }

  // Public Parameter //

  /* Electron Kinematic Parameters */
  int beamenergy_mode( void ){ return beam_energy_mode; }

  int beam_particle(void){ return beam_particle_id; }

  /* Electron Beam Waveform File */
  int ElectronBeam_Waveform_flag(void){ return electronbeam_waveform_flag; }  
  char *ElectronBeam_Waveform_File(void){ return electronbeam_waveform_file; }
  double ElectronBeam_Waveform_dt(void){ return electronbeam_waveform_dt; }
  int ElectronBeam_Waveform_new_ntbin(void){ return electronbeam_waveform_new_ntbin; }
    
  /* OUTPUT FILE Selectrion */
  int Outoputfileflag_detectorplane( void ){ return outoputfileflag_detectorplane; }
  int Outoputfileflag_secondparticle( void ){ return outoputfileflag_secondparticle; }
  int Outoputfileflag_dedx( void ){ return outoputfileflag_dedx; }
  int Outoputfileflag_faradaycup( void ){ return faradaycup_flag; }  
  int Outoputfileflag_faradaycup4( void ){ return faradaycup4_flag;}
  int Outoputfileflag_energydeposit( void ){ return outoputfileflag_energydeposit; }  

  /* Beam Energy */ 
  double primaryenergy( int i ){ return primary_energy[i]; }

  double primaryenergy_crystalball_mean( int i ){ return crystalball_mean[i]; }
  double primaryenergy_crystalball_sigma( int i ){ return crystalball_sigma[i]; }
  double primaryenergy_crystalball_alpha( int i ){ return crystalball_alpha[i]; }
  double primaryenergy_crystalball_n( int i ){ return crystalball_n[i]; }
  double primaryenergy_crystalball_xmin( int i ){ return crystalball_xmin[i]; }
  double primaryenergy_crystalball_xmax( int i ){ return crystalball_xmax[i]; }

  int beam_emmitance_mode(void){ return beamemmitance_mode; }

  // Beam injection position ( Center position )
  double injection_positionx( int i ){ return injection_position_x[i]; }
  double injection_positiony( int i ){ return injection_position_y[i]; }
  double injection_positionz( int i ){ return injection_position_z[i]; }
  //                         ( Sigma = Spread ) 
  double injection_positionsigma( int i ){ return injection_position_sigma[i]; }

  double injection_directionx( int i ){ return injection_direction_x[i]; }
  double injection_directiony( int i ){ return injection_direction_y[i]; }
  double injection_directionz( int i ){ return injection_direction_z[i]; }

  double injection_transverse_phi0x( int i ){ return injection_transverse_phi0_x[i]; }
  double injection_transverse_phi0y( int i ){ return injection_transverse_phi0_y[i]; }

  // shift of the center of beam position 
  double injection_position_shiftx( int i ){ return injection_position_shift_x[i]; }
  double injection_position_sigmax( int i ){ return injection_position_sigma_x[i]; }
  double injection_position_shifty( int i ){ return injection_position_shift_y[i]; }
  double injection_position_sigmay( int i ){ return injection_position_sigma_y[i]; }

  // Beam Emmitance 
  // i=0 : horizantal
  // Beam Emmitance 
  // i=0 : horizantal 
  //  =1 : vertival 
  //     g: gamma
  //     b: beta
  //     a: alpha
  //     e: emmitance
  double injection_beamemmitance_g( int i, int j ){ return injection_beam_emmitance_g[i][j]; }
  double injection_beamemmitance_b( int i, int j ){ return injection_beam_emmitance_b[i][j]; }
  double injection_beamemmitance_a( int i, int j ){ return injection_beam_emmitance_a[i][j]; }
  double injection_beamemmitance_e( int i, int j ){ return injection_beam_emmitance_e[i][j]; }

  // Atmosphere Condition
  // Air Compositin 
  double Air_composition_nitrogen(int i){ return air_composition_N2[i]; }  
  double Air_composition_oxygen(int i){ return air_composition_O2[i]; }  
  double Air_composition_argon(int i){ return air_composition_Ar[i]; }  
  double Air_composition_carbon_dioxide(int i){ return air_composition_CO2[i]; }  
  // Air Condition
  double Air_temperature(int i){ return air_temperature[i]; }
  double Air_pressure(int i){ return air_pressure[i]; }
  double Air_relative_humidity(int i){ return air_relative_humidity[i]; }

  // Geomagnetic Field @ BRM ( ELS )
  int Geomagnetic_onoff(void){ return geomagnetic_onoff; }
  
  double Geomagnetic_theta(int i){ return geomagnetic_theta[i]; }
  double Geomagnetic_declination(int i){ return geomagnetic_declination[i]; } 
  double Geomagnetic_horizontal(int i){ return geomagnetic_horizontal[i]; }
  double Geomagnetic_vertical(int i){ return geomagnetic_vertical[i]; }

  // Cu Collimter Diameter
  double CuCollimeterDiameter(void){ return cu_collimeter_d; }

  // Slit Width
  double SlitWidth(void){ return slit_width; }

  // -- Beam Dump Flag
  // Faraday Cup 
  //  = 1 : Faraday Cup 1 
  //  = 4 : Faraday Cup 4  
  //  = 5 : Faraday Cup 5
  //  = other : No Faraday Cup
  int Faradaycup_flag(void){ return faradaycup_flag; }

  // made by B.K.Shin-san, added in 2013.04.23
  double Faradaycup5_height (void ){ return faradaycup5_height ;  }

  // Fraday Cup 4 
  // Not be used any more, in 2013.04.23
  int Faradaycup4_flag(void){ return faradaycup4_flag; }

  double R_Faradaycup4_front1(void){ return r_faradaycup4_front1; }
  double Length_Faradaycup4_front1(void){ return length_faradaycup4_front1; }
  double Positionz_Faradaycup4_font1(void){ return positionz_faradaycup4_front1; }
  double R_Faradaycup4_front2(void){ return r_faradaycup4_front2; }
  double Length_Faradaycup4_front2(void){ return length_faradaycup4_front2; }
  double Positionz_Faradaycup4_font2(void){ return positionz_faradaycup4_front2; }

  double R_Faradaycup4_copper(void){ return r_faradaycup4_copper; }
  double Length_Faradaycup4_copper(void){ return length_faradaycup4_copper; }
  double Positionz_Faradaycup4_copper(void){ return positionz_faradaycup4_copper; }

  double Positionz_alongbeamline_Faradaycup4(void){ return position_alongbeamline_faradaycup4; }
  double Positionz_transverse_Faradaycup4(void){ return position_transverse_faradaycup4; }

  // Screen Monitor 3
  int Screen_monitor3_flag(void){ return screen_monitor3_flag; }
  double Positionz_ScreenMonitor3(void){ return positionz_screen_monitor3; }
  
  // Beam Attenuator
  int Beam_attenuator_flag(void){ return beam_attenuator_flag; }
  int Beam_attenuator_material(void){ return beam_attenuator_material; }   
  double Length_Beam_Attenuator(void){ return length_beam_attenuator; }
  
  // Pb Collimator added in 2012.09.24
  int Pb_collimator_flag(void){ return pb_collimator_flag; }
  double Dinner_Pb_Collimator(void){ return dinner_pb_collimator; }
  double Douter_Pb_Collimator(void){ return douter_pb_collimator; }
  double Length_Pb_Collimator(void){ return length_pb_collimator; }

  // Virtual Test Chamber added in 2011.12.13
  int Virtual_chamber_flag(void){ return virtual_chamber_flag; }

  double R_virtual_chamber_injectionhole(void){ return r_virtual_chamber_injectionhole; }  
  double R_virtual_chamber_innerhole(void){ return r_virtual_chamber_innerhole; }
  double Length_virtual_chamber(void){ return length_virtual_chamber; }
  double Gap_virtual_chamber(void){ return gap_virtual_chamber; }
  double Positionz_virtual_chamber(void){ return positionz_virtual_chamber; }

  double VC_temperature(void){ return vc_temperature; }
  double VC_pressure(void){ return vc_pressure; }

  // -- QM Magnetic Field Flag
  int QM_Magnetic_field_flag(void){ return qm_magnetic_field_flag; }
  double QM1_Magnetic_field(void){ return qm1_magnetic_field; }  
  double QM2_Magnetic_field(void){ return qm2_magnetic_field; }  

  // -- BM Magnetic Field Flag
  int BM_Magnetic_field_flag(void){ return bm_magnetic_field_flag; }
  double BM_Magnetic_field(void){ return bm_magnetic_field; }

  // -- Alignment Parameters
  double CuCollimeter_shiftz(void){ return cu_collimeter_shiftz; }
  double Vertical_beamline_shift(int i){ return vertical_beamline_shift[i]; }
  double BMYork_Coil_shift(int i){ return bmyork_coil_shift[i]; }
  double Slit_Collimeter_shift(int i){ return slit_collimeter_shift[i]; } 
  double Faradaycup_shift(int i){ return faradaycup_shift[i]; }
  double Coverbox_Hole_shift(int i){ return coverbox_hole_shift[i]; }

  // -- Physics Proccess Flag
  int Photo_nuclear_process_flag(void){ return photo_nuclear_process_flag; }  

  // -- Cerenkov Process Flag
  int Cerenkov_process_flag(void){ return cerenkov_process_flag; }

  double Cutlenght_particle(void){ return cutlenght_particle; };

  // -- Energy Deposit Mesh 
  double delta_PathLength(void){ return delta_pathlength;}

  double Position_of_mirror_6_BRM_coordinate(int i){ return position_of_mirror_6_BRM_coordinate[i]; }
  double Angle_ELS_mirror_6(void){ return angle_ELS_mirror_6; }
  double Distance_ELS_from_mirror6(void){ return distance_ELS_from_mirror6; }

  double Angle_North_Yaxis_BRM(void){return angle_north_yaxis_brm;}

  double Mesh_dX(void){ return mesh_dx; }
  double Mesh_dY(void){ return mesh_dy; }
  double Mesh_dZ(void){ return mesh_dz; }
  double Mesh_Xregion(void){ return mesh_xregion; }
  double Mesh_Yregion(void){ return mesh_yregion; }
  double Mesh_Zregion(void){ return mesh_zregion; }

private:

  int skipspaceline( char *c );
  int skipcomment( char *c );
  int read_parameterID( char* c );

  void read_parameter( int paramID, int paramN, char* c );
  void setup_parameters(void);
  void dump_parameter(void);
  void read_line( char *c );

  // Primary Electron //
  int beam_energy_mode;

  int beam_particle_id;

  int electronbeam_waveform_flag;
  char electronbeam_waveform_file[1024];
  double electronbeam_waveform_dt;
  int electronbeam_waveform_new_ntbin;

  int outoputfileflag_detectorplane;
  int outoputfileflag_secondparticle;
  int outoputfileflag_dedx;
  int outoputfileflag_energydeposit;

  double primary_energy[4];

  double crystalball_mean[4];
  double crystalball_sigma[4];
  double crystalball_alpha[4];
  double crystalball_n[4];
  double crystalball_xmin[4];
  double crystalball_xmax[4];

  int beamemmitance_mode;

  double injection_position_x[4];
  double injection_position_y[4];
  double injection_position_z[4];

  double injection_position_sigma[4];

  double injection_direction_x[4];
  double injection_direction_y[4];
  double injection_direction_z[4];

  double injection_transverse_phi0_x[4];
  double injection_transverse_phi0_y[4];

  double injection_position_shift_x[4];
  double injection_position_sigma_x[4];
  double injection_position_shift_y[4];
  double injection_position_sigma_y[4];

  double injection_beam_emmitance_g[2][4];
  double injection_beam_emmitance_b[2][4];
  double injection_beam_emmitance_a[2][4];
  double injection_beam_emmitance_e[2][4];

  // Private Detector Parameter //

  // World Geometory 
  double size_of_world[3];

  double elssite_ground_height[4];
  
  // Atmosphere Condition
  double air_composition_N2[4]; 
  double air_composition_O2[4];
  double air_composition_Ar[4];
  double air_composition_CO2[4];
  
  double air_temperature[4];
  double air_pressure[4];
  double air_relative_humidity[4];

  // Geomagnetic Field @ BRM ( ELS )                                     
  double geomagnetic_onoff;
  double geomagnetic_theta[4];
  double geomagnetic_declination[4];
  double geomagnetic_horizontal[4];
  double geomagnetic_vertical[4]; 
 
  // Cu Collimter Diameter
  double cu_collimeter_d; 

  // Slit width
  double slit_width;

  // -- Beam Dump Flag ... and geometory parameters 
  int faradaycup_flag;

  // Faraday Cup 4
  int faradaycup4_flag;

  double r_faradaycup4_front1;
  double length_faradaycup4_front1;
  double positionz_faradaycup4_front1;

  double r_faradaycup4_front2;
  double length_faradaycup4_front2;
  double positionz_faradaycup4_front2;

  double r_faradaycup4_copper;
  double length_faradaycup4_copper;
  double positionz_faradaycup4_copper;

  double position_alongbeamline_faradaycup4;
  double position_transverse_faradaycup4;

  // Faraday Cup 5 
  double faradaycup5_height;

  //  Screen Monitor 3
  int screen_monitor3_flag;
  double positionz_screen_monitor3;

  // Beam Attenuator
  int beam_attenuator_flag;
  int beam_attenuator_material;
  double length_beam_attenuator;

  // Pb Collimator added in 2012.09.24
  int pb_collimator_flag;
  double dinner_pb_collimator;
  double douter_pb_collimator;
  double length_pb_collimator;

  // Virtual Test Chamber added in 2011.12.13                
  int virtual_chamber_flag;
  double r_virtual_chamber_injectionhole; 
  double r_virtual_chamber_innerhole;
  double length_virtual_chamber;
  double gap_virtual_chamber;
  double positionz_virtual_chamber;

  double vc_temperature;
  double vc_pressure;

  // -- Magnetic Field Flag
  int qm_magnetic_field_flag;
  double qm1_magnetic_field;
  double qm2_magnetic_field;

  int bm_magnetic_field_flag;
  double bm_magnetic_field;

  // -- Alignment Parameters
  double cu_collimeter_shiftz; 
  double vertical_beamline_shift[2];
  double bmyork_coil_shift[3];
  double slit_collimeter_shift[2];
  double faradaycup_shift[3];
  double coverbox_hole_shift[2];

  // -- Physics Proccess Flag
  int photo_nuclear_process_flag; 
  int cerenkov_process_flag;
  
  double cutlenght_particle;

  // -- Energy Deposit Mesh                                                     
  double delta_pathlength;

  double position_of_mirror_6_BRM_coordinate[3];
  double angle_ELS_mirror_6;
  double distance_ELS_from_mirror6;

  double angle_north_yaxis_brm;

  double mesh_dx;
  double mesh_dy;
  double mesh_dz;
  double mesh_xregion;
  double mesh_yregion;
  double mesh_zregion;

};
//----------------------------------------------------------------------
#endif
