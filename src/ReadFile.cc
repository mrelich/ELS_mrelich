//===================================================================
//    SetUp Parameters
//      Author    : T.Shibata
//
//      Creation  : 2006.10.18
//    Last Update : 2010.11.09
//    Last Update : 2011.05.15
//    Last Update : 2011.12.12
//    Last Update : 2011.12.13
//    Last Update : 2011.12.14
//    Last Update : 2012.01.07
//    Last Update : 2012.01.20  # added to Geomagnetic Filed
//    Last Update : 2012.01.31
//    Last Update : 2012.09.24
//    Last Update : 2012.10.10
//===================================================================
#include "ReadFile.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
Plist::Plist(){

  // Set Up Parameter Initialized  --------------

  elssite_ground_height[0] = 0.;
  elssite_ground_height[1] = 0.;
  elssite_ground_height[2] = 0.;
  elssite_ground_height[3] = 1;

  outoputfileflag_detectorplane     = 0;
  outoputfileflag_secondparticle    = 0;
  outoputfileflag_dedx              = 0; 
  outoputfileflag_energydeposit     = 0; 
  
  beam_energy_mode  = 1;

  beam_particle_id = -1;
 
  electronbeam_waveform_flag = 0;
  electronbeam_waveform_dt = 1.;
  electronbeam_waveform_new_ntbin = 10000;

  primary_energy[0] = 40.0;
  primary_energy[1] = 0.01;
  primary_energy[2] = 0.01;
  primary_energy[3] = 1;

  crystalball_mean[0]= 40.0;  
  crystalball_mean[1]= 0.003;
  crystalball_mean[2]= 0.003;
  crystalball_mean[3]= 1;

  crystalball_sigma[0]= 0.2409;
  crystalball_sigma[1]= 0.003;
  crystalball_sigma[2]= 0.003;
  crystalball_sigma[3]= 1;

  crystalball_alpha[0]= 1.8209;
  crystalball_alpha[1]= 0.09;
  crystalball_alpha[2]= 0.09;
  crystalball_alpha[3]= 1;

  crystalball_n[0] = 0.67931;
  crystalball_n[1] = 0.16;
  crystalball_n[2] = 0.16;
  crystalball_n[3] = 1;

  crystalball_xmin[0] = 34.0;
  crystalball_xmin[1] = 0.0;
  crystalball_xmin[2] = 0.0;
  crystalball_xmin[3] = 1;
   
  crystalball_xmax[0] = 41.0;
  crystalball_xmax[1] = 0.0;
  crystalball_xmax[2] = 0.0;
  crystalball_xmax[3] = 1;

  beamemmitance_mode  = 0;

  injection_position_x[0] = 5597.0;
  injection_position_x[1] = 0.0;
  injection_position_x[2] = 0.0;
  injection_position_x[3] = 1;

  injection_position_y[0] = -1718.2;
  injection_position_y[1] = 0.0;
  injection_position_y[2] = 0.0;
  injection_position_y[3] = 1;

  injection_position_z[0] = 1300.3;
  injection_position_z[1] = 0.0;
  injection_position_z[2] = 0.0;
  injection_position_z[3] = 1;

  injection_position_sigma[0] = 0.0; //2.5;
  injection_position_sigma[1] = 0.0;
  injection_position_sigma[2] = 0.0;
  injection_position_sigma[3] = 1;

  injection_direction_x[0] = 1.0;
  injection_direction_x[1] = 0.0;
  injection_direction_x[2] = 0.0;
  injection_direction_x[3] = 1;
 
  injection_direction_y[0] = 0.0;
  injection_direction_y[1] = 0.0;
  injection_direction_y[2] = 0.0;
  injection_direction_y[3] = 1;

  injection_direction_z[0] = 0.0;
  injection_direction_z[1] = 0.0;
  injection_direction_z[2] = 0.0;
  injection_direction_z[3] = 1;

  injection_transverse_phi0_x[0]=0.0;
  injection_transverse_phi0_x[1]=0.0;
  injection_transverse_phi0_x[2]=0.0;
  injection_transverse_phi0_x[3]=1;

  injection_transverse_phi0_y[0]=0.0;
  injection_transverse_phi0_y[1]=0.0;
  injection_transverse_phi0_y[2]=0.0;
  injection_transverse_phi0_y[3]=1;

  injection_position_shift_x[0] = 0.0;
  injection_position_shift_x[1] = 0.0;
  injection_position_shift_x[2] = 0.0;
  injection_position_shift_x[3] = 1;

  injection_position_sigma_x[0] = 0.0;
  injection_position_sigma_x[1] = 0.0;
  injection_position_sigma_x[2] = 0.0;
  injection_position_sigma_x[3] = 1;

  injection_position_shift_y[0] = 0.0;
  injection_position_shift_y[1] = 0.0;
  injection_position_shift_y[2] = 0.0;
  injection_position_shift_y[3] = 1;

  injection_position_sigma_y[0] = 0.0;
  injection_position_sigma_y[1] = 0.0;
  injection_position_sigma_y[2] = 0.0;
  injection_position_sigma_y[3] = 1;

  injection_beam_emmitance_g[0][0] = 19.42179;
  injection_beam_emmitance_b[0][0] = 0.386771;
  injection_beam_emmitance_a[0][0] = -2.55182;
  injection_beam_emmitance_e[0][0] = 1.12408E-7;

  injection_beam_emmitance_g[1][0] = 36.47994;
  injection_beam_emmitance_b[1][0] = 0.36142;
  injection_beam_emmitance_a[1][0] = -3.4906;
  injection_beam_emmitance_e[1][0] = 2.13585E-7;

  injection_beam_emmitance_g[0][1]=injection_beam_emmitance_g[0][2]=0.0; 
  injection_beam_emmitance_b[0][1]=injection_beam_emmitance_b[0][2]=0.0;
  injection_beam_emmitance_a[0][1]=injection_beam_emmitance_a[0][2]=0.0;
  injection_beam_emmitance_e[0][1]=injection_beam_emmitance_e[0][2]=0.0;

  injection_beam_emmitance_g[0][3]=1.;
  injection_beam_emmitance_b[0][3]=1.;
  injection_beam_emmitance_a[0][3]=1.;
  injection_beam_emmitance_e[0][3]=1.;

  injection_beam_emmitance_g[1][1]=injection_beam_emmitance_g[1][2]=0.0; 
  injection_beam_emmitance_b[1][1]=injection_beam_emmitance_b[1][2]=0.0;
  injection_beam_emmitance_a[1][1]=injection_beam_emmitance_a[1][2]=0.0;
  injection_beam_emmitance_e[1][1]=injection_beam_emmitance_e[1][2]=0.0;

  injection_beam_emmitance_g[1][3]=1.;
  injection_beam_emmitance_b[1][3]=1.;
  injection_beam_emmitance_a[1][3]=1.;
  injection_beam_emmitance_e[1][3]=1.;

  size_of_world[0] = 100.0;
  size_of_world[1] = 100.0;
  size_of_world[2] = 100.0;

  air_composition_N2[0]=78.084;
  air_composition_O2[0]=20.946;
  air_composition_Ar[0]=0.930;
  air_composition_CO2[0]=0.040;

  air_composition_N2[1]=air_composition_N2[2]=1.0;
  air_composition_O2[1]=air_composition_O2[2]=1.0;
  air_composition_Ar[1]=air_composition_Ar[2]=1.0;
  air_composition_CO2[1]=air_composition_CO2[2]=1.0;
  air_composition_N2[3]=1;
  air_composition_O2[3]=1;
  air_composition_Ar[3]=1;
  air_composition_CO2[3]=1;

  air_temperature[0]=15.0;
  air_pressure[0]=860.0;
  air_relative_humidity[0]=0.0;
  air_temperature[1]=air_temperature[2]=1.0;
  air_pressure[1]=air_pressure[2]=1.0;
  air_relative_humidity[1]=air_relative_humidity[2]=1.0;
  air_temperature[3]=1;
  air_pressure[3]=1;
  air_relative_humidity[3]=1;

  geomagnetic_onoff=0;
  for( int i=1; i<3; i++ ){ 
    geomagnetic_theta[i]=0.0;
    geomagnetic_declination[i]=0.0;
    geomagnetic_horizontal[i]=0.0;
    geomagnetic_vertical[i]=0.0;
  }

  cu_collimeter_d = 10.0;

  slit_width = 50.0;
 
  faradaycup_flag  = 0; 
  faradaycup4_flag = 0;

  r_faradaycup4_front1 = 90./2.0;
  length_faradaycup4_front1 = 0.1;
  positionz_faradaycup4_front1 = 15.7;
   
  r_faradaycup4_front2= 75./2.0;
  length_faradaycup4_front2 = 0.1;
  positionz_faradaycup4_front2= 20.8;

  r_faradaycup4_copper = 70.0/2.0;
  length_faradaycup4_copper = 60.0;
  positionz_faradaycup4_copper = 24.0;

  position_alongbeamline_faradaycup4=0.0;
  position_transverse_faradaycup4=0.0;  

  faradaycup5_height = 100.;

  screen_monitor3_flag = 0;
  positionz_screen_monitor3 = 159.5;
  
  beam_attenuator_flag = 0;
  beam_attenuator_material = 0;
  length_beam_attenuator = 11.0;

  pb_collimator_flag = 0;
  dinner_pb_collimator = 10.0;
  douter_pb_collimator = 55.0;
  length_pb_collimator = 25.0;

  virtual_chamber_flag = 0;
  r_virtual_chamber_injectionhole = 20.0; // unit=mm
  r_virtual_chamber_innerhole = 100.0;    // unit=mm
  length_virtual_chamber = 500.0;         // unit=mm
  gap_virtual_chamber = 94.0;             // unit=mm
  positionz_virtual_chamber = 5000.0;     // unit=mm

  vc_temperature = 15.0;
  vc_pressure = 1013.0;

  qm_magnetic_field_flag = 0;
  qm1_magnetic_field = 0.0;
  qm2_magnetic_field = 0.0;

  bm_magnetic_field_flag = 1;
  bm_magnetic_field      = 0.0;

  cu_collimeter_shiftz = 0.0;

  vertical_beamline_shift[0] = 0.0;
  vertical_beamline_shift[1] = 0.0; 

  bmyork_coil_shift[0] = 0.0;
  bmyork_coil_shift[1] = 0.0;
  bmyork_coil_shift[2] = 0.0;

  slit_collimeter_shift[0]=0.0;
  slit_collimeter_shift[1]=0.0;

  faradaycup_shift[0] = 0.0;
  faradaycup_shift[1] = 0.0;
  faradaycup_shift[2] = 0.0;

  coverbox_hole_shift[0] = 0.0;
  coverbox_hole_shift[1] = 0.0;
 
  photo_nuclear_process_flag = 0;
  cerenkov_process_flag = 0;

  cutlenght_particle = 1.0; //mm

  delta_pathlength = 0.10; // m

  position_of_mirror_6_BRM_coordinate[0] = -2.166254023; // unit=m
  position_of_mirror_6_BRM_coordinate[1] = 12.654034; // unit=m
  position_of_mirror_6_BRM_coordinate[2] = 0.0;   // 5.4494 m
 
  angle_ELS_mirror_6 = 9.0;   // deg
    
  distance_ELS_from_mirror6 = 100.0;  // m

  angle_north_yaxis_brm = 69.39492;

  mesh_dx = 0.1; // m
  mesh_dy = 0.1; // m
  mesh_dz = 0.1; // m

  mesh_xregion = 200.0; // m
  mesh_yregion = 200.0; // m
  mesh_zregion = 200.0; // m
  //---------------------------------------------

}
//------------------------------------------------
Plist::~Plist(){}
//------------------------------------------------
void Plist::readfile(char *setupfile)
{
  ifstream fin;
  fin.open(setupfile);

  cout <<  setupfile << endl;

  char c[256];
  while( fin.getline(c, sizeof c) ){

    if( !skipspaceline(c) ) continue;

    string s=c;
    if(s.empty()){
      continue;    
    }else{
      read_line(c);
    }
  }

  setup_parameters();
  dump_parameter();

  fin.close();

  return;
}
//------------------------------------------------
int Plist::skipspaceline( char *c )
{
  string s=c;
  int rv(0);

  string::size_type index;
  index = s.find_first_not_of(" ");
  if( index != string::npos ) rv=1;

  return rv;
}
//------------------------------------------------
int Plist::skipcomment( char *c )
{
  string s=c;
  int rv(0);

  string::size_type index1,index2,index3;

  index1 = s.find("*");
  index2 = s.find("#");
  index3 = s.find("//");
  if( index1 != string::npos ||
      index2 != string::npos ||
      index3 != string::npos ) rv=1;

  return rv;
}
//-------------------------------------------------
int Plist::read_parameterID( char* c )
{
  int param_ID(-1);

  if( strcasecmp( c, "beam_energy_mode" )         == 0 ) param_ID = 1;

  if( strcasecmp( c, "beam_particle_id" )         == 0 ) param_ID = 2;

  if( strcasecmp( c, "outoputfileflag_detectorplane")     == 0 ) param_ID = 5;
  if( strcasecmp( c, "outoputfileflag_secondparticle")    == 0 ) param_ID = 6;
  if( strcasecmp( c, "outoputfileflag_dedx")              == 0 ) param_ID = 7;
  if( strcasecmp( c, "outoputfileflag_energydeposit")     == 0 ) param_ID = 8;
  
  if( strcasecmp( c, "primary_energy")            == 0 ) param_ID = 11;

  if( strcasecmp( c, "crystalball_mean")          == 0 ) param_ID = 21;
  if( strcasecmp( c, "crystalball_sigma")         == 0 ) param_ID = 22;
  if( strcasecmp( c, "crystalball_alpha")         == 0 ) param_ID = 23;
  if( strcasecmp( c, "crystalball_n")             == 0 ) param_ID = 24;
  if( strcasecmp( c, "crystalball_xmin")          == 0 ) param_ID = 25;
  if( strcasecmp( c, "crystalball_xmax")          == 0 ) param_ID = 26;

  if( strcasecmp( c, "beam_emmitance_mode"  )     == 0 ) param_ID = 100;

  if( strcasecmp( c, "electronbeam_waveform_flag"      ) == 0 ) param_ID = 101;
  if( strcasecmp( c, "electronbeam_waveform_file"      ) == 0 ) param_ID = 102;
  if( strcasecmp( c, "electronbeam_waveform_dt"        ) == 0 ) param_ID = 103;
  if( strcasecmp( c, "electronbeam_waveform_new_ntbin" ) == 0 ) param_ID = 104;

  if( strcasecmp( c, "injection_position_x" )     == 0 ) param_ID = 111;
  if( strcasecmp( c, "injection_position_y" )     == 0 ) param_ID = 112;
  if( strcasecmp( c, "injection_position_z" )     == 0 ) param_ID = 113;

  if( strcasecmp( c, "injection_position_sigma" ) == 0 ) param_ID = 114;

  if( strcasecmp( c, "injection_direction_x" )    == 0 ) param_ID = 121;
  if( strcasecmp( c, "injection_direction_y" )    == 0 ) param_ID = 122;
  if( strcasecmp( c, "injection_direction_z" )    == 0 ) param_ID = 123;

  if( strcasecmp( c, "injection_transverse_phi0_x" ) == 0 ) param_ID = 125;
  if( strcasecmp( c, "injection_transverse_phi0_y" ) == 0 ) param_ID = 126;

  if( strcasecmp( c, "injection_position_shift_x") == 0 ) param_ID = 131;
  if( strcasecmp( c, "injection_position_sigma_x") == 0 ) param_ID = 132;
  if( strcasecmp( c, "injection_position_shift_y") == 0 ) param_ID = 133;
  if( strcasecmp( c, "injection_position_sigma_y") == 0 ) param_ID = 134; 

  if( strcasecmp( c, "beam_emmitance_h_gamma" )      == 0 ) param_ID = 151;
  if( strcasecmp( c, "beam_emmitance_h_beta" )       == 0 ) param_ID = 152;
  if( strcasecmp( c, "beam_emmitance_h_alpha" )      == 0 ) param_ID = 153;
  if( strcasecmp( c, "beam_emmitance_h_emmitance" )   == 0 ) param_ID = 154;

  if( strcasecmp( c, "beam_emmitance_v_gamma" )      == 0 ) param_ID = 161;
  if( strcasecmp( c, "beam_emmitance_v_beta" )       == 0 ) param_ID = 162;
  if( strcasecmp( c, "beam_emmitance_v_alpha" )      == 0 ) param_ID = 163;
  if( strcasecmp( c, "beam_emmitance_v_emmitance" )   == 0 ) param_ID = 164;
 
  // Atmosphere Condition
  if( strcasecmp( c, "air_composition_N2" )    == 0 ) param_ID = 501;
  if( strcasecmp( c, "air_composition_O2" )    == 0 ) param_ID = 502;
  if( strcasecmp( c, "air_composition_Ar" )    == 0 ) param_ID = 503;
  if( strcasecmp( c, "air_composition_CO2" )   == 0 ) param_ID = 504;

  if( strcasecmp( c, "air_temperature" )         == 0 ) param_ID = 511;
  if( strcasecmp( c, "air_pressure" )            == 0 ) param_ID = 512;
  if( strcasecmp( c, "air_relative_humidity" )   == 0 ) param_ID = 513;

  // Geomagnetic Field @ BRM ( ELS )                                  
  if(  strcasecmp( c, "geomagnetic_onoff" )      == 0 ) param_ID = 700;
  if(  strcasecmp( c, "geomagnetic_theta" )       ==0 ) param_ID = 701;    
  if(  strcasecmp( c, "geomagnetic_declination" ) ==0 ) param_ID = 702;
  if(  strcasecmp( c, "geomagnetic_horizontal"  ) ==0 ) param_ID = 703;
  if(  strcasecmp( c, "geomagnetic_vertical"    ) ==0 ) param_ID = 704;

  // World Geometory
  if(  strcasecmp( c, "size_of_world" ) == 0 )  param_ID = 1000;

  // ELS Site Ground Height
  if( strcasecmp( c, "elssite_ground_height" ) == 0 ) param_ID = 1500;

  // Cu Collimter Diameter
  if( strcasecmp( c, "cu_collimeter_d" ) == 0 ) param_ID = 2000;

  // Slit width
  if( strcasecmp( c, "slit_width" ) == 0 ) param_ID = 2500;

  // -- Beam Dump Flag
  if( strcasecmp( c, "faradaycup_flag" ) == 0 ) param_ID = 3000;

  if( strcasecmp( c, "faradaycup2_copper_flag" ) == 0 ){}           
  if( strcasecmp( c, "r_faradaycup2_copper" ) == 0 ){}              
  if( strcasecmp( c, "length_faradaycup2_copper" ) == 0 ){}         
  if( strcasecmp( c, "positionz_faradaycup2" ) == 0 ){}             
  if( strcasecmp( c, "positionz_faradaycup2_copper" ) == 0 ){}      
  if( strcasecmp( c, "position_alongbeamline_faradaycup2" ) == 0 ){}
  if( strcasecmp( c, "position_transverse_faradaycup2" ) == 0 ){}   
  if( strcasecmp( c, "length_faradaycup2_carbon" ) == 0 ){}         

  if( strcasecmp( c, "faradaycup4_flag" ) == 0 )             param_ID = 3010;

  if( strcasecmp( c, "r_faradaycup4_front1" ) == 0 )         param_ID = 3011;
  if( strcasecmp( c, "length_faradaycup4_front1" ) == 0 )    param_ID = 3012;
  if( strcasecmp( c, "positionz_faradaycup4_front1" ) == 0 ) param_ID = 3013;
  if( strcasecmp( c, "r_faradaycup4_front2" ) == 0 )         param_ID = 3014;
  if( strcasecmp( c, "length_faradaycup4_front2" ) == 0 )    param_ID = 3015;
  if( strcasecmp( c, "positionz_faradaycup4_front2" ) == 0 ) param_ID = 3016;
  if( strcasecmp( c, "r_faradaycup4_copper" ) == 0 )         param_ID = 3017;
  if( strcasecmp( c, "length_faradaycup4_copper" ) == 0 )    param_ID = 3018;
  if( strcasecmp( c, "positionz_faradaycup4_copper" ) == 0 ) param_ID = 3019;

  if( strcasecmp( c, "position_alongbeamline_faradaycup4" ) == 0 ) param_ID = 3020;
  if( strcasecmp( c, "position_transverse_faradaycup4" ) == 0 )    param_ID = 3021;  

  if( strcasecmp( c, "faradaycup5_height" ) == 0 )    param_ID = 3022;

  if( strcasecmp( c, "screen_monitor3_flag"     ) == 0 ) param_ID = 3025;  
  if( strcasecmp( c, "positionz_screen_monitor3") == 0 ) param_ID = 3026;

  if( strcasecmp( c, "beam_attenuator_flag"     ) == 0 ) param_ID = 3030;
  if( strcasecmp( c, "beam_attenuator_material" ) == 0 ) param_ID = 3031;       
  if( strcasecmp( c, "length_beam_attenuator"   ) == 0 ) param_ID = 3032;

  // Pb Collimator added in 2012.09.24
  if( strcasecmp( c, "pb_collimator_flag"       ) == 0 ) param_ID = 3040;
  if( strcasecmp( c, "dinner_pb_collimator"     ) == 0 ) param_ID = 3041;
  if( strcasecmp( c, "douter_pb_collimator"     ) == 0 ) param_ID = 3042;
  if( strcasecmp( c, "length_pb_collimator"     ) == 0 ) param_ID = 3043;  

  // Virtual Test Chamber added in 2011.12.13
  if( strcasecmp( c, "virtual_chamber_flag"            ) == 0 ) param_ID = 3090;
  if( strcasecmp( c, "r_virtual_chamber_injectionhole" ) == 0 ) param_ID = 3091;
  if( strcasecmp( c, "r_virtual_chamber_innerhole"     ) == 0 ) param_ID = 3092;
  if( strcasecmp( c, "length_virtual_chamber"          ) == 0 ) param_ID = 3093;
  if( strcasecmp( c, "gap_virtual_chamber"             ) == 0 ) param_ID = 3094;
  if( strcasecmp( c, "positionz_virtual_chamber"       ) == 0 ) param_ID = 3095;  
  if( strcasecmp( c, "vc_temperature"                  ) == 0 ) param_ID = 3096;
  if( strcasecmp( c, "vc_pressure"                     ) == 0 ) param_ID = 3097; 

  // -- Magnetic Field Flag
  if( strcasecmp( c, "qm_magnetic_field_flag" ) == 0 ) param_ID = 3100; 
  if( strcasecmp( c, "qm1_magnetic_field" ) == 0 ) param_ID = 3101; 
  if( strcasecmp( c, "qm2_magnetic_field" ) == 0 ) param_ID = 3102; 

  if( strcasecmp( c, "bm_magnetic_field_flag" ) == 0 ) param_ID = 3200; 
  if( strcasecmp( c, "bm_magnetic_field"      ) == 0 ) param_ID = 3201;

  // -- Alignment Parameters  
  if( strcasecmp( c, "cu_collimeter_shiftz" ) == 0 ) param_ID = 3501;

  if( strcasecmp( c, "vertical_beamline_shift_x" ) == 0 ) param_ID = 3511;
  if( strcasecmp( c, "vertical_beamline_shift_y" ) == 0 ) param_ID = 3512;

  if( strcasecmp( c, "bmyork_coil_shift_x" ) == 0 ) param_ID = 3521;
  if( strcasecmp( c, "bmyork_coil_shift_y" ) == 0 ) param_ID = 3522;
  if( strcasecmp( c, "bmyork_coil_shift_z" ) == 0 ) param_ID = 3523;

  if( strcasecmp( c, "slit_collimeter_shift_0" ) == 0 ) param_ID = 3531;
  if( strcasecmp( c, "slit_collimeter_shift_1" ) == 0 ) param_ID = 3532;
  
  if( strcasecmp( c, "faradaycup_shift_x" ) == 0 ) param_ID = 3541;
  if( strcasecmp( c, "faradaycup_shift_y" ) == 0 ) param_ID = 3542;
  if( strcasecmp( c, "faradaycup_shift_z" ) == 0 ) param_ID = 3543;

  if( strcasecmp( c, "coverbox_hole_shift_x") == 0 ) param_ID = 3551;
  if( strcasecmp( c, "coverbox_hole_shift_y") == 0 ) param_ID = 3552;
  
  // -- Physics Proccess Flag
  if( strcasecmp( c, "photo_nuclear_process_flag" )  == 0 ) param_ID = 4000; 
  if( strcasecmp( c, "cerenkov_process_flag" )  == 0 ) param_ID = 4001; 

  if( strcasecmp( c, "cutlenght_particle" ) == 0 ) param_ID = 4010;

  // -- Energy Deposit Mesh
  if( strcasecmp( c, "delta_pathlength" )  == 0 ) param_ID = 5001;

  if( strcasecmp( c, "position_of_mirror_6_BRM_coordinate_x" )  == 0 ) param_ID = 5010;
  if( strcasecmp( c, "position_of_mirror_6_BRM_coordinate_y" )  == 0 ) param_ID = 5011;
  if( strcasecmp( c, "position_of_mirror_6_BRM_coordinate_z" )  == 0 ) param_ID = 5012;

  if( strcasecmp( c, "angle_ELS_mirror_6" )         == 0 ) param_ID = 5020;
  if( strcasecmp( c, "distance_ELS_from_mirror6" )  == 0 ) param_ID = 5021;

  if( strcasecmp( c, "angle_north_yaxis_brm" )  == 0 ) param_ID = 5022;

  if( strcasecmp( c, "mesh_dx" )  == 0 ) param_ID = 5030;
  if( strcasecmp( c, "mesh_dy" )  == 0 ) param_ID = 5031;
  if( strcasecmp( c, "mesh_dz" )  == 0 ) param_ID = 5032;

  if( strcasecmp( c, "mesh_xregion" )  == 0 ) param_ID = 5035;
  if( strcasecmp( c, "mesh_yregion" )  == 0 ) param_ID = 5036;
  if( strcasecmp( c, "mesh_zregion" )  == 0 ) param_ID = 5037;

  return param_ID;
}
//-------------------------------------------------
void Plist::read_parameter( int paramID, int paramN, char* c )
{
  
  // Primary Electron //
  if( paramID == 1   ) beam_energy_mode            = atoi(c);
 
  if( paramID == 2   ) beam_particle_id           = atoi(c); 

  if( paramID == 5   ) outoputfileflag_detectorplane  = atoi(c);
  if( paramID == 6   ) outoputfileflag_secondparticle = atoi(c);
  if( paramID == 7   ) outoputfileflag_dedx           = atoi(c);
  if( paramID == 8   ) outoputfileflag_energydeposit  = atoi(c);

  if( paramID == 11  ) primary_energy[paramN]      = atof(c);
  
  if( paramID == 21 ) crystalball_mean[paramN]     = atof(c);
  if( paramID == 22 ) crystalball_sigma[paramN]    = atof(c);
  if( paramID == 23 ) crystalball_alpha[paramN]    = atof(c);
  if( paramID == 24 ) crystalball_n[paramN]        = atof(c);
  if( paramID == 25 ) crystalball_xmin[paramN]     = atof(c);
  if( paramID == 26 ) crystalball_xmax[paramN]     = atof(c);

  if( paramID == 100 ) beamemmitance_mode          = atoi(c);

  if( paramID == 101 ) electronbeam_waveform_flag  = atoi(c);
  if( paramID == 102 ){ sprintf( electronbeam_waveform_file, "%s", c ); }
  if( paramID == 103 ) electronbeam_waveform_dt        = atof(c);
  if( paramID == 104 ) electronbeam_waveform_new_ntbin = atoi(c); 

  if( paramID == 111 ) injection_position_x[paramN] = atof(c);
  if( paramID == 112 ) injection_position_y[paramN] = atof(c);
  if( paramID == 113 ) injection_position_z[paramN] = atof(c);

  if( paramID == 114 ) injection_position_sigma[paramN] = atof(c);

  if( paramID == 121 ) injection_direction_x[paramN] = atof(c);
  if( paramID == 122 ) injection_direction_y[paramN] = atof(c);
  if( paramID == 123 ) injection_direction_z[paramN] = atof(c);

  if( paramID == 125 ) injection_transverse_phi0_x[paramN] = atof(c);
  if( paramID == 126 ) injection_transverse_phi0_y[paramN] = atof(c);

  if( paramID == 131 ) injection_position_shift_x[paramN] = atof(c);
  if( paramID == 132 ) injection_position_sigma_x[paramN] = atof(c);
  if( paramID == 133 ) injection_position_shift_y[paramN] = atof(c);
  if( paramID == 134 ) injection_position_sigma_y[paramN] = atof(c);

  if( paramID == 151 ) injection_beam_emmitance_g[0][paramN] = atof(c);
  if( paramID == 152 ) injection_beam_emmitance_b[0][paramN] = atof(c);
  if( paramID == 153 ) injection_beam_emmitance_a[0][paramN] = atof(c);
  if( paramID == 154 ) injection_beam_emmitance_e[0][paramN] = atof(c);

  if( paramID == 161 ) injection_beam_emmitance_g[1][paramN] = atof(c);
  if( paramID == 162 ) injection_beam_emmitance_b[1][paramN] = atof(c);
  if( paramID == 163 ) injection_beam_emmitance_a[1][paramN] = atof(c);
  if( paramID == 164 ) injection_beam_emmitance_e[1][paramN] = atof(c);

  // Atmosphere Condition
  if( paramID == 501 ) air_composition_N2[paramN] = atof(c);
  if( paramID == 502 ) air_composition_O2[paramN] = atof(c);
  if( paramID == 503 ) air_composition_Ar[paramN] = atof(c);
  if( paramID == 504 ) air_composition_CO2[paramN] = atof(c);

  if( paramID == 511 ) air_temperature[paramN] = atof(c);
  if( paramID == 512 ) air_pressure[paramN] = atof(c);
  if( paramID == 513 ) air_relative_humidity[paramN] = atof(c);

  // Geomagnetic Field @ BRM ( ELS )
  if( paramID == 700 ) geomagnetic_onoff               = atoi(c);
  if( paramID == 701 ) geomagnetic_theta[paramN]       = atof(c);
  if( paramID == 702 ) geomagnetic_declination[paramN] = atof(c);
  if( paramID == 703 ) geomagnetic_horizontal[paramN]  = atof(c);
  if( paramID == 704 ) geomagnetic_vertical[paramN]    = atof(c);

  // World Geometory
  if( paramID == 1000 ) size_of_world[paramN] = atof(c);

  // ELS Site Ground Height    
  if( paramID == 1500 ) elssite_ground_height[paramN] = atof(c);

  // Cu Collimter Diameter
  if( paramID == 2000 ) cu_collimeter_d = atof(c);

  // Slit width
  if( paramID == 2500 ) slit_width = atof(c);

  // -- Beam Dump Flag ( Faraday Cup )
  if( paramID == 3000 ) faradaycup_flag = atoi(c);

  if( paramID == 3010 ) faradaycup4_flag             = atoi(c);
  if( paramID == 3011 ) r_faradaycup4_front1         = atof(c);
  if( paramID == 3012 ) length_faradaycup4_front1    = atof(c);
  if( paramID == 3013 ) positionz_faradaycup4_front1 = atof(c);
  if( paramID == 3014 ) r_faradaycup4_front2         = atof(c);
  if( paramID == 3015 ) length_faradaycup4_front2    = atof(c);
  if( paramID == 3016 ) positionz_faradaycup4_front2 = atof(c);
  if( paramID == 3017 ) r_faradaycup4_copper         = atof(c);
  if( paramID == 3018 ) length_faradaycup4_copper    = atof(c);
  if( paramID == 3019 ) positionz_faradaycup4_copper = atof(c);

  if( paramID == 3020 ) position_alongbeamline_faradaycup4 = atof(c);
  if( paramID == 3021 ) position_transverse_faradaycup4    = atof(c);  

  if( paramID == 3022 ) faradaycup5_height = atof(c);

  if( paramID == 3025 ) screen_monitor3_flag         = atoi(c);
  if( paramID == 3026 ) positionz_screen_monitor3    = atof(c);

  if( paramID == 3030 ) beam_attenuator_flag         = atoi(c);
  if( paramID == 3031 ) beam_attenuator_material     = atoi(c);
  if( paramID == 3032 ) length_beam_attenuator       = atof(c); 

  // Pb Collimator added in 2012.09.24
  if( paramID == 3040 ) pb_collimator_flag     = atoi(c);
  if( paramID == 3041 ) dinner_pb_collimator   = atof(c);
  if( paramID == 3042 ) douter_pb_collimator   = atof(c);
  if( paramID == 3043 ) length_pb_collimator   = atof(c);

  // Virtual Test Chamber added in 2011.12.13
  if( paramID == 3090 ) virtual_chamber_flag            = atoi(c); 
  if( paramID == 3091 ) r_virtual_chamber_injectionhole = atof(c);
  if( paramID == 3092 ) r_virtual_chamber_innerhole     = atof(c);
  if( paramID == 3093 ) length_virtual_chamber          = atof(c);
  if( paramID == 3094 ) gap_virtual_chamber             = atof(c);
  if( paramID == 3095 ) positionz_virtual_chamber       = atof(c);    
  if( paramID == 3096 ) vc_temperature                  = atof(c);
  if( paramID == 3097 ) vc_pressure                     = atof(c);
  
  // -- Magnetic Field Flag
  if( paramID == 3100 ) qm_magnetic_field_flag = atoi(c);
  if( paramID == 3101 ) qm1_magnetic_field = atof(c);
  if( paramID == 3102 ) qm2_magnetic_field = atof(c);

  if( paramID == 3200 ) bm_magnetic_field_flag = atoi(c);
  if( paramID == 3201 ) bm_magnetic_field      = atof(c);

  // -- Alignment Parameters            
  if( paramID == 3501 ) cu_collimeter_shiftz   = atof(c);
   
  if( paramID == 3511 ) vertical_beamline_shift[0] = atof(c);
  if( paramID == 3512 ) vertical_beamline_shift[1] = atof(c);

  if( paramID == 3521 ) bmyork_coil_shift[0] = atof(c);
  if( paramID == 3522 ) bmyork_coil_shift[1] = atof(c);
  if( paramID == 3523 ) bmyork_coil_shift[2] = atof(c);
 
  if( paramID == 3531 ) slit_collimeter_shift[0] = atof(c);
  if( paramID == 3532 ) slit_collimeter_shift[1] = atof(c);

  if( paramID == 3541 ) faradaycup_shift[0] = atof(c);
  if( paramID == 3542 ) faradaycup_shift[1] = atof(c);
  if( paramID == 3543 ) faradaycup_shift[2] = atof(c);
  
  if( paramID == 3551 ) coverbox_hole_shift[0] = atof(c);
  if( paramID == 3552 ) coverbox_hole_shift[1] = atof(c);

  // -- Physics Proccess Flag
  if( paramID == 4000 ) photo_nuclear_process_flag = atoi(c); 
  if( paramID == 4001 ) cerenkov_process_flag = atoi(c); 

  if( paramID == 4010 ) cutlenght_particle = atof(c); 

  // -- Energy Deposit Mesh        
  if( paramID == 5001 ) delta_pathlength = atof(c);

  if( paramID == 5010 ) position_of_mirror_6_BRM_coordinate[0] = atof(c);
  if( paramID == 5011 ) position_of_mirror_6_BRM_coordinate[1] = atof(c);
  if( paramID == 5012 ) position_of_mirror_6_BRM_coordinate[2] = atof(c);

  if( paramID == 5020 ) angle_ELS_mirror_6         = atof(c);
  if( paramID == 5021 ) distance_ELS_from_mirror6  = atof(c);
  if( paramID == 5022 ) angle_north_yaxis_brm  = atof(c);

   
  if( paramID == 5030 ) mesh_dx = atof(c);
  if( paramID == 5031 ) mesh_dy = atof(c);
  if( paramID == 5032 ) mesh_dz = atof(c);

  if( paramID == 5035 ) mesh_xregion = atof(c);
  if( paramID == 5036 ) mesh_yregion = atof(c);
  if( paramID == 5037 ) mesh_zregion = atof(c);

  return;
}
//----------------------------------------------------------------------
void Plist::setup_parameters(void)
{
  return;
}
//----------------------------------------------------------------------
void Plist::dump_parameter(void)
{
  return;
}
//----------------------------------------------------------------------
void Plist::read_line( char *c )
{
  int i(0);
  char *token;
  token = strtok( c , " " );

  if( skipcomment(token)==1 ) return;
  i=read_parameterID(token);

  int j(0);
  while( token != NULL ){
    token = strtok( NULL , " " );
    if( token!=NULL ) {
      if( skipcomment(token)==1 ) break;
      read_parameter(i,j,token);
      j++;
    }
  }

  return;
}
//----------------------------------------------------------------------
