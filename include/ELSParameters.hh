//===================================================================
//    ELS Parameter Class : ELSParameters.hh
//      Author    : T.Shibata
//
//      Creation  : 2009.11.30
//    Last Update : 2010.11.09
//    Last Update : 2011.01.24
//    Last Update : 2011.02.03
//    Last Update : 2011.04.12
//    Last Update : 2011.05.15
//    Last Update : 2011.12.11
//    Last Update : 2011.12.12
//    Last Update : 2011.12.13
//    Last Update : 2011.12.14
//    Last Update : 2011.12.25
//    Last Update : 2012.01.07
//    Last Update : 2012.01.25
//    Last Update : 2012.05.01
//    Last Update : 2012.09.24
//    Last Update : 2012.10.02
//===================================================================
#include <iostream>
#include <fstream>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
 
#include <CLHEP/Vector/ThreeVector.h>
#include "G4ThreeVector.hh"
 
#include "ReadFile.hh"

#ifndef ELSParameters_h
#define ELSParameters_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class Plist;
class ELSParameters{
public:
  // Construction
  ELSParameters();
  // Deconstruction
  ~ELSParameters();

  void els_geometory_initialize(Plist *plist);

  void els_alignment_setting(Plist *plist);

  void els_parameter_initialize(Plist *plist);

  // Return Concrete Pad ( G4Box )
  double sizeConcretepad(int i){ return size_of_concretepad[i]; };
  double positionConcretepad(int i){ return position_of_concretepad[i]; };
  double rotationConcretepad(int i){ return rotation_of_concretepad[i]; };

  // Return ELS Container ( G4Box )
  double sizeOuterELSContainer(int i){ return outer_size_of_ELScontainer[i]; };
  double sizeInnerELSContainer(int i){ return inner_size_of_ELScontainer[i]; };
  double positionOuterELSContainer(int i){ return position_of_ELS_outer[i]; };
  double positionInnerELSContainer(int i){ return position_of_ELS_inner[i]; };
  double rotationOuterELSContainer(int i){ return rotation_of_ELS_outer[i]; };
  double rotationInnerELSContainer(int i){ return rotation_of_ELS_inner[i]; };
  double positionELSContainer(int i){ return position_of_ELS[i]; };
  double rotationELSContainer(int i){ return rotation_of_ELS[i]; };

  double sizeELSinjectionhole(int i){ return size_of_ELS_injection_hole[i]; };
  double positionELSinjectionhole(int i){ return position_of_ELS_injection_hole[i]; };

  double sizeOuterCoverbox(int i){ return size_of_coverbox_outer[i]; };
  double sizeInnerCoverbox(int i){ return size_of_coverbox_inner[i]; };
  double positionOuterCoverbox(int i){ return position_of_coverbox_outer[i]; };
  double positionInnerCoverbox(int i){ return position_of_coverbox_inner[i]; };

  double sizeCoverboxBoard(int i){ return size_of_coverbox_board[i]; };
  double rminCoverboxBoardHole(void){ return rmin_coverbox_boardhole; };
  double rmaxCoverboxBoardHole(void){ return rmax_coverbox_boardhole; };
  double lCoverboxBoardHole(void){ return l_coverbox_boardhole; };
  double positionCoverboxBoard(int i){ return position_of_coverbox_board[i]; };
  double positionCoverboxBoardHole(int i){ return position_of_coverbox_boardhole[i]; };

  // Return Base A & B & C ( G4Box )
  double sizeBase(int i, int j){ return size_of_base[i][j]; };
  double positionBase(int i, int j){ return position_of_base[i][j]; };
  double rotationBase(int i, int j){ return rotation_of_base[i][j]; };

  // Return BeamLine Component ( G4Tube )
  double numBLcomp(void){ return NumBLComp; };
  double rminBL(int i){ return RMinBL[i]; };
  double rmaxBL(int i){ return RMaxBL[i]; };
  double lBL(int i){ return L_BL[i]; };
  double rminBLff(int i){ return RMinBL_FF[i]; };
  double rmaxBLff(int i){ return RMaxBL_FF[i]; };
  double lBLff(int i){ return L_BL_FF[i]; };
  double rminBLbf(int i){ return RMinBL_BF[i]; };
  double rmaxBLbf(int i){ return RMaxBL_BF[i]; };
  double lBLbf(int i){ return L_BL_BF[i]; };
  double positionBLff(int i, int j){ return position_BL_FF[i][j]; };
  double positionBL(int i, int j){ return position_BL[i][j]; };
  double positionBLbf(int i, int j){ return position_BL_BF[i][j];}; 
  double rotationBL(int i, int j){ return rotation_BL[i][j]; };

  // Blank Frange
  double rminEGUNBlank(void){ return RMin_EGUNBlank; };
  double rmaxEGUNBlank(void){ return RMax_EGUNBlank; };
  double lEGUNBlank(void){ return L_EGUNBlank; };
  double positionEGUNBlank(int i){ return position_EGUNBlank[i]; };
  double rotationEGUNBlank(int i){ return rotation_EGUNBlank[i]; };

  // EGUN Anode
  double rminAnode(void){ return RMin_Anode; };
  double rmaxAnode(void){ return RMax_Anode; };
  double lAnode(void){ return L_Anode; };
  double positionAnode(int i){ return position_Anode[i];};
  double rotationAnode(int i){ return rotation_Anode[i]; };

  // Return BeamLine Component ( Screen Monitor  )
  // parts-0 ( Frange )
  double rminSMff( int i ){ if( i==0 ) return RMin_SM1_FF;
                            if( i==1 ) return RMin_SM2_FF;
			    return 0;};
  double rmaxSMff( int i ){ if( i==0 ) return RMax_SM1_FF;
                            if( i==1 ) return RMax_SM2_FF;
			    return 0;};
  double lSMff( int i ){ if( i==0 ) return L_SM1_FF; 
                         if( i==1 ) return L_SM2_FF; 
			 return 0;};
  double positionSMff( int i, int j ){ if( i==0 ) return position_SM1_FF[j]; 
                                       if( i==1 ) return position_SM2_FF[j]; 
				       return 0;}; 
  double rminSMbf( int i ){ if( i==0 ) return RMin_SM1_BF;
                            if( i==1 ) return RMin_SM2_BF; 
			    return 0;};
  double rmaxSMbf( int i ){ if( i==0 ) return RMax_SM1_BF;
                            if( i==1 ) return RMax_SM2_BF; 
			    return 0;};
  double lSMbf( int i ){ if( i==0 ) return L_SM1_BF; 
                         if( i==1 ) return L_SM2_BF; 
			 return 0;};
  double positionSMbf( int i, int j ){ if( i==0 ) return position_SM1_BF[j]; 
                                       if( i==1 ) return position_SM2_BF[j]; 
				       return 0;};
  double rotationSMf( int i, int j ){ if( i==0 ) return rotation_SM1_F[j]; 
                                      if( i==1 ) return rotation_SM2_F[j]; 
				      return 0;};
  // parts-1
  double rminSMparts1( int i ){ if( i==0 ) return Rmin_SM1_parts1; 
                               if( i==1 ) return Rmin_SM2_parts1; 
			       return 0; };
  double rmaxSMparts1( int i ){ if( i==0 ) return Rmax_SM1_parts1;
                               if( i==1 ) return Rmax_SM2_parts1; 
			       return 0; };
  double lSMparts1( int i ){ if( i==0 ) return L_SM1_parts1;
                            if( i==1 ) return L_SM2_parts1; 
			    return 0; };
  double positionSMparts1( int i, int j ){ if( i==0 ) return position_SM1_parts1[j];
                                           if( i==1 ) return position_SM2_parts1[j]; 
					   return 0; }; 
  // parts-2
  double rminSMparts2( int i ){ if( i==0 ) return Rmin_SM1_parts2;
                                if( i==1 ) return Rmin_SM2_parts2; 
				return 0; };
  double rmaxSMparts2( int i ){ if( i==0 ) return Rmax_SM1_parts2;
                                if( i==1 ) return Rmax_SM2_parts2; 
				return 0; };
  double lSMparts2( int i ){ if( i==0 ) return L_SM1_parts2;
                             if( i==1 ) return L_SM2_parts2; 
			     return 0;};
  double positionSMparts2( int i, int j ){ if( i==0 ) return position_SM1_parts2[j]; 
                                           if( i==1 ) return position_SM2_parts2[j]; 
					   return 0;};
  // parts-3
  double rminSMparts3( int i ){ if( i==0 ) return Rmin_SM1_parts3;
                                if( i==1 ) return Rmin_SM2_parts3; 
				return 0;}
  double rmaxSMparts3( int i ){ if( i==0 ) return Rmax_SM1_parts3;
                                if( i==1 ) return Rmax_SM2_parts3; 
				return 0;}
  double lSMparts3( int i ){ if( i==0 ) return L_SM1_parts3;
                             if( i==1 ) return L_SM2_parts3; 
			     return 0;};
  double positionSMparts3( int i, int j ){ if( i==0 ) return position_SM1_parts3[j];
                                           if( i==1 ) return position_SM2_parts3[j]; 
					   return 0;};
  double rotationSMcomp( int i, int j ){ if( i==0 ) return rotation_SM1_comp[j]; 
                                         if( i==1 ) return rotation_SM2_comp[j]; 
					 return 0;};
  double rotationSMduct( int i, int j ){ if( i==0 ) return rotation_SM1_duct[j]; 
                                         if( i==1 ) return rotation_SM2_duct[j]; 
					 return 0;};
  double rotationSMbody( int i, int j ){ if( i==0 ) return rotation_SM1_body[j]; 
                                         if( i==1 ) return rotation_SM2_body[j]; 
					 return 0;};

  // Return BeamLine Component ( 90-deg Bending Magnet Duct  )

  //  Cu-Collimeter
  double sizeCuCollimeter( int i ){ return size_of_CuCollimeter[i]; };
  double positionCuCollimeter( int i ){ return position_CuCollimeter[i]; };
  double rotationCuCollimeter( int i ){ return rotation_CuCollimeter[i]; };

  double rminCuCollimeterhole(void){ return RMin_CuCollimeter_hole; };
  double rmaxCuCollimeterhole(void){ return RMax_CuCollimeter_hole; };
  double lCuCollimeterhole(void){ return L_CuCollimeter_hole; };

  double positionCuCollimeterhole( int i ){ return position_CuCollimeter_hole[i]; };
  double rotationCuCollimeterhole( int i ){ return rotation_CuCollimeter_hole[i]; };

  // BM Duct Frange
  double rminBMf(void){ return Rmin_V85_BM; };
  double rmaxBMf(void){ return Rmax_V85_BM; };
  double lBMf(void){ return L_V85_BM; };  
  double sizeInnerBMf(int i){ return size_of_V85_BM_IN[i]; };
  double positionBMf(int i, int j ){  return position_BM_F[i][j]; };
  double rotationBMf(int i, int j ){ return rotation_BM_F[i][j]; };
  
  // BM Duct parts-1-3
  double sizeOuterBML( int i, int j ){ if(i==0) return size_of_BM_parts_1_outer[j]; 
                                       if(i==1) return size_of_BM_parts_2_outer[j];
                                       if(i==2) return size_of_BM_parts_3_outer[j]; 
				       return 0;};
  double sizeInnerBML( int i, int j ){ if(i==0) return size_of_BM_parts_1_inner[j];
                                       if(i==1) return size_of_BM_parts_2_inner[j];
                                       if(i==2) return size_of_BM_parts_3_inner[j]; 
				       return 0;};
  double positionOuterBML( int i, int j ){ if(i==0) return position_of_BM_parts_1_outer[j]; 
                                           if(i==1) return position_of_BM_parts_2_outer[j];
                                           if(i==2) return position_of_BM_parts_3_outer[j]; 
					   return 0;}
  double positionInnerBML( int i, int j ){ if(i==0) return position_of_BM_parts_1_inner[j];
                                           if(i==1) return position_of_BM_parts_2_inner[j];
                                           if(i==2) return position_of_BM_parts_3_inner[j]; 
					   return 0;}

  // BM Duct parts-1-3 Duct vacuum
  double sizeBMLvac( int i, int j ){  if(i==0) return size_of_BM_parts_1_vac[j];
                                      if(i==1) return size_of_BM_parts_2_vac[j];
                                      if(i==2) return size_of_BM_parts_3_vac[j]; 
				      return 0;};
  double positionBMLvac( int i, int j ){ if(i==0) return position_of_BM_parts_1_vac[j];
                                         if(i==1) return position_of_BM_parts_2_vac[j];
                                         if(i==2) return position_of_BM_parts_3_vac[j]; 
					 return 0;}

  // BM Duct parts-4(1-8)
  double sizeBMparts4m(int i, int j ){ if( i==0 ) return size_of_BM_parts_4_1[j]; 
                                       if( i==1 ) return size_of_BM_parts_4_2[j]; 
				       return 0;};
  double positionBMparts4m(int i, int j ){ if( i==0 )return position_of_BM_parts_4_1[j];
                                           if( i==1 ) return position_of_BM_parts_4_2[j]; 
					   return 0;};
  double rminBMparts43(void){ return Rmin_BM_parts_4_3; };
  double rmaxBMparts43(void){ return Rmax_BM_parts_4_3; };
  double lBMparts43(void){ return L_BM_parts_4_3; };
  double positionBMparts43(int i){ return position_of_BM_parts_4_3[i]; };
  double rotationBMparts43(int i){ return rotation_of_BM_parts_4_3[i]; };

  double sizeBMparts44(int i){ return size_of_BM_parts_4_4[i]; };
  double positionBMparts44(int i){  return position_of_BM_parts_4_4[i]; };

  double rminBMparts45(void){ return Rmin_BM_parts_4_5; };
  double rmaxBMparts45(void){ return Rmax_BM_parts_4_5; };
  double lBMparts45(void){ return L_BM_parts_4_5; };
  double positionBMparts45(int i){ return position_of_BM_parts_4_5[i]; };
  double rotationBMparts45(int i){ return rotation_of_BM_parts_4_5[i]; };

  double rminBMparts46(void){ return Rmin_BM_parts_4_6; };
  double rmaxBMparts46(void){ return Rmax_BM_parts_4_6; };
  double lBMparts46(void){ return L_BM_parts_4_6; };
  double phistartBMparts46(void){ return Phi_start_BM_parts_4_6; };
  double phistopBMparts46(void){ return Phi_end_BM_parts_4_6; };
  double positionBMparts46(int i){ return position_of_BM_parts_4_6[i]; };
  double rotationBMparts46(int i){ return rotation_of_BM_parts_4_6[i]; };
   
  double sizeBMparts47(int i){ return size_of_BM_parts_4_7[i]; };
  double positionBMparts47(int i){  return position_of_BM_parts_4_7[i]; };

  double rminBMparts48(void){ return Rmin_BM_parts_4_8; };
  double rmaxBMparts48(void){ return Rmax_BM_parts_4_8; };
  double lBMparts48(void){ return L_BM_parts_4_8; };
  double phistartBMparts48(void){ return Phi_start_BM_parts_4_8; };
  double phistopBMparts48(void){ return Phi_end_BM_parts_4_8; };
  double positionBMparts48(int i){ return position_of_BM_parts_4_8[i]; };
  double rotationBMparts48(int i){ return rotation_of_BM_parts_4_8[i]; };

  double sizeBMparts49(int i){ return size_of_BM_parts_4_9[i]; };
  double positionBMparts49(int i){ return position_of_BM_parts_4_9[i]; };
     
  // BM Main Duct vacuum  
  double sizeBMparts51(int i){ return size_of_BM_parts_5_1[i]; };
  double positionBMparts51(int i){  return position_of_BM_parts_5_1[i]; };

  double rminBMparts52(void){ return Rmin_BM_parts_5_2; };
  double rmaxBMparts52(void){ return Rmax_BM_parts_5_2; };
  double lBMparts52(void){ return L_BM_parts_5_2; };
  double positionBMparts52(int i){ return position_of_BM_parts_5_2[i]; };
  double rotationBMparts52(int i){ return rotation_of_BM_parts_5_2[i]; };

  double sizeBMparts53(int i){ return size_of_BM_parts_5_3[i]; };
  double positionBMparts53(int i){ return position_of_BM_parts_5_3[i]; };

  double rminBMparts54(void){ return Rmin_BM_parts_5_4; };
  double rmaxBMparts54(void){ return Rmax_BM_parts_5_4; };
  double lBMparts54(void){ return L_BM_parts_5_4; };
  double positionBMparts54(int i){ return position_of_BM_parts_5_4[i]; };
  double rotationBMparts54(int i){ return rotation_of_BM_parts_5_4[i]; };
   
  // York & Coil
  double sizeBMYork(int i, int j ){ return size_of_BMYork[i][j]; };
  double positionBMYork(int i, int j ){ return position_of_BMYork[i][j]; };

  double sizeBMCoilOuter(int i ){ return size_of_BMCoil_Outer[i]; };
  double sizeBMCoilInner(int i ){ return size_of_BMCoil_Inner[i]; };
  double positionBMCoil(int i, int j ){ return position_of_BMCoil[i][j]; };

  // Return BeamLine Component ( SLIT  )
  // parts-0 ( Frange )
  double rminSLITff( void ){ return RMin_SLIT_FF; }
  double rmaxSLITff( void ){ return RMax_SLIT_FF; }
  double lSLITff( void ){ return L_SLIT_FF; }
  double positionSLITff( int i ){  return position_SLIT_FF[i]; };
  double rminSLITbf(void){ return RMin_SLIT_BF; }
  double rmaxSLITbf(void){ return RMax_SLIT_BF; }
  double lSLITbf(void){ return L_SLIT_BF; }                       
  double positionSLITbf( int i ){ return position_SLIT_BF[i]; };
  double rotationSLITf( int i ){ return rotation_SLIT_F[i]; }
  // parts-1
  double rminSLITparts1(void){ return Rmin_SLIT_parts1; }
  double rmaxSLITparts1(void){ return Rmax_SLIT_parts1; }
  double lSLITparts1(void){ return L_SLIT_parts1; }
  double positionSLITparts1(int i){  return position_SLIT_parts1[i]; };
  // parts-2
  double rminSLITparts2(void){ return Rmin_SLIT_parts2; }
  double rmaxSLITparts2(void){ return Rmax_SLIT_parts2; }
  double lSLITparts2(void){ return L_SLIT_parts2; }
  double positionSLITparts2(int i){ return position_SLIT_parts2[i]; };
  // parts-3
  double rminSLITparts3(void){ return Rmin_SLIT_parts3; }
  double rmaxSLITparts3(void){ return Rmax_SLIT_parts3; }
  double lSLITparts3(void){ return L_SLIT_parts3; }
  double positionSLITparts3(int i){ return position_SLIT_parts3[i]; };
  double rotationSLITcomp(int i){ return rotation_SLIT_comp[i]; } 
  double rotationSLITduct(int i){ return rotation_SLIT_duct[i]; }
  double rotationSLITbody(int i){ return rotation_SLIT_body[i]; }

  // collimeter
  double sizeCollimeter_Vac(int i){ return size_of_collimeter_vac[i]; }
  double positionCollimeter_Vac(int i, int j){ return position_of_collimeter_vac[i][j]; };

  double sizeCollimeter(int i){ return size_of_collimeter[i]; }
  double positionCollimeter(int i, int j){ return position_of_collimeter[i][j]; };

  // Blank Frange
  double rminBlank(void){ return RMin_Blank; };
  double rmaxBlank(void){ return RMax_Blank; };
  double lBlank(void){ return L_Blank; };
  double positionBlank(int i){ return position_Blank[i]; };
  double rotationBlank(int i){ return rotation_Blank[i]; };

  // Beam Window 
  double rminWindow(void){ return RMin_Window; };
  double rmaxWindow(void){ return RMax_Window; };
  double lWindow(void){ return L_Window; };
  double positionWindow(int i){ return position_Window[i]; };
  double rotationWindow(int i){ return rotation_Window[i]; };

  double rminTiWindow(void){ return RMin_TiWindow; };
  double rmaxTiWindow(void){ return RMax_TiWindow; };
  double lTiWindow(void){ return L_TiWindow; };
  double positionTiWindow(int i){ return position_TiWindow[i]; };
  double rotationTiWindow(int i){ return rotation_TiWindow[i]; };

  double lWindow_vac(void){ return L_Window_vac; };
  double position_z_Window_vac(void){ return position_Window_vac_z; };

  // Top Plate
  double sizeTopPlate( int i ){ return size_TopPlate[i]; };
  double positionTopPlate(int i){ return position_TopPlate[i]; };
  double rotationTopPlate(int i){ return rotation_TopPlate[i]; };
 
  double rmin_TopPlate_hole(void){ return RMin_TopPlate_hole; }
  double rmax_TopPlate_hole(void){ return RMax_TopPlate_hole; }
  double l_TopPlate_hole(void){ return L_TopPlate_hole; }
  double positionTopPlate_hole(int i){ return position_TopPlate_hole[i]; };
  double rotationTopPlate_hole(int i){ return rotation_TopPlate_hole[i]; };

  // Faraday Cup
  // Al Cylinder
  // parts 1-3
  double rminFDAlCylinder(int i){return RMin_FDAlCylinder[i]; };
  double rmaxFDAlCylinder(int i){return RMax_FDAlCylinder[i]; };
  double lFDAlCylinder(int i){return L_FDAlCylinder[i]; };
  double positionFDAlCylinder(int i, int j){ return position_FDAlCylinder[i][j]; };
  double rotationFDAlCylinder(int i, int j){ return rotation_FDAlCylinder[i][j]; };

  // FD Body
  // parts 1-3 ( Pb ), parts 4 ( C )
  double rminFDBody(int i){return RMin_FDBody[i]; };
  double rmaxFDBody(int i){return RMax_FDBody[i]; };
  double lFDBody(int i){return L_FDBody[i]; };
  double positionFDBody(int i, int j){ return position_FDBody[i][j]; };
  double rotationFDBody(int i, int j){ return rotation_FDBody[i][j]; };

  // FC2-> FC4 in 2012.05.01
  // Front : Carbon 
  //       --> Front1/2 thin shield ( 1 and 2 ) Cu or Al or Any other material ?  in 2012.05.01
  // Creat Faraday Cup 4 in 2012.10.02  
  // Geometry was fixed except position                      

   // Ground Layer (G4Polycone)
  int numZPlanesGNDLayerFC4(void){ return NumZPlanesGNDLayerFC4;}
  double irGNDLayerFC4(int i){ return IrGNDLayerFC4[i]; }
  double orGNDLayerFC4(int i){ return OrGNDLayerFC4[i]; }
  double zGNDLayerFC4(int i){ return ZGNDLayerFC4[i]; }
  double positionGNDLayerFC4(int i){ return position_GNDLayerFC4[i]; }
  double rotationGNDLayerFC4(int i){ return rotation_GNDLayerFC4[i]; }

   // Shield Layer (G4Polycone)
  int numZPlanesShieldLayerFC4(void){ return NumZPlanesShieldLayerFC4;}
  double irShieldLayerFC4(int i){ return IrShieldLayerFC4[i]; }
  double orShieldLayerFC4(int i){ return OrShieldLayerFC4[i]; }
  double zShieldLayerFC4(int i){ return ZShieldLayerFC4[i]; }
  double positionShieldLayerFC4(int i){ return position_ShieldLayerFC4[i]; }
  double rotationShieldLayerFC4(int i){ return rotation_ShieldLayerFC4[i]; }

   // Body : Copper (G4Polycone)
  int numZPlanesBodyFC4(void){ return NumZPlanesBodyFC4;}
  double irBodyFC4(int i){ return IrBodyFC4[i]; }
  double orBodyFC4(int i){ return OrBodyFC4[i]; }
  double zBodyFC4(int i){ return ZBodyFC4[i]; }
  double positionBodyFC4(int i){ return position_BodyFC4[i]; }
  double rotationBodyFC4(int i){ return rotation_BodyFC4[i]; }
  
  // Isolator-1 : Macor   
  double rminIso1FC4(void){ return RMin_Iso1FC4; };
  double rmaxIso1FC4(void){ return RMax_Iso1FC4; };
  double lIso1FC4(void){ return L_Iso1FC4; };
  double positionIso1FC4(int i){ return position_Iso1FC4[i]; };
  double rotationIso1FC4(int i){ return rotation_Iso1FC4[i]; };

  // Isolator-2 : Macor 
  double rminIso2FC4(void){ return RMin_Iso2FC4; };
  double rmaxIso2FC4(void){ return RMax_Iso2FC4; };
  double lIso2FC4(void){ return L_Iso2FC4; };
  double positionIso2FC4(int i){ return position_Iso2FC4[i]; };
  double rotationIso2FC4(int i){ return rotation_Iso2FC4[i]; };

  // Top Plate : A5052 
  double sizeofTopPlateFC4(int i){ return size_of_TopPlateFC4[i]; }
  double positionTopPlateFC4(int i){ return position_TopPlateFC4[i]; }
  double rotationTopPlateFC4(int i){ return rotation_TopPlateFC4[i]; }

  //------------------------------------------------------
  // Faraday Cup 5 
  // made by B.K.Shin-san, in 2013.04.23 
  //------------------------------------------------------
  // FC 5 add at 130405
  // body : Copper (G4 Polycon) No 05
  int numZPlanesBodyFC5(void){ return NumZPlanesBodyFC5;}
  double irBodyFC5(int i){ return IrBodyFC5[i]; }
  double orBodyFC5(int i){ return OrBodyFC5[i]; }
  double zBodyFC5(int i){ return ZBodyFC5[i]; }
  double positionBodyFC5(int i){ return position_BodyFC5[i]; }
  double rotationBodyFC5(int i){ return rotation_BodyFC5[i]; }

  // FC 5 add at 130405
  // Barrel Suppoter No 04
  int numZPlanesBSFC5(void){ return NumZPlanesBSFC5;}
  double irBSFC5(int i){ return IrBSFC5[i]; }
  double orBSFC5(int i){ return OrBSFC5[i]; }
  double zBSFC5(int i){ return ZBSFC5[i]; }
  double positionBSFC5(int i){ return position_BSFC5[i]; }
  double rotationBSFC5(int i){ return rotation_BSFC5[i]; }

  // Barrel Suppoter Vaccum  No 04
  double rmaxBSVFC5(){ return rMaxBSVFC5; }
  double rminBSVFC5(){ return rMinBSVFC5; }
  double hBSVFC5(){ return HBSVFC5; }
  double positionBSVFC5(int i){ return position_BSVFC5[i]; }
  double rotationBSVFC5(int i){ return rotation_BSVFC5[i]; }

  // Barrel Suppoter Vaccum bottom  No 04
  double rmaxBSVBFC5(){ return rMaxBSVBFC5; }
  double rminBSVBFC5(){ return rMinBSVBFC5; }
  double hBSVBFC5(){ return HBSVBFC5; }
  double positionBSVBFC5(int i){ return position_BSVBFC5[i]; }
  double rotationBSVBFC5(int i){ return rotation_BSVBFC5[i]; }

  // Titan shield  No 07
  int numZPlanesTSFC5(void){ return NumZPlanesTSFC5;}
  double irTSFC5(int i){ return IrTSFC5[i]; }
  double orTSFC5(int i){ return OrTSFC5[i]; }
  double zTSFC5(int i){ return ZTSFC5[i]; }
  double positionTSFC5(int i){ return position_TSFC5[i]; }
  double rotationTSFC5(int i){ return rotation_TSFC5[i]; }

  // Titan Sheild Vaccum   No 07
  double brmaxTSVFC5(){ return brMaxTSVFC5; }
  double brminTSVFC5(){ return brMinTSVFC5; }
  double trmaxTSVFC5(){ return trMaxTSVFC5; }
  double trminTSVFC5(){ return trMinTSVFC5; }
  double hTSVFC5(){ return HTSVFC5; }
  double positionTSVFC5(int i){ return position_TSVFC5[i]; }
  double rotationTSVFC5(int i){ return rotation_TSVFC5[i]; }

  // Bottom Supporter No 06 (BotS)
  int numZPlanesBotSFC5(void){ return NumZPlanesBotSFC5;}
  double irBotSFC5(int i){ return IrBotSFC5[i]; }
  double orBotSFC5(int i){ return OrBotSFC5[i]; }
  double zBotSFC5(int i){ return ZBotSFC5[i]; }
  double positionBotSFC5(int i){ return position_BotSFC5[i]; }
  double rotationBotSFC5(int i){ return rotation_BotSFC5[i]; }

  // Top Sorppoter No 15 (TopS)
  int numZPlanesTopSFC5(void){ return NumZPlanesTopSFC5;}
  double irTopSFC5(int i){ return IrTopSFC5[i]; }
  double orTopSFC5(int i){ return OrTopSFC5[i]; }
  double zTopSFC5(int i){ return ZTopSFC5[i]; }
  double positionTopSFC5(int i){ return position_TopSFC5[i]; }
  double rotationTopSFC5(int i){ return rotation_TopSFC5[i]; }

  // Top Suppoter Vaccum bottom
  double rmaxTopSVFC5(){ return rMaxTopSVFC5; }
  double rminTopSVFC5(){ return rMinTopSVFC5; }
  double hTopSVFC5(){ return HTopSVFC5; }
  double positionTopSVFC5(int i){ return position_TopSVFC5[i]; }
  double rotationTopSVFC5(int i){ return rotation_TopSVFC5[i]; }

  // CF90 Frange  No  02  (CF90)
  int numZPlanesCF90FC5(void){ return NumZPlanesCF90FC5;}
  double irCF90FC5(int i){ return IrCF90FC5[i]; }
  double orCF90FC5(int i){ return OrCF90FC5[i]; }
  double zCF90FC5(int i){ return ZCF90FC5[i]; }
  double positionCF90FC5(int i){ return position_CF90FC5[i]; }
  double rotationCF90FC5(int i){ return rotation_CF90FC5[i]; }

  // Top Surpporter  No 06b (bTopS)
  int numZPlanesbTopSFC5(void){ return NumZPlanesbTopSFC5;}
  double irbTopSFC5(int i){ return IrbTopSFC5[i]; }
  double orbTopSFC5(int i){ return OrbTopSFC5[i]; }
  double zbTopSFC5(int i){ return ZbTopSFC5[i]; }
  double positionbTopSFC5(int i){ return position_bTopSFC5[i]; }
  double rotationbTopSFC5(int i){ return rotation_bTopSFC5[i]; }

  // Copper Chamber No 01 (CP)
  int numZPlanesCPFC5(void){ return NumZPlanesCPFC5;}
  double irCPFC5(int i){ return IrCPFC5[i]; }
  double orCPFC5(int i){ return OrCPFC5[i]; }
  double zCPFC5(int i){ return ZCPFC5[i]; }
  double positionCPFC5(int i){ return position_CPFC5[i]; }
  double rotationCPFC5(int i){ return rotation_CPFC5[i]; }

  // Top Plate (TP)
  double xTPFC5(){ return XTPFC5; }
  double yTPFC5(){  return YTPFC5; }
  double zTPFC5(){  return ZTPFC5; }
  double positionTPFC5(int i){ return position_TPFC5[i]; }
  double rotationTPFC5(int i){ return rotation_TPFC5[i]; }

  // Feed Through
  double rmaxFTFC5(){ return rMaxFTFC5; }
  double rminFTFC5(){ return rMinFTFC5; }
  double hFTFC5(){ return HFTFC5; }
  double positionFTFC5(int i){ return position_FTFC5[i]; }
  double rotationFTFC5(int i){ return rotation_FTFC5[i]; }

  // Feed Through 2
  double rmaxFT2FC5(){ return rMaxFT2FC5; }
  double rminFT2FC5(){ return rMinFT2FC5; }
  double hFT2FC5(){ return HFT2FC5; }
  double positionFT2FC5(int i){ return position_FT2FC5[i]; }
  double rotationFT2FC5(int i){ return rotation_FT2FC5[i]; }

  // Vaccum duct
  double rmaxVDFC5(){ return rMaxVDFC5; }
  double rminVDFC5(){ return rMinVDFC5; }
  double hVDFC5(){ return HVDFC5; }
  double positionVDFC5(int i){ return position_VDFC5[i]; }
  double rotationVDFC5(int i){ return rotation_VDFC5[i]; }

  // Vaccum duct Vacuum
  double rmaxVDVFC5(){ return rMaxVDVFC5; }
  double rminVDVFC5(){ return rMinVDVFC5; }
  double hVDVFC5(){ return HVDVFC5; }
  double positionVDVFC5(int i){ return position_VDVFC5[i]; }
  double rotationVDVFC5(int i){ return rotation_VDVFC5[i]; }

  //-----------------------------------------------------


  // Screen Monitor 3
  double sizeScreenMonitor3(int i){ return size_of_ScreenMonitor3[i]; };
  double positionScreenMonitor3(int i){  return position_ScreenMonitor3[i]; };
  double rotationScreenMonitor3(int i){  return rotation_ScreenMonitor3[i]; };

  // Beam Attenuator
  double rminBeamAttenuator(void){ return RMin_BeamAttenuator; }
  double rmaxBeamAttenuator(void){ return RMax_BeamAttenuator; }
  double lBeamAttenuator(void){ return L_BeamAttenuator; }
  double positionBeamAttenuator(int i){ return position_BeamAttenuator[i]; }
  double rotationBeamAttenuator(int i){ return rotation_BeamAttenuator[i]; }  

  // Pb Collimator added in 2012.09.24
  double rminPbCollimator(void){ return RMin_PbCollimator; }
  double rmaxPbCollimator(void){ return RMax_PbCollimator; }
  double lPbCollimator(void){ return L_PbCollimator; }
  double positionPbCollimator(int i){ return position_PbCollimator[i]; }
  double rotationPbCollimator(int i){ return rotation_PbCollimator[i]; }
  
  // Virtual Test Chamber added in 2011.12.13
   // Cylinder( Material ) parts 1-5
  double rminVirtualChamberCylinder(int i){return RMin_VirtualChamberCylinder[i]; };
  double rmaxVirtualChamberCylinder(int i){return RMax_VirtualChamberCylinder[i]; };
  double lVirtualChamberCylinder(int i){return L_VirtualChamberCylinder[i]; };
  double positionVirtualChamberCylinder(int i, int j){ return position_VirtualChamberCylinder[i][j]; };
  double rotationVirtualChamberCylinder(int i, int j){ return rotation_VirtualChamberCylinder[i][j]; };

  // Target : added in 2013.04.17, to check backscattering effect
  double rminVirtualChamberTarget(void){ return RMin_VirtualChamberTarget; }
  double rmaxVirtualChamberTarget(void){ return RMax_VirtualChamberTarget; }
  double lVirtualChamberTarget(void){ return L_VirtualChamberTarget;       }
  double positionVirtualChamberTarget(int i){ return position_VirtualChamberTarget[i]; }
  double rotationVirtualChamberTarget(int i){ return rotation_VirtualChamberTarget[i]; }

   // Body : Inner side ( air ) parts 1-4
  double rminVirtualChamberBody(int i){ return RMin_VirtualChamberBody[i]; };
  double rmaxVirtualChamberBody(int i){ return RMax_VirtualChamberBody[i]; };
  double lVirtualChamberBody(int i){ return L_VirtualChamberBody[i]; };
  double positionVirtualChamberBody(int i, int j){ return position_VirtualChamberBody[i][j]; };
  double rotationVirtualChamberBody(int i, int j){ return rotation_VirtualChamberBody[i][j]; };  

  // ICE methods
  double positionIce(int i){ return ice_position[i]; };
  double sizeIce(int i){ return ice_size[i]; };
  double rotationIce(int i){ return ice_rotation[i]; };


  // Magnet Field
  double BeamLinePoistion(int i){ return position_of_beam_start[i]; };
  double QuadrupoleManget1FiledResion(int i, int j){ return quadrupole_magnet1_field_region[i][j]; };
  double QuadrupoleManget2FiledResion(int i, int j){ return quadrupole_magnet2_field_region[i][j]; };

  double BendingMangetFiledResion(int i, int j){ return bending_magnet_field_region[i][j]; };

  // Beam Injection Position
  double Beam_injection_position(int i){ return beam_injection_position[i]; }

  //---------------------------------------
  // Pb Block & Concrete Block
  //---------------------------------------

  void shielding_parameter_initialize();

  // Pb Tower
  double sizePbTower(int i, int j){ return size_of_PbTower[i][j]; };
  double positionPbTower(int i, int j){ return position_of_PbTower[i][j]; };

  // BM Duct Pb Block
  // parts-1 
  double rminPbBlockBMparts1(int i){ return Rmin_PbBlock_BM_parts1[i]; };
  double rmaxPbBlockBMparts1(int i){ return Rmax_PbBlock_BM_parts1[i]; };
  double lPbBlockparts1(int i){ return L_PbBlock_BM_parts1[i]; };
  double phistartPbBlockBMparts1(int i){ return Phi_start_PbBlock_BM_parts1[i]; };
  double phideltaPbBlockBMparts1(int i){ return Phi_delta_PbBlock_BM_parts1[i]; };
  double positionPbBlockBMparts1(int i, int j){ return position_of_PbBlock_BM_parts1[i][j]; };
  double rotationPbBlockBMparts1(int i, int j){ return rotation_of_PbBlock_BM_parts1[i][j]; };
  // parts-2
  double sizePbBlockBMparts2(int i, int j){ return size_of_PbBlock_BM_parts2[i][j]; };
  double positionPbBlockBMparts2(int i, int j){ return position_of_PbBlock_BM_parts2[i][j]; };
  // parts-3
  double sizePbBlockBMparts3(int i, int j){ return size_of_PbBlock_BM_parts3[i][j]; };
  double positionPbBlockBMparts3(int i, int j){ return position_of_PbBlock_BM_parts3[i][j]; };
  // Concrete Blocks
  int numCon(void){return numConcrete; };
  double sizeConcreteBlock(int i, int j){ return size_of_ConcreteBlock[i][j]; };
  double positionConcreteBlock(int i, int j){ return position_of_ConcreteBlock[i][j]; };

  //---------------------------------------
  // Detector Planes
  //---------------------------------------
  void detectorplane_parameter_initialize();

  int numDP(void){return numDetectorPlane; };
  double sizeDetectorPlane(int i, int j){ return size_of_DetectorPlane[i][j]; };
  double positionDetectorPlane(int i, int j){ return position_of_DetectorPlane[i][j]; };

  //DUMP
  void dump_ELS_geometory(void);

private:

  //-------------------------------------------------
  // ELS Fence &
    // Concrete pad &
  //  ELS ( 40ft HighCube ) Container // unit = mm
  // ------------------------------------------------

  // Geometory
  double size_of_concretepad[3];
  double distance_concretepad_ELS[3];

  /* This parameter means the bias height of ELS site ground 
     From survey results */
  double elssite_ground_height;
  
  double outer_size_of_ELScontainer[3];
  double inner_size_of_ELScontainer[3];
 
  double thick_of_ELScontainer_wall_forward;
  double thick_of_ELScontainer_wall_backward;
  double thick_of_ELScontainer_wall_side;
  double thick_of_ELScontainer_wall_roof;
  double thick_of_ELScontainer_wall_bottom;

  double position_bias_ELScontainer_wall_forward;
  double position_bias_ELScontainer_wall_backward;
  double position_bias_ELScontainer_wall_side;

  // Position & Rotation
  double position_of_concretepad[3];
  double rotation_of_concretepad[3];

  double position_of_ELS[3];
   double position_of_ELS_outer[3];
   double position_of_ELS_inner[3];
  double rotation_of_ELS[3];
   double rotation_of_ELS_outer[3];
   double rotation_of_ELS_inner[3];

  double size_of_ELS_injection_hole[3];
  double position_of_ELS_injection_hole[3];

  //----------------------------------------
  // Cover Box
  //----------------------------------------
  double size_of_coverbox_outer[3];
  double thickness_of_coverbox;
  double size_of_coverbox_inner[3];

  double thickness_of_coverbox_board;
  double size_of_coverbox_board[3];
  double rmin_coverbox_boardhole;
  double rmax_coverbox_boardhole;
  double l_coverbox_boardhole;

  double position_of_coverbox_outer[3];
  double position_of_coverbox_inner[3];

  double zbias_of_coverbox_board;
  double position_of_coverbox_board[3];
  double position_of_coverbox_boardhole[3];

  /* alignment of cover box hole 
     [0] ... shift along x-axis
     [1] ... shift along y-axis */
  double alignment_coverbox_boardhole_shift[2];

  /* alignment of vertical beam line 
     [0] ... shift along x-axis 
     [1] ... shift along y-axis */
  double alignment_vertical_beamline_shift[2];


  //-------------------------------------------------
  //  Base A & B & C
  //-------------------------------------------------

  // Geometory
  double size_of_baseA[3];
  double size_of_baseB[3];
  double size_of_baseC[3];
  
  double space_bwt_base_wall_forward;
  double space_bwt_base_wall_backward;
  double space_bwt_base_wall_side;

  // Position
  double position_of_baseA[3];
  double position_of_baseB[3];
  double position_of_baseC[3];

  double rotation_of_baseA[3];
  double rotation_of_baseB[3];
  double rotation_of_baseC[3];

  double size_of_base[3][3];
  double position_of_base[3][3];
  double rotation_of_base[3][3];

  //-------------------------------------------------
  //  Beam Line Component
  //-------------------------------------------------

  double position_of_beam_start[3];

  double Total_Beam_Length_H;
  double Total_Beam_Length_H2;
  double Total_Beam_Length_V;  

  // Frange Geometory

  // Frange ICF203
  double RMin_ICF203;
  double RMax_ICF203;
  double L_ICF203;

  // Frange ICF070(1)
  double RMin_ICF070_1;
  double RMax_ICF070_1;
  double L_ICF070_1;

  // Frange ICF070(2)
  double RMin_ICF070_2;
  double RMax_ICF070_2;
  double L_ICF070_2;

  // Frange V85
  double RMin_V85;
  double RMax_V85;
  double L_V85;

  // Beam Line Components 

  // No.1 EGUN
  double RMin_EGUN_FF;
  double RMax_EGUN_FF;
  double L_EGUN_FF;

  double RMin_EGUN;
  double RMax_EGUN;
  double L_EGUN;

  double RMin_EGUN_BF;
  double RMax_EGUN_BF;
  double L_EGUN_BF;

  double position_EGUN_FF[3];
  double position_EGUN[3];
  double position_EGUN_BF[3];
  double rotation_EGUN[3];

  //  Blank Frange    ICF203
  double RMin_EGUNBlank;
  double RMax_EGUNBlank;
  double L_EGUNBlank;

  double position_EGUNBlank[3];
  double rotation_EGUNBlank[3];
  
  //  Anode
  double RMin_Anode;
  double RMax_Anode;
  double L_Anode;

  double position_Anode[3];
  double rotation_Anode[3];

  // No.2 ML         G4Tube     Fragne ICF203-ICF070
  double RMin_ML_FF;
  double RMax_ML_FF;
  double L_ML_FF;

  double RMin_ML;
  double RMax_ML;
  double L_ML;

  double RMin_ML_BF;
  double RMax_ML_BF;
  double L_ML_BF;

  double position_ML_FF[3];
  double position_ML[3];
  double position_ML_BF[3];
  double rotation_ML[3]; 

  // No.3 ICF070-GV  G4Tube     Fragne ICF070-ICF070
  double RMin_ICF070GV_FF;
  double RMax_ICF070GV_FF;
  double L_ICF070GV_FF;

  double RMin_ICF070GV;
  double RMax_ICF070GV;
  double L_ICF070GV;

  double RMin_ICF070GV_BF;
  double RMax_ICF070GV_BF;
  double L_ICF070GV_BF;

  double position_ICF070GV_FF[3];
  double position_ICF070GV[3];
  double position_ICF070GV_BF[3];
  double rotation_ICF070GV[3]; 

  // No.4 ICF070-V85 duct G4Tube   Fragne ICF070-V85
  double RMin_ICF070V85Duct_FF;
  double RMax_ICF070V85Duct_FF;
  double L_ICF070V85Duct_FF;

  double RMin_ICF070V85Duct;
  double RMax_ICF070V85Duct;
  double L_ICF070V85Duct;

  double RMin_ICF070V85Duct_BF;
  double RMax_ICF070V85Duct_BF;
  double L_ICF070V85Duct_BF;

  double position_ICF070V85Duct_FF[3];
  double position_ICF070V85Duct[3];
  double position_ICF070V85Duct_BF[3];
  double rotation_ICF070V85Duct[3];

  // No.5 CoreMonitor1 G4Tube      Fragne V85-V85
  double RMin_CM1_FF;
  double RMax_CM1_FF;
  double L_CM1_FF;

  double RMin_CM1;
  double RMax_CM1;
  double L_CM1;

  double RMin_CM1_BF;
  double RMax_CM1_BF;
  double L_CM1_BF;

  double position_CM1_FF[3];
  double position_CM1[3];
  double position_CM1_BF[3];
  double rotation_CM1[3];

  // No.6 PB+B-Tube G4Tube         Fragne V85-V85
  double RMin_PBBTube_FF;
  double RMax_PBBTube_FF;
  double L_PBBTube_FF;

  double RMin_PBBTube;
  double RMax_PBBTube;
  double L_PBBTube;

  double RMin_PBBTube_BF;
  double RMax_PBBTube_BF;
  double L_PBBTube_BF;

  double position_PBBTube_FF[3];
  double position_PBBTube[3];
  double position_PBBTube_BF[3];
  double rotation_PBBTube[3];

  // No.7 CoreMonitor2 G4Tube      Fragne V85-V85
  double RMin_CM2_FF;
  double RMax_CM2_FF;
  double L_CM2_FF;

  double RMin_CM2;
  double RMax_CM2;
  double L_CM2;

  double RMin_CM2_BF;
  double RMax_CM2_BF;
  double L_CM2_BF;

  double position_CM2_FF[3];
  double position_CM2[3];
  double position_CM2_BF[3];
  double rotation_CM2[3];

  // No.8 2m-Tube G4Tube           Fragne V85-V85
  double RMin_2mTube_FF;
  double RMax_2mTube_FF;
  double L_2mTube_FF;

  double RMin_2mTube;
  double RMax_2mTube;
  double L_2mTube;

  double RMin_2mTube_BF;
  double RMax_2mTube_BF;
  double L_2mTube_BF;

  double position_2mTube_FF[3];
  double position_2mTube[3];
  double position_2mTube_BF[3];
  double rotation_2mTube[3];

  // No.9 Drift Duct G4Tube        Fragne V85-V85
  double RMin_V85V85Duct_FF;
  double RMax_V85V85Duct_FF;
  double L_V85V85Duct_FF;

  double RMin_V85V85Duct;
  double RMax_V85V85Duct;
  double L_V85V85Duct;

  double RMin_V85V85Duct_BF;
  double RMax_V85V85Duct_BF;
  double L_V85V85Duct_BF;

  double position_V85V85Duct_FF[3];
  double position_V85V85Duct[3];
  double position_V85V85Duct_BF[3];
  double rotation_V85V85Duct[3];

  // No.10 V85-GV G4Tube           Fragne V85-V85
  double RMin_V85GV_FF;
  double RMax_V85GV_FF;
  double L_V85GV_FF;

  double RMin_V85GV;
  double RMax_V85GV;
  double L_V85GV;

  double RMin_V85GV_BF;
  double RMax_V85GV_BF;
  double L_V85GV_BF;

  double position_V85GV_FF[3];
  double position_V85GV[3];
  double position_V85GV_BF[3];
  double rotation_V85GV[3];

  // No.11 QM G4Tube               Fragne V85-V85
  double RMin_QM_FF;
  double RMax_QM_FF;
  double L_QM_FF;

  double RMin_QM;
  double RMax_QM;
  double L_QM;

  double RMin_QM_BF;
  double RMax_QM_BF;
  double L_QM_BF;

  double position_QM_FF[3];
  double position_QM[3];
  double position_QM_BF[3];
  double rotation_QM[3];

  // No.12 ScreenMonitor1 Fragne V85-V85
  double L_SM1;
  // parts-0 ( Frange )
  double RMin_SM1_FF;
  double RMax_SM1_FF;
  double L_SM1_FF;
  double RMin_SM1_BF;
  double RMax_SM1_BF;
  double L_SM1_BF; 
  double position_SM1_FF[3];
  double position_SM1_BF[3];
  double rotation_SM1_F[3];
  // parts-1
  double Rmin_SM1_parts1;
  double Rmax_SM1_parts1;
  double L_SM1_parts1;
  double position_SM1_parts1[3];
  // parts-2
  double Rmin_SM1_parts2;
  double Rmax_SM1_parts2;
  double L_SM1_parts2;
  double position_SM1_parts2[3];
  // parts-3
  double Rmin_SM1_parts3;
  double Rmax_SM1_parts3;
  double L_SM1_parts3;
  double position_SM1_parts3[3];

  double rotation_SM1_comp[3];
  double rotation_SM1_duct[3];
  double rotation_SM1_body[3];

  // No.13 CoreMonitor3 G4Tube     Fragne V85-V85
  double RMin_CM3_FF;
  double RMax_CM3_FF;
  double L_CM3_FF;

  double RMin_CM3;
  double RMax_CM3;
  double L_CM3;

  double RMin_CM3_BF;
  double RMax_CM3_BF;
  double L_CM3_BF;

  double position_CM3_FF[3];
  double position_CM3[3];
  double position_CM3_BF[3];
  double rotation_CM3[3];

  // No.14 BM                      Fragne V85-V85-V85
  // Total ... 4 + 9 + 3 = 16 parts 

  // No.14.5 Cu-Collimeter G4         Fragne V85-V85
  double L_CuCollimeter;
  double W_CuCollimeter;
  double H_CuCollimeter;
  double size_of_CuCollimeter[3];
  double position_CuCollimeter[3];
  double rotation_CuCollimeter[3];

  double RMin_CuCollimeter_hole;
  double RMax_CuCollimeter_hole; 
  double L_CuCollimeter_hole;
  double position_CuCollimeter_hole[3];
  double rotation_CuCollimeter_hole[3];

     // alignment parameter in collimator hole
     // alignment is on only vertical axis(z)
     double alignment_position_CuCollimeter_hole[3];
          
  // parts-0 : BM-Duct Frange ( G4Subtraction : G4Tubs - G4Box )
  // Three Franges 
  // 0 : Insident frange
  // 1 : Straight out frange
  // 2 : Bending out frange
  double Rmin_V85_BM;
  double Rmax_V85_BM;
  double L_V85_BM;
  double size_of_V85_BM_IN[3];  
  double position_BM_F[3][3];
  double rotation_BM_F[3][3];

  double L1_BM;
  double LH_BM;
  double LV_BM;

  // parts-1 : Straight Duct ( Insident part )
  double size_of_BM_parts_1_outer[3];
  double size_of_BM_parts_1_inner[3];
  double position_of_BM_parts_1_outer[3];
  double position_of_BM_parts_1_inner[3];

  double size_of_BM_parts_1_vac[3];
  double position_of_BM_parts_1_vac[3];

  // parts-2 : Straight Duct ( Straight Outt part )
  double size_of_BM_parts_2_outer[3];
  double size_of_BM_parts_2_inner[3];
  double position_of_BM_parts_2_outer[3];
  double position_of_BM_parts_2_inner[3];

  double size_of_BM_parts_2_vac[3];
  double position_of_BM_parts_2_vac[3];

  // parts-3 : Straight Duct ( Bending Outt part )
  double size_of_BM_parts_3_outer[3];
  double size_of_BM_parts_3_inner[3];
  double position_of_BM_parts_3_outer[3];
  double position_of_BM_parts_3_inner[3];
 
  double size_of_BM_parts_3_vac[3];
  double position_of_BM_parts_3_vac[3];

  // parts-4(1) : Main Body ( Side Board left )
  double size_of_BM_parts_4_1[3];
  double position_of_BM_parts_4_1[3];
  
  // parts-4(2) : Main Body ( Side Board right )
  double size_of_BM_parts_4_2[3];
  double position_of_BM_parts_4_2[3];

  // parts-4(3) : Main Body ( Side Board Cut Volume 1 )
  double Rmin_BM_parts_4_3;
  double Rmax_BM_parts_4_3;
  double L_BM_parts_4_3;
  double position_of_BM_parts_4_3[3];
  double rotation_of_BM_parts_4_3[3];

  // parts-4(4) : Main Body ( Side Board Cut Volume 2 )
  double size_of_BM_parts_4_4[3];
  double position_of_BM_parts_4_4[3];

  // parts-4(5) : Main Body ( Side Board Cut Volume 3 )
  double Rmin_BM_parts_4_5;
  double Rmax_BM_parts_4_5;
  double L_BM_parts_4_5;
  double position_of_BM_parts_4_5[3];
  double rotation_of_BM_parts_4_5[3];
  
  // parts-4(6) : Main Body ( Side Board  )
  double Rmin_BM_parts_4_6;
  double Rmax_BM_parts_4_6;
  double L_BM_parts_4_6;
  double Phi_start_BM_parts_4_6;
  double Phi_end_BM_parts_4_6;
  double position_of_BM_parts_4_6[3];
  double rotation_of_BM_parts_4_6[3];

  // parts-4(7) : Main Body ( Side Board  )
  double size_of_BM_parts_4_7[3];
  double position_of_BM_parts_4_7[3];

  // parts-4(8) : Main Body ( Side Board  )
  double Rmin_BM_parts_4_8;
  double Rmax_BM_parts_4_8;
  double L_BM_parts_4_8;
  double Phi_start_BM_parts_4_8;
  double Phi_end_BM_parts_4_8;
  double position_of_BM_parts_4_8[3];
  double rotation_of_BM_parts_4_8[3];

  // parts-4(9) : Main Body ( Bottom Board  )
  double size_of_BM_parts_4_9[3];
  double position_of_BM_parts_4_9[3];

  //Vacuum region of BM-Duct
    // parts-5(1) : Main Body
    double size_of_BM_parts_5_1[3];
    double position_of_BM_parts_5_1[3]; 
     
    // parts-5(2) : Main Body ( Cut Volume 1 )
    double Rmin_BM_parts_5_2;
    double Rmax_BM_parts_5_2;
    double L_BM_parts_5_2;
    double position_of_BM_parts_5_2[3];
    double rotation_of_BM_parts_5_2[3];

    // parts-5(3) : Main Body ( Cut Volume 2 )
    double size_of_BM_parts_5_3[3];
    double position_of_BM_parts_5_3[3];

    // parts-5(4) : Main Body ( Cut Volume 3 )
    double Rmin_BM_parts_5_4;
    double Rmax_BM_parts_5_4;
    double L_BM_parts_5_4;
    double position_of_BM_parts_5_4[3];
    double rotation_of_BM_parts_5_4[3];

  // BM York & Coil
  double size_of_BMYork[5][3];
  double position_of_BMYork[5][3];

  double size_of_BMCoil_Outer[3];
  double size_of_BMCoil_Inner[3];
  double position_of_BMCoil[2][3];

     /* alignment of BM York & Coil
      0 : shift along x-axis
      1 : shift along y-axis
      2 : shift along z-axis
    */
     double  alignment_position_BMYork_Coil[3];


  // No.15 V85-GV2 G4Tube           Fragne V85-V85
  double RMin_V85GV2_FF;
  double RMax_V85GV2_FF;
  double L_V85GV2_FF;

  double RMin_V85GV2;
  double RMax_V85GV2;
  double L_V85GV2;

  double RMin_V85GV2_BF;
  double RMax_V85GV2_BF;
  double L_V85GV2_BF;

  double position_V85GV2_FF[3];
  double position_V85GV2[3];
  double position_V85GV2_BF[3];
  double rotation_V85GV2[3];

  // No.16 T-Duct G4Tube           Fragne V85-V85
  double RMin_TDuct_FF;
  double RMax_TDuct_FF;
  double L_TDuct_FF;

  double RMin_TDuct;
  double RMax_TDuct;
  double L_TDuct;

  double RMin_TDuct_BF;
  double RMax_TDuct_BF;
  double L_TDuct_BF;

  double position_TDuct_FF[3];
  double position_TDuct[3];
  double position_TDuct_BF[3];
  double rotation_TDuct[3];

  // No.17 Blank Frange            Fragne V85-Blank
  double RMin_Blank;
  double RMax_Blank;
  double L_Blank;

  double position_Blank[3];
  double rotation_Blank[3];

  //==================================
  // Vertical Beam Line
  //==================================

  // No.18 DriftTube G4Tube        Fragne V85-V85
  double RMin_V85V85Duct2_FF;
  double RMax_V85V85Duct2_FF;
  double L_V85V85Duct2_FF;

  double RMin_V85V85Duct2;
  double RMax_V85V85Duct2;
  double L_V85V85Duct2;

  double RMin_V85V85Duct2_BF;
  double RMax_V85V85Duct2_BF;
  double L_V85V85Duct2_BF;

  double position_V85V85Duct2_FF[3];
  double position_V85V85Duct2[3];
  double position_V85V85Duct2_BF[3];
  double rotation_V85V85Duct2[3];

  // No.19 Slit                    Fragne V85-V85
  double L_SLIT;
  // parts-0 ( Frange )
  double RMin_SLIT_FF;
  double RMax_SLIT_FF;
  double L_SLIT_FF;
  double RMin_SLIT_BF;
  double RMax_SLIT_BF;
  double L_SLIT_BF; 
  double position_SLIT_FF[3];
  double position_SLIT_BF[3];
  double rotation_SLIT_F[3];
  // parts-1
  double Rmin_SLIT_parts1;
  double Rmax_SLIT_parts1;
  double L_SLIT_parts1;
  double position_SLIT_parts1[3];
  // parts-2
  double Rmin_SLIT_parts2;
  double Rmax_SLIT_parts2;
  double L_SLIT_parts2;
  double position_SLIT_parts2[3];
  // parts-3
  double Rmin_SLIT_parts3;
  double Rmax_SLIT_parts3;
  double L_SLIT_parts3;
  double position_SLIT_parts3[3];

  double rotation_SLIT_comp[3];
  double rotation_SLIT_duct[3];
  double rotation_SLIT_body[3];

  // Ta Collimeter
  double size_of_collimeter_vac[3];
  double position_of_collimeter_vac[2][3];

  double size_of_collimeter[3];
  double position_of_collimeter[2][3];

     // alignment parameter in slit
     //   collimators position ( we can shift right collimator and 
     //                             left collimator respectively. ) 
     //  shift is only along x-axis ... straight beam line 
     /*
              vertical beam line ( z-axis )
            FD  <--     ^    --> 
           East         |         West
                        |
                    [1] | [0]...collimator(Ta)    
                        | 
                        |                  
        ---->---------->    --> horizontal beam direction (X-axis)                               
      */
     double slit_width;
     double alignment_position_Slit_Collimeter[2];

  // No.20 Drift Duct G4Tube       Fragne V85-V85
  double RMin_V85V85Duct3_FF;
  double RMax_V85V85Duct3_FF;
  double L_V85V85Duct3_FF;

  double RMin_V85V85Duct3;
  double RMax_V85V85Duct3;
  double L_V85V85Duct3;

  double RMin_V85V85Duct3_BF;
  double RMax_V85V85Duct3_BF;
  double L_V85V85Duct3_BF;

  double position_V85V85Duct3_FF[3];
  double position_V85V85Duct3[3];
  double position_V85V85Duct3_BF[3];
  double rotation_V85V85Duct3[3];

  // No.22 Screen Monitor 2     Fragne V85-V85
  double L_SM2;
  // parts-0 ( Frange )
  double RMin_SM2_FF;
  double RMax_SM2_FF;
  double L_SM2_FF;
  double RMin_SM2_BF;
  double RMax_SM2_BF;
  double L_SM2_BF; 
  double position_SM2_FF[3];
  double position_SM2_BF[3];
  double rotation_SM2_F[3];
  // parts-1
  double Rmin_SM2_parts1;
  double Rmax_SM2_parts1;
  double L_SM2_parts1;
  double position_SM2_parts1[3];
  // parts-2
  double Rmin_SM2_parts2;
  double Rmax_SM2_parts2;
  double L_SM2_parts2;
  double position_SM2_parts2[3];
  // parts-3
  double Rmin_SM2_parts3;
  double Rmax_SM2_parts3;
  double L_SM2_parts3;
  double position_SM2_parts3[3];

  double rotation_SM2_comp[3];
  double rotation_SM2_duct[3];
  double rotation_SM2_body[3];

  // No.22 CoreMonitor4 G4Tube     Fragne V85-V85
  double RMin_CM4_FF;
  double RMax_CM4_FF;
  double L_CM4_FF;

  double RMin_CM4;
  double RMax_CM4;
  double L_CM4;

  double RMin_CM4_BF;
  double RMax_CM4_BF;
  double L_CM4_BF;

  double position_CM4_FF[3];
  double position_CM4[3];
  double position_CM4_BF[3];
  double rotation_CM4[3];

  // No.23 Beam Window G4Tube      Fragne V85
  double RMin_Window;
  double RMax_Window;
  double L_Window;

  double position_Window[3];
  double rotation_Window[3];

  //  No.23.5 TiWindow
  double RMin_TiWindow;
  double RMax_TiWindow;
  double L_TiWindow;

  double position_TiWindow[3];
  double rotation_TiWindow[3];

  double L_Window_vac;
  double position_Window_vac_z;

  // No.23.6 Top Plate
  double size_TopPlate[3];
  double position_TopPlate[3];
  double rotation_TopPlate[3];

  double RMin_TopPlate_hole;
  double RMax_TopPlate_hole;
  double L_TopPlate_hole;
  double position_TopPlate_hole[3];
  double rotation_TopPlate_hole[3];

  // No.24 Faraday Cup
  double distance_fc;  
  double first_position_v;

  // Al Cylinder 
  // parts 1-3
  double RMin_FDAlCylinder[3];
  double RMax_FDAlCylinder[3];
  double L_FDAlCylinder[3];
  double position_FDAlCylinder[3][3];
  double rotation_FDAlCylinder[3][3];

  // FD Body
  // parts 1-3 ( Pb ), parts 4 ( C )   
  double delta_fc;
  double RMin_FDBody[4];
  double RMax_FDBody[4];
  double L_FDBody[4];
  double position_FDBody[4][3];
  double rotation_FDBody[4][3];

  /* Faraday Cup Alignment */
  double alignment_position_Faradaycup[3];

  // FC2-> FC4 in 2012.05.01 
  // Front : Carbon
  //       --> Front1/2 thin shield ( 1 and 2 ) Cu or Al or Any other material ?  in 2012.05.01      
  // Creat Faraday Cup 4 in 2012.10.02                                                                                             
  // Geometry was fixed except position    

   // Ground Layer (G4Polycone)
   int NumZPlanesGNDLayerFC4;
   double IrGNDLayerFC4[8];
   double OrGNDLayerFC4[8];
   double ZGNDLayerFC4[8];
   double position_GNDLayerFC4[3];
   double rotation_GNDLayerFC4[3];
  
   // Shield Layer (G4Polycone)
   int NumZPlanesShieldLayerFC4;
   double IrShieldLayerFC4[14];
   double OrShieldLayerFC4[14];
   double ZShieldLayerFC4[14];
   double position_ShieldLayerFC4[3];
   double rotation_ShieldLayerFC4[3];

   // Body : Copper (G4Polycone)
   int NumZPlanesBodyFC4;
   double IrBodyFC4[6];
   double OrBodyFC4[6];
   double ZBodyFC4[6];
   double position_BodyFC4[3];
   double rotation_BodyFC4[3];
  
  // Isolator-1 : Macor   
  double RMin_Iso1FC4;
  double RMax_Iso1FC4;
  double L_Iso1FC4;
  double position_Iso1FC4[3];
  double rotation_Iso1FC4[3];

  // Isolator-2 : Macor 
  double RMin_Iso2FC4;
  double RMax_Iso2FC4;
  double L_Iso2FC4;
  double position_Iso2FC4[3];
  double rotation_Iso2FC4[3];

  // Top Plate : A5052 
  double size_of_TopPlateFC4[3];
  double position_TopPlateFC4[3];
  double rotation_TopPlateFC4[3];


  //-------------------------------------------
  // Faraday Cup 5   
  // made by B.K.Shin-san, added in 2013.04.23
  //-------------------------------------------
  //Parts 5  Body : Copper (G4Polycone)
  int NumZPlanesBodyFC5;
  double IrBodyFC5[4];
  double OrBodyFC5[4];
  double ZBodyFC5[4];
  double position_BodyFC5[3];
  double rotation_BodyFC5[3];

  //Parts No 4  Berrel Supporter  (G4Polycone)
  int NumZPlanesBSFC5;
  double IrBSFC5[8];
  double OrBSFC5[8];
  double ZBSFC5[8];

  double position_BSFC5[3];
  double rotation_BSFC5[3];

  // Vaccum for Berrel Supporter (G4Polycone);
  double rMaxBSVFC5;
  double rMinBSVFC5;
  double HBSVFC5;
  double  position_BSVFC5[3];
  double  rotation_BSVFC5[3];

  // Vaccum for Berrel Supporter Bottom  (G4Polycone);
  double rMaxBSVBFC5;
  double rMinBSVBFC5;
  double HBSVBFC5;
  double  position_BSVBFC5[3];
  double  rotation_BSVBFC5[3];

  //Parts No 7(??)  Titan sheild    (G4Polycone)
  int NumZPlanesTSFC5;
  double IrTSFC5[12];
  double OrTSFC5[12];
  double ZTSFC5[12];
  double position_TSFC5[3];
  double rotation_TSFC5[3];

  //Parts No 7(Vac) Titan sheild Vaccum (G4Cons)
  double brMaxTSVFC5;
  double brMinTSVFC5;
  double trMaxTSVFC5;
  double trMinTSVFC5;
  double HTSVFC5;
  double position_TSVFC5[3];
  double rotation_TSVFC5[3];

  //Parts No 06 Bottom Supporter    (G4Polycone)
  int NumZPlanesBotSFC5;
  double IrBotSFC5[4];
  double OrBotSFC5[4];
  double ZBotSFC5[4];

  double position_BotSFC5[3];
  double rotation_BotSFC5[3];

  //Parts No 15  Top Supporter    (G4Polycone)
  int NumZPlanesTopSFC5;
  double IrTopSFC5[6];
  double OrTopSFC5[6];
  double ZTopSFC5[6];

  double position_TopSFC5[3];
  double rotation_TopSFC5[3];

  // Vaccum for Berrel Supporter Bottom  (G4Polycone);
  double rMaxTopSVFC5;
  double rMinTopSVFC5;
  double HTopSVFC5;
  double  position_TopSVFC5[3];
  double  rotation_TopSVFC5[3];

  //Parts No 01&02  CF90 Frange
  int NumZPlanesCF90FC5;
  double IrCF90FC5[6];
  double OrCF90FC5[6];
  double ZCF90FC5[6];

  double position_CF90FC5[5];
  double rotation_CF90FC5[5];

  //Parts No 06b   TopSupporter
  int NumZPlanesbTopSFC5;
  double IrbTopSFC5[4];
  double OrbTopSFC5[4];
  double ZbTopSFC5[4];

  double position_bTopSFC5[5];
  double rotation_bTopSFC5[5];

  //Parts No 01 Copper chamber
  int NumZPlanesCPFC5;
  double IrCPFC5[10];
  double OrCPFC5[10];
  double ZCPFC5[10];

  double position_CPFC5[3];
  double rotation_CPFC5[3];

  //Parts No 12 Top Plate TP
  double XTPFC5;
  double YTPFC5;
  double ZTPFC5;

  double position_TPFC5[3];
  double rotation_TPFC5[3];

  // FeedThrough No11  (G4 Tubs);
  double rMaxFTFC5;
  double rMinFTFC5;
  double HFTFC5;
  double  position_FTFC5[3];
  double  rotation_FTFC5[3];

  // FeedThrough No11  (G4 Tubs);
  double rMaxFT2FC5;
  double rMinFT2FC5;
  double HFT2FC5;
  double  position_FT2FC5[3];
  double  rotation_FT2FC5[3];

  // Vaccum Duct  Vacuum (G4 Tubs);
  double rMaxVDVFC5;
  double rMinVDVFC5;
  double HVDVFC5;
  double  position_VDVFC5[3];
  double  rotation_VDVFC5[3];

  // Vaccum Duct  (G4 Tubs);
  double rMaxVDFC5;
  double rMinVDFC5;
  double HVDFC5;
  double  position_VDFC5[3];
  double  rotation_VDFC5[3];



  // Screen Monitor 3
  double size_of_ScreenMonitor3[3];
  double position_ScreenMonitor3[3];
  double rotation_ScreenMonitor3[3];

  // Beam Attenuator 
  double RMin_BeamAttenuator;
  double RMax_BeamAttenuator;
  double L_BeamAttenuator;
  double position_BeamAttenuator[3];
  double rotation_BeamAttenuator[3];

  // Pb Collimator added in 2012.09.24
  double RMin_PbCollimator;
  double RMax_PbCollimator;
  double L_PbCollimator;
  double position_PbCollimator[3];
  double rotation_PbCollimator[3];

  // Virtual Test Chamber                            
  // Cylinder( Material ) parts 1-5    
  double RMin_VirtualChamberCylinder[5];
  double RMax_VirtualChamberCylinder[5];
  double L_VirtualChamberCylinder[5];
  double position_VirtualChamberCylinder[5][3];
  double rotation_VirtualChamberCylinder[5][3];

  // Target : added in 2013.04.17, to check backscattering effect
  double RMin_VirtualChamberTarget;
  double RMax_VirtualChamberTarget;
  double L_VirtualChamberTarget;
  double position_VirtualChamberTarget[3];
  double rotation_VirtualChamberTarget[3];

  // Body : Inner side ( No conductor ) parts 1-4 
  double RMin_VirtualChamberBody[4];
  double RMax_VirtualChamberBody[4];
  double L_VirtualChamberBody[4];
  double position_VirtualChamberBody[4][3];
  double rotation_VirtualChamberBody[4][3];


  //--------------------------------------
  // ICE information
  //--------------------------------------
  double ice_position[3];
  double ice_size[3];
  double ice_rotation[3];


  //--------------------------------------
  // G4Tube List
  //--------------------------------------
  int NumBLComp;
  double RMinBL_FF[100];
  double RMaxBL_FF[100];
  double L_BL_FF[100];
  double RMinBL[100];
  double RMaxBL[100];
  double L_BL[100];
  double RMinBL_BF[100];
  double RMaxBL_BF[100];
  double L_BL_BF[100];
  double position_BL_FF[100][3];
  double position_BL[100][3];
  double position_BL_BF[100][3];
  double rotation_BL[100][3];
  void inputBeamLineG4TubeList(void);


  // Magnetic Field
  double quadrupole_magnet1_field_region[3][2];
  double quadrupole_magnet2_field_region[3][2];

  double bending_magnet_field_region[3][2];

  //----------------------------------------
  // Beam injection Position
  // Z=0;
  //----------------------------------------
  double beam_injection_position[3];

  //---------------------------------------
  // Pb Block 
  //---------------------------------------

  // PB Block tower 
  //  0: Forward
  //  1: Backward
  //  2: Side
  //  3: Side
  //  4: Straight Line 

  double standard_pb_L;
  double standard_pb_W;
  double standard_pb_H;

  double additional_pb_hight;

  double size_of_PbTower[5][3];
  double position_of_PbTower[5][3];
  
  // BM Duct Pb Block
  //  0: Forward
  //  1: Backward
  // parts-1
  double Rmin_PbBlock_BM_parts1[2];
  double Rmax_PbBlock_BM_parts1[2];
  double L_PbBlock_BM_parts1[2];
  double Phi_start_PbBlock_BM_parts1[2];
  double Phi_delta_PbBlock_BM_parts1[2];
  double position_of_PbBlock_BM_parts1[2][3];
  double rotation_of_PbBlock_BM_parts1[2][3];
  // parts-2
  double size_of_PbBlock_BM_parts2[2][3];
  double position_of_PbBlock_BM_parts2[2][3]; 
  // parts-3
  double size_of_PbBlock_BM_parts3[2][3];
  double position_of_PbBlock_BM_parts3[2][3];

  //---------------------------------------
  // Concrete Block
  //---------------------------------------
  int numConcrete;

  double thick_concrete;
  double height_concrete;
  double posiz_concrete;
  double size_of_ConcreteBlock[50][3];
  double position_of_ConcreteBlock[50][3];


  //---------------------------------------
  // Detector Planes
  //---------------------------------------
  double thin_dp;
  int numDetectorPlane;
  double size_of_DetectorPlane[200][3];
  double position_of_DetectorPlane[200][3];

  

};
#endif
