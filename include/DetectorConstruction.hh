//===================================================================
//    DetectorConstruction.hh 
//      Author    : T.Shibata
//      Creation  : 2009.11.30
//      Last Update ; 2010.11.01
//      Last Update ; 2011.02.03
//      Last Update ; 2011.05.01
//      Last Update : 2011.12.12               
//      Last Update : 2011.12.13               
//      Last Update : 2011.12.14 
//      Last Update : 2011.12.25
//      Last Update : 2012.01.07
//      Last Update : 2012.05.01
//      Last Update : 2012.09.25
//      Last Update : 2012.10.02
//===================================================================
#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4SystemOfUnits.hh"

#include "ReadFile.hh"

#include "GeoMagneticField.hh"
#include "ELSMagneticField.hh"
#include "ELSParameters.hh"
#include "AtmosfericFormula.hh"

#include "RayleighFormula.hh"
#include "MieFormula.hh"

class G4Box;
class G4Tubs;
class G4Cons;
class G4Polycone;
class G4Polyhedra;
class G4Sphere;
class G4Trd;
class G4UnionSolid;
class G4SubtractionSolid;
class G4IntersectionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4AssemblyVolume;

class ELSParameters;

class ELSMagnetic;
class GeoMagnetic;

class AtmForm;

class Rayleigh;
class Mie;

// User Detector Construction 
class DetectorConstruction : public G4VUserDetectorConstruction
{                                                                                     
  public:
 
  DetectorConstruction();
  ~DetectorConstruction();

  public: 
  G4VPhysicalVolume* Construct();

  void setup_detector_parameter(Plist *plist);

  void setTargetMaterial (G4String);

  ELSMagnetic *elsmagnetic;

  ELSParameters *elsparam;
  
  AtmForm * atmform;

  GeoMagnetic *geomagneticfield;

  Rayleigh *rayleigh;
  Mie *mie;

  private:

  ofstream OutputFileAir;

  // Air Condition Parameters
  G4double Air_composition_nitrogen[4];
  G4double Air_composition_oxygen[4];
  G4double Air_composition_argon[4];
  G4double Air_composition_carbon_dioxide[4];
 
  G4double Air_temperature[4];
  G4double Air_pressure[4];
  G4double Air_relative_humidity[4];

  G4double saturated_vapor_pressure;
  G4double vapor_pressure;
  G4double vapor_pressure_rate;

  G4double air_mass_density;

  G4double air_composition_nitrogen;
  G4double air_composition_oxygen;
  G4double air_composition_argon;
  G4double air_composition_carbon_dioxide;
  G4double air_composition_water;

  // Air Refraction index and Attenuation length, added in 2012.09.24
  double lambda_min;
  double lambda_max;
  vector<double> Mlambda;
  vector<double> RefractionIndexAir;
  vector<double> RaylieghAttenuation;
  vector<double> MieAttenuation;
  double Mie_forward_g;
  double Mie_barckward_g;
  double Mie_r;

  // GeoMagneticField Parameters 
  G4int    geomagnetic_onoff;
  G4double geomagnetic_theta[4];
  G4double geomagnetic_declination[4];
  G4double geomagnetic_horizontal[4];
  G4double geomagnetic_vertical[4];

  // Plist Parameters
  G4double CuCollimeter_Diameter;
  G4double Collimeter_Width;
  
  G4int Faradaycup_flag;

  G4int Faradaycup4_flag;

  G4int ScreenMonitor3_flag;

  G4int Cerenkov_flag;
  
  G4int BeamAttenuator_flag;
  G4int BeamAttenuator_material;

  G4int PbCollimator_flag;

  G4int Virtual_chamber_flag;
  G4double vc_temperature;
  G4double vc_pressure;

  G4double vc_mass_density;

  G4int QM_Magnetic_field_flag;
   G4double QM1_Magnetic_field;
   G4double QM2_Magnetic_field;

  G4int BM_Magnetic_field_flag; 
   G4double BM_Magnetic_field;

  // World = Experimental Hole
  G4Box*             solidWorld;   
  G4LogicalVolume*   logicWorld;   
  G4VPhysicalVolume* physiWorld;   

  G4double SizeOfWorld[3];   

  // Ground Plane ( for Kill all of Particles )
  G4Box*             solidGroundPlane;   
  G4LogicalVolume*   logicGroundPlane;   
  G4VPhysicalVolume* physiGroundPlane;

  G4double GroundPlane_Thickness;
  G4double SizeOfGroundPlane[3];
  G4double PositionOfGroundPlane[3];
  G4double RotationOfGroundPlane[3];

  // ConcretePad
  G4Box*             solidConcretepad;
  G4LogicalVolume*   logicConcretepad;
  G4VPhysicalVolume* physiConcretepad;

  // ELS Container
  G4SubtractionSolid* solidELSContainer;
  G4Box*              solidELSContainer_outer;   
  G4Box*              solidELSContainer_inner;    

  G4Box*              solidELSContainer_hole;
  G4SubtractionSolid* solidELSContainer_withhole;
   
  G4LogicalVolume*   logicELSContainer;   
  G4VPhysicalVolume* physiELSContainer;   

  // ELS Cover Box
  G4Box*             solidCoverBox_outer; 
  G4Box*             solidCoverBox_inner; 
  G4SubtractionSolid* solidsolidCoverBox;

  G4LogicalVolume*   logicCoverBox;
  G4VPhysicalVolume* physiCoverBox;

  G4Box*             solidCoverBox_board;
  G4Tubs*            solidCoverBox_boardhole;
  G4SubtractionSolid* solidCoverBox_board_withhole;

  G4LogicalVolume*   logicCoverBoxBoard;
  G4VPhysicalVolume* physiCoverBoxBoard;

  // Base A & B & C in ELS Container
  G4Box*             solidBase[3];   
  G4LogicalVolume*   logicBase[3];
  G4VPhysicalVolume* physiBase[3];

  // Beam Line 
  // G4Tubs 
  G4Tubs*            solidBL_tubs[30];
  G4LogicalVolume*   logicBL_tubs[30];   
  G4VPhysicalVolume* physiBL_tubs[30];

  G4Tubs*            solidBL_vac_tubs[30];
  G4LogicalVolume*   logicBL_vac_tubs[30];   
  G4VPhysicalVolume* physiBL_vac_tubs[30];

  G4Tubs*            solidBL_FF_tubs[30];
  G4LogicalVolume*   logicBL_FF_tubs[30];
  G4VPhysicalVolume* physiBL_FF_tubs[30];

  G4Tubs*            solidBL_FF_vac_tubs[30];
  G4LogicalVolume*   logicBL_FF_vac_tubs[30];
  G4VPhysicalVolume* physiBL_FF_vac_tubs[30];

  G4Tubs*            solidBL_BF_tubs[30];
  G4LogicalVolume*   logicBL_BF_tubs[30];
  G4VPhysicalVolume* physiBL_BF_tubs[30];

  G4Tubs*            solidBL_BF_vac_tubs[30];
  G4LogicalVolume*   logicBL_BF_vac_tubs[30];
  G4VPhysicalVolume* physiBL_BF_vac_tubs[30];

  // EGUNBlank Frange
  G4Tubs*            solidEGUNBlankFrange;
  G4LogicalVolume*   logicEGUNBlankFrange;
  G4VPhysicalVolume* physiEGUNBlankFrange;

  // Anode Frange
  G4Tubs*            solidAnode;
  G4LogicalVolume*   logicAnode;
  G4VPhysicalVolume* physiAnode;

  // Screen Monitor
  G4Tubs*            solidSM_FF_tubs[2];
  G4LogicalVolume*   logicSM_FF_tubs[2];
  G4VPhysicalVolume* physiSM_FF_tubs[2];

  G4Tubs*            solidSM_FF_vac_tubs[2];
  G4LogicalVolume*   logicSM_FF_vac_tubs[2];
  G4VPhysicalVolume* physiSM_FF_vac_tubs[2];

  G4Tubs*            solidSM_BF_tubs[2];
  G4LogicalVolume*   logicSM_BF_tubs[2];
  G4VPhysicalVolume* physiSM_BF_tubs[2];

  G4Tubs*            solidSM_BF_vac_tubs[2];
  G4LogicalVolume*   logicSM_BF_vac_tubs[2];
  G4VPhysicalVolume* physiSM_BF_vac_tubs[2];

  G4Tubs*            solidSM_parts1_tubs[2];
  G4Tubs*            solidSM_parts12_tubs[2];
  G4Tubs*            solidSM_parts13_tubs[2]; // use for vacuum
  G4Tubs*            solidSM_parts2_tubs[2];
  G4Tubs*            solidSM_parts3_tubs[2];

  G4SubtractionSolid* solidSM_duct[2];
  G4LogicalVolume*    logicSM_duct[2];
  G4VPhysicalVolume*  physiSM_duct[2];

  G4SubtractionSolid* solidSM_duct_vac[2];
  G4LogicalVolume*    logicSM_duct_vac[2];
  G4VPhysicalVolume*  physiSM_duct_vac[2];

  G4SubtractionSolid* solidSM_body[2];
  G4LogicalVolume*    logicSM_body[2];
  G4VPhysicalVolume*  physiSM_body[2];

  G4Tubs*            solidSM_body_vac[2];
  G4LogicalVolume*   logicSM_body_vac[2];
  G4VPhysicalVolume* physiSM_body_vac[2];

  // BM Duct 
  // Cu Collimter
  G4Box*              solidCuCollimeter_body;
  G4Tubs*             solidCuCollimeter_hole;
  G4Tubs*             solidCuCollimeter_vacuum;

  G4SubtractionSolid* solidCuCollimeter;
  G4LogicalVolume*    logicCuCollimeter;
  G4LogicalVolume*    logicCuCollimeter_vacuum;
  G4VPhysicalVolume*  physiCuCollimeter;
  G4VPhysicalVolume*  physiCuCollimeter_vacuum;

  //  BM Duct Frange
  G4Tubs*             solidBMF_outer;
  G4Box*              solidBMF_inner;
  G4SubtractionSolid* solidBMF;
  G4LogicalVolume*    logicBMF;
  G4VPhysicalVolume*  physiBMF[3];

  G4Box*              solidBMF_vac;
  G4LogicalVolume*    logicBMF_vac;
  G4VPhysicalVolume*  physiBMF_vac[3];
  // BM Duct parts-1-3
  G4Box*             solidBMLduct_outer[3];
  G4Box*             solidBMLduct_inner[3];
  G4SubtractionSolid* solidBMLduct[3];
  G4LogicalVolume*   logicBMLduct[3];
  G4VPhysicalVolume* physiBMLduct[3];

  G4Box*             solidBMLduct_vac[3];
  G4LogicalVolume*   logicBMLduct_vac[3];
  G4VPhysicalVolume* physiBMLduct_vac[3];
  // BM Duct parts-4(1-9)
  // parts 4-(1)-(5)
  G4Box*              solidBMmainBody[2];
  G4Tubs*             solidBMmainSub1;
  G4Box*              solidBMmainSub2;
  G4Tubs*             solidBMmainSub3;
  G4SubtractionSolid* solidBMmain1[2];  
  G4SubtractionSolid* solidBMmain2[2];
  G4SubtractionSolid* solidBMmain[2];
  G4LogicalVolume*    logicBMmain[2];
  G4VPhysicalVolume*  physiBMmain[2];
  // parts 4-(6)
  G4Tubs*             solidBMSub1;
  G4LogicalVolume*    logicBMSub1;
  G4VPhysicalVolume*  physiBMSub1;
  // parts 4-(7)
  G4Box*              solidBMSub2;
  G4LogicalVolume*    logicBMSub2;
  G4VPhysicalVolume*  physiBMSub2;
  // parts 4-(8)
  G4Tubs*             solidBMSub3;
  G4LogicalVolume*    logicBMSub3;
  G4VPhysicalVolume*  physiBMSub3;
  // parts 4-(9)
  G4Box*              solidBMSub4;
  G4LogicalVolume*    logicBMSub4;
  G4VPhysicalVolume*  physiBMSub4;

  // BM Duct vacuum ( parts5-(1)-(4)
  G4Box*              solidBMmainBody_vac;
  G4Tubs*             solidBMmainSub1_vac;
  G4Box*              solidBMmainSub2_vac;
  G4Tubs*             solidBMmainSub3_vac;  
  G4SubtractionSolid* solidBMmain1_vac;
  G4SubtractionSolid* solidBMmain2_vac;
  G4SubtractionSolid* solidBMmain_vac;
  G4LogicalVolume*    logicBMmain_vac;
  G4VPhysicalVolume*  physiBMmain_vac;

  // BM York & Coil
  G4Box*              solidBMYork[5];
  G4LogicalVolume*    logicBMYork[5];
  G4VPhysicalVolume*  physiBMYork[5];
 
  G4Box*              solidBMCoilOuter;
  G4Box*              solidBMCoilInner;
  G4SubtractionSolid* solidBMCoil;
  G4LogicalVolume*    logicBMCoil;
  G4VPhysicalVolume*  physiBMCoil[2];

  // SLIT
  G4Tubs*            solidSLIT_FF_tubs;
  G4LogicalVolume*   logicSLIT_FF_tubs;
  G4VPhysicalVolume* physiSLIT_FF_tubs;

  G4Tubs*            solidSLIT_FF_vac_tubs;
  G4LogicalVolume*   logicSLIT_FF_vac_tubs;
  G4VPhysicalVolume* physiSLIT_FF_vac_tubs;

  G4Tubs*            solidSLIT_BF_tubs;
  G4LogicalVolume*   logicSLIT_BF_tubs;
  G4VPhysicalVolume* physiSLIT_BF_tubs;

  G4Tubs*            solidSLIT_BF_vac_tubs;
  G4LogicalVolume*   logicSLIT_BF_vac_tubs;
  G4VPhysicalVolume* physiSLIT_BF_vac_tubs;

  G4Tubs*            solidSLIT_parts1_tubs;
  G4Tubs*            solidSLIT_parts12_tubs;
  G4Tubs*            solidSLIT_parts13_tubs; // use for vacuum
  G4Tubs*            solidSLIT_parts2_tubs;
  G4Tubs*            solidSLIT_parts3_tubs;

  G4SubtractionSolid *solidSLIT_duct;
  G4LogicalVolume*   logicSLIT_duct;
  G4VPhysicalVolume* physiSLIT_duct;

  G4SubtractionSolid* solidSLIT_duct_vac;
  G4LogicalVolume*    logicSLIT_duct_vac;
  G4VPhysicalVolume*  physiSLIT_duct_vac;

  G4SubtractionSolid* solidSLIT_body;
  G4LogicalVolume*   logicSLIT_body;
  G4VPhysicalVolume* physiSLIT_body;

  G4Box*              solidCollimeter_vac;
  G4Tubs*             solidSLIT_body_vac_0;
  G4SubtractionSolid* solidSLIT_body_vac_1;
  
  G4SubtractionSolid* solidSLIT_body_vac;
  G4LogicalVolume*    logicSLIT_body_vac;
  G4VPhysicalVolume*  physiSLIT_body_vac;

  G4Box*             solidCollimeter[2];
  G4LogicalVolume*   logicCollimeter[2];
  G4VPhysicalVolume* physiCollimeter[2];

  // Blank Frange
  G4Tubs*            solidBlankFrange;
  G4LogicalVolume*   logicBlankFrange;
  G4VPhysicalVolume* physiBlankFrange;

  // Output Frange & Window
  G4Tubs*            solidWindow;
  G4LogicalVolume*   logicWindow;
  G4VPhysicalVolume* physiWindow;

  G4Tubs*            solidWindow_vac;
  G4LogicalVolume*   logicWindow_vac;
  G4VPhysicalVolume* physiWindow_vac;

  G4Tubs*            solidTiWindow;
  G4LogicalVolume*   logicTiWindow;
  G4VPhysicalVolume* physiTiWindow;

  // Top Plate
  G4Box*              solidTopPlate_body;
  G4Tubs*             solidTopPlate_hole;
  G4SubtractionSolid* solidTopPlate;
  G4LogicalVolume*    logicTopPlate;
  G4VPhysicalVolume*  physiTopPlate;

  // Faraday Cup
  // Al Cylinder
  G4Tubs*            solidFDAlCylinder[3];
  G4LogicalVolume*   logicFDAlCylinder[3];
  G4VPhysicalVolume* physiFDAlCylinder[3];
  // FD Body
  G4Tubs*            solidFDBody[4];
  G4LogicalVolume*   logicFDBody[4];
  G4VPhysicalVolume* physiFDBody[4];

  //--------------------------------------------------------------------------------------------
  // FC2-> FC4 in 2012.05.01                                                                                    
  // Front : Carbon                                                                                             
  //       --> Front1/2 thin shield ( 1 and 2 ) Cu or Al or Any other material ?  in 2012.05.01
  // Creat Faraday Cup 4 in 2012.10.02                                                                                             
  // Geometry was fixed except position 

  // Ground Layer (G4Polycone)
  G4Polycone*        solidGNDLayerFC4;
  G4LogicalVolume*   logicGNDLayerFC4;
  G4VPhysicalVolume* physiGNDLayerFC4;
  int numZPlanesGNDLayerFC4;
  G4double irGNDLayerFC4[8];
  G4double orGNDLayerFC4[8];
  G4double zGNDLayerFC4[8];

  // Shield Layer (G4Polycone)
  G4Polycone*        solidShieldLayerFC4;
  G4LogicalVolume*   logicShieldLayerFC4;
  G4VPhysicalVolume* physiShieldLayerFC4;
  int numZPlanesShieldLayerFC4;
  G4double irShieldLayerFC4[14];
  G4double orShieldLayerFC4[14];
  G4double zShieldLayerFC4[14];

  // Body : Copper (G4Polycone) 
  G4Polycone*        solidBodyFC4;
  G4LogicalVolume*   logicBodyFC4;
  G4VPhysicalVolume* physiBodyFC4;
  int numZPlanesBodyFC4;
  G4double irBodyFC4[6];
  G4double orBodyFC4[6];
  G4double zBodyFC4[6];

  // Isolator-1 : Macor 
  G4Tubs*            solidIso1FC4;
  G4LogicalVolume*   logicIso1FC4;
  G4VPhysicalVolume* physiIso1FC4;

  // Isolator-2 : Macor
  G4Tubs*            solidIso2FC4;
  G4LogicalVolume*   logicIso2FC4;
  G4VPhysicalVolume* physiIso2FC4;

  // Top Plate : A5052 
  G4Box*             solidTopPlateFC4;
  G4LogicalVolume*   logicTopPlateFC4;
  G4VPhysicalVolume* physiTopPlateFC4;



  //======================================
  // Faraday Cup 5 -- made by B.K.Shin-san  
  // added in 2013.04.23 

  // Body(BeamDump) : Copper (G4Polycone)                                                                                                         
  G4Polycone*        solidBodyFC5;
  G4LogicalVolume*   logicBodyFC5;
  G4VPhysicalVolume* physiBodyFC5;
  int numZPlanesBodyFC5;
  G4double irBodyFC5[4];
  G4double orBodyFC5[4];
  G4double zBodyFC5[4];
 
  // BS : No.04 Barrel Supporter (G4Polycone)                                                                                           
  G4Polycone*        solidBSFC5;
  G4LogicalVolume*   logicBSFC5;
  G4VPhysicalVolume* physiBSFC5;
  int numZPlanesBSFC5;
  G4double irBSFC5[8];
  G4double orBSFC5[8];
  G4double zBSFC5[8];
  
  // BSV : No.04 Vac Barrel Supporter Vac (G4Tubs)                                                                                      
  G4Tubs*        solidBSVFC5;
  G4LogicalVolume*   logicBSVFC5;
  G4VPhysicalVolume* physiBSVFC5;


  // BSVB : No.04 Vac Barrel Supporter bottom Vac (G4Tubs)                                                                              
  G4Tubs*        solidBSVBFC5;
  G4LogicalVolume*   logicBSVBFC5;
  G4VPhysicalVolume* physiBSVBFC5;


  // TS : No.07(??) Titan shield  (G4Polycone)                                                                                          
  G4Polycone*        solidTSFC5;
  G4LogicalVolume*   logicTSFC5;
  G4VPhysicalVolume* physiTSFC5;
  int numZPlanesTSFC5;
  G4double irTSFC5[12];
  G4double orTSFC5[12];
  G4double zTSFC5[12];
  
  // TSV : No.07 Titan shieid Vac (G4Cons)  
  G4Cons*            solidTSVFC5;
  G4LogicalVolume*   logicTSVFC5;
  G4VPhysicalVolume* physiTSVFC5;

  // BotS : No.06 Bottom Supporter
  G4Polycone*        solidBotSFC5;
  G4LogicalVolume*   logicBotSFC5;
  G4VPhysicalVolume* physiBotSFC5;
  int numZPlanesBotSFC5;
  G4double irBotSFC5[4];
  G4double orBotSFC5[4];
  G4double zBotSFC5[4];

  // TopS : No.15 Bottom Supporter  
  G4Polycone*        solidTopSFC5;
  G4LogicalVolume*   logicTopSFC5;
  G4VPhysicalVolume* physiTopSFC5;
  int numZPlanesTopSFC5;
  G4double irTopSFC5[6];
  G4double orTopSFC5[6];
  G4double zTopSFC5[6];

  // TopSV : No.15 top Supporter bottom Vac (G4Tubs)
  G4Tubs*        solidTopSVFC5;
  G4LogicalVolume*   logicTopSVFC5;
  G4VPhysicalVolume* physiTopSVFC5;

  // CF90: No.01 & 02 CF90 
  G4Polycone*        solidCF90FC5;
  G4LogicalVolume*   logicCF90FC5;
  G4VPhysicalVolume* physiCF90FC5;
  int numZPlanesCF90FC5;
  G4double irCF90FC5[6];
  G4double orCF90FC5[6];
  G4double zCF90FC5[6];

  // Top supporter bTopS : No 06b  
  G4Polycone*        solidbTopSFC5;
  G4LogicalVolume*   logicbTopSFC5;
  G4VPhysicalVolume* physibTopSFC5;
  int numZPlanesbTopSFC5;
  G4double irbTopSFC5[4];
  G4double orbTopSFC5[4];
  G4double zbTopSFC5[4];

  // Copper Chamber  CP  : No 01     
  G4Polycone*        solidCPFC5;
  G4LogicalVolume*   logicCPFC5;
  G4VPhysicalVolume* physiCPFC5;
  int numZPlanesCPFC5;
  G4double irCPFC5[10];
  G4double orCPFC5[10];
  G4double zCPFC5[10];

  // Top  Plate   : No 12
  G4Box*        solidTPFC5;
  G4LogicalVolume*   logicTPFC5;
  G4VPhysicalVolume* physiTPFC5;

  // FT FeedThrough No 11  
  G4Tubs*        solidFTFC5;
  G4LogicalVolume*   logicFTFC5;
  G4VPhysicalVolume* physiFTFC5;

  // FT2 FeedThrough No 11 
  G4Tubs*        solidFT2FC5;
  G4LogicalVolume*   logicFT2FC5;
  G4VPhysicalVolume* physiFT2FC5;

  // VD : Vaccum duct No  
  G4Tubs*        solidVDFC5;
  G4LogicalVolume*   logicVDFC5;
  G4VPhysicalVolume* physiVDFC5;

  // VDV : Vaccum duct VacuumNo   
  G4Tubs*        solidVDVFC5;
  G4LogicalVolume*   logicVDVFC5;
  G4VPhysicalVolume* physiVDVFC5;

  //======================================


  // Screen Monitor3 
  G4Box*             solidSM3;
  G4LogicalVolume*   logicSM3;
  G4VPhysicalVolume*   physiSM3;

  // Beam Attenuator
  G4Tubs*            solidBeamAttenuator;
  G4LogicalVolume*   logicBeamAttenuator;
  G4VPhysicalVolume* physiBeamAttenuator;

  // Pb Collimator added in 2012.09.25
  G4Tubs*            solidPbCollimator;
  G4LogicalVolume*   logicPbCollimator;
  G4VPhysicalVolume* physiPbCollimator;

  // Virtual Test Chamber added in 2011.12.13                                 
  //  -->  modified in 2013.04.17, to check backscattering effect  
  //                     all of gemetory parameters are fixed           
  //  
  // Cylinder
  G4Tubs*            solidVirtualChamberCylinder[5];
  G4LogicalVolume*   logicVirtualChamberCylinder[5];
  G4VPhysicalVolume* physiVirtualChamberCylinder[5];
 
  // Target 
  G4Tubs*            solidVirtualChamberTarget;
  G4LogicalVolume*   logicVirtualChamberTarget;
  G4VPhysicalVolume* physiVirtualChamberTarget;

  // Vacuum Region
  G4Tubs*             solidVirtualChamberCylinder_vacbase;
  G4Tubs*             solidVirtualChamberCylinder_vactarget;
  G4SubtractionSolid* solidVirtualChamberCylinder_vac;
  G4LogicalVolume*    logicVirtualChamberCylinder_vac; 
  G4VPhysicalVolume*  physiVirtualChamberCylinder_vac;

  // Body  // -> not be used in 2013.04.17  
  G4Tubs*            solidVirtualChamberBody[4];
  G4LogicalVolume*   logicVirtualChamberBody[4];
  G4VPhysicalVolume* physiVirtualChamberBody[4];

  //---------------------------------
  // Pb Block
  //---------------------------------

  G4Box*             solidPbTower[5];   
  G4LogicalVolume*   logicPbTower[5];   
  G4VPhysicalVolume* physiPbTower[5];
  
  G4Tubs*            solidPbBlockBMparts1[2];
  G4LogicalVolume*   logicPbBlockBMparts1[2];
  G4VPhysicalVolume* physiPbBlockBMparts1[2];

  G4Box*             solidPbBlockBMparts2[2];
  G4LogicalVolume*   logicPbBlockBMparts2[2];
  G4VPhysicalVolume* physiPbBlockBMparts2[2];

  G4Box*             solidPbBlockBMparts3[2];
  G4LogicalVolume*   logicPbBlockBMparts3[2];
  G4VPhysicalVolume* physiPbBlockBMparts3[2];

  // Concrete Block
  G4Box*             solidConcreteBlock[20];   
  G4LogicalVolume*   logicConcreteBlock[20];   
  G4VPhysicalVolume* physiConcreteBlock[20];

  G4double SizeOfConcreteBlock[20][3];
  G4double PositionOfConcreteBlock[20][3];
  G4double RotationOfConcreteBlock[20][3];

  //=========================================
  // Detector Plane
  //=========================================
  G4int Outoputfileflag_detectorplane;

  // Particle Counting for Shower Development 
  G4Box*             solidDetectorPlane[200];   
  G4LogicalVolume*   logicDetectorPlane[200];   
  G4VPhysicalVolume* physiDetectorPlane[200];
  

  //========================================
  // added in 2013.04.22 
  //  source was made by Drmitri-san 
  //========================================
  int Flag_DmitriFC;

  G4Material* fAir;    // The air we breathe
  G4Material* fCopper; // Copper used for making the faraday cup and the copper plates
  G4Material* fTitanium; // Titanium used for making the tianium plates

  // The Faraday Cup
  G4Tubs* fSolidFaradayCup;
  G4LogicalVolume* fLogicFaradayCup;
  G4VPhysicalVolume* fPhysiFaradayCup;

  // The upper titanium plate
  G4Tubs* fSolidUpperTiPlate;
  G4LogicalVolume* fLogicUpperTiPlate;
  G4VPhysicalVolume* fPhysUpperTiPlate;

  // The upper copper plate
  G4Tubs* fSolidUpperCuPlate;
  G4LogicalVolume* fLogicUpperCuPlate;
  G4VPhysicalVolume* fPhysUpperCuPlate;

  // Count air plate 
  G4Tubs* fSolidCountAirPlate;
  G4LogicalVolume* fLogicCountAirPlate;
  G4VPhysicalVolume* fPhysCountAirPlate;

  // The bottom copper plate
  G4Tubs* fSolidBottomCuPlate;
  G4LogicalVolume* fLogicBottomCuPlate;
  G4VPhysicalVolume* fPhysBottomCuPlate;

  // The bottom titanium plate
  G4Tubs* fSolidBottomTiPlate;
  G4LogicalVolume* fLogicBottomTiPlate;
  G4VPhysicalVolume* fPhysBottomTiPlate;


  //=======================================
  // FD Dummy Mirror , added in 2013.09.10 
  //=======================================
  G4Sphere *solidMirrorSphere;   
  G4Tubs   *solidMirrorTubs;
  G4IntersectionSolid *solidMirror;
  G4LogicalVolume *logicMirror;
  G4VPhysicalVolume *physiMirror6;
  G4VPhysicalVolume *physiMirror7;

  //=======================================
  // ICE addition
  //=======================================
  G4Box* solidIce;
  G4LogicalVolume* logicIce;
  G4VPhysicalVolume* physiIce;

};
#endif
