//===================================================================
//      FD Detector Constructor
//      Author    : T.Shibata
//      Creation  : 2009.12.01
//      Last Update : 2010.11.11
//      Last Update : 2011.05.15
//      Last Update : 2011.12.12
//      Last Update : 2011.12.13
//      Last Update : 2011.12.14
//      Last Update : 2011.12.21
//      Last Update : 2011.12.25
//      Last Update : 2012.01.07
//      Last Update : 2012.01.20
//      Last Update : 2012.01.27
//      Last Update : 2012.05.01
//      Last Update : 2012.09.24,25
//      Last Update : 2012.10.02
//===================================================================
#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"      
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"

#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "PhysicsParameters.hh"

#include <CLHEP/Vector/ThreeVector.h>
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "DetectorPlaneName.hh"

#include <iomanip>

#include "stdlib.h"
#include "stdio.h"

DetectorConstruction::DetectorConstruction( )
{
  OutputFileAir.open("geant4-air-condition.dat");

  elsmagnetic = new ELSMagnetic();
  geomagneticfield = new GeoMagnetic();
  elsparam = new ELSParameters();
  atmform = new AtmForm();

  rayleigh = new Rayleigh();
  mie = new Mie();
}

DetectorConstruction::~DetectorConstruction()
{
  OutputFileAir.close();
  
  delete elsmagnetic;
  delete geomagneticfield;
  delete elsparam;  
  delete atmform;

  delete rayleigh;
  delete mie;

}

void DetectorConstruction::setup_detector_parameter(Plist *plist)
{

  cout << " ===== SetUp Detector Parameter ======= " << endl;
  
  for( int i=0; i<3; i++ ){
      // World Geometory
      SizeOfWorld[i] = plist->SizeOfWorld(i);
  }

  for( int i=0; i<2; i++ ){
      // Ground Geometory
      SizeOfGroundPlane[i] = SizeOfWorld[i];
  }
  SizeOfGroundPlane[2]=plist->ELSSiteGroundHeight(0); 

  Outoputfileflag_detectorplane = plist->Outoputfileflag_detectorplane();
  
  GroundPlane_Thickness = 1.0; // unit = mm
  SizeOfGroundPlane[2] = GroundPlane_Thickness;  

  // Air Condition Parameters
  for( int i=0; i<4; i++ ){
    Air_composition_nitrogen[i]       = plist->Air_composition_nitrogen(i);
    Air_composition_oxygen[i]         = plist->Air_composition_oxygen(i);
    Air_composition_argon[i]          = plist->Air_composition_argon(i);
    Air_composition_carbon_dioxide[i] = plist->Air_composition_carbon_dioxide(i);

    Air_temperature[i]       = plist->Air_temperature(i); 
    Air_pressure[i]          = plist->Air_pressure(i);
    Air_relative_humidity[i] = plist->Air_relative_humidity(i);
  }

  // GeoMagneticField Parameters
  geomagnetic_onoff    = plist->Geomagnetic_onoff();
  geomagnetic_theta[0] =  plist->Angle_North_Yaxis_BRM() - plist->Angle_ELS_mirror_6();
  for( int i=0; i<4; i++ ){
       //geomagnetic_theta[i]       =  plist->Geomagnetic_theta(i);
       geomagnetic_declination[i] = plist->Geomagnetic_declination(i);
       geomagnetic_horizontal[i]  = plist->Geomagnetic_horizontal(i);
       geomagnetic_vertical[i]    = plist->Geomagnetic_vertical(i);
  }

  CuCollimeter_Diameter = plist->CuCollimeterDiameter()*mm;
  Collimeter_Width    = plist->SlitWidth()*mm;

  Faradaycup_flag      = plist->Faradaycup_flag();
  Faradaycup4_flag     = plist->Faradaycup4_flag();

  ScreenMonitor3_flag  = plist->Screen_monitor3_flag();

  Cerenkov_flag        = plist->Cerenkov_process_flag();

  BeamAttenuator_flag  = plist->Beam_attenuator_flag();
  BeamAttenuator_material = plist->Beam_attenuator_material();    

  PbCollimator_flag =  plist->Pb_collimator_flag(); // added in 2012.09.25

  Virtual_chamber_flag = plist->Virtual_chamber_flag();
  vc_temperature = plist->VC_temperature();
  vc_pressure    = plist->VC_pressure(); 


  QM_Magnetic_field_flag = plist->QM_Magnetic_field_flag();
    QM1_Magnetic_field = plist->QM1_Magnetic_field();
    QM2_Magnetic_field = plist->QM2_Magnetic_field();
  BM_Magnetic_field_flag = plist->BM_Magnetic_field_flag();  

  BM_Magnetic_field = plist->BM_Magnetic_field();
   
  elsparam->els_geometory_initialize(plist);

  elsmagnetic->elsmagneticfield_initialize(elsparam);
  elsmagnetic->elsmagneticfield_onoff( QM_Magnetic_field_flag,
                                       QM1_Magnetic_field,
				       QM2_Magnetic_field,
                                       BM_Magnetic_field_flag,
                                       BM_Magnetic_field );

  geomagneticfield->Geomagneticfield_onoff(geomagnetic_onoff);
  geomagneticfield->InputGeomagneticInfomation( geomagnetic_theta[0],
						geomagnetic_declination[0],
						geomagnetic_horizontal[0],
						geomagnetic_vertical[0] );

  //========================================
  // added in 2013.04.22
  //  source was made by Drmitri-san
  //========================================
  Flag_DmitriFC=0;
  //

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //-------------------------------------------------------------------------
  // Magnetic field
  //-------------------------------------------------------------------------
  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
    {

     G4FieldManager* fieldELSMgr
       = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldELSMgr->SetDetectorField(elsmagnetic);
     fieldELSMgr->CreateChordFinder(elsmagnetic);

     /*
     G4FieldManager* fieldGeoMgr
         = G4TransportationManager::GetTransportationManager()->GetFieldManager();
     fieldGeoMgr->SetDetectorField(geomagneticfield);
     fieldGeoMgr->CreateChordFinder(geomagneticfield);
     */

   fieldIsInitialized = true;
  }

  //--------- Define the NIST material -----
  // source is made by Dmitri-san
  G4String name;
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1);
  // Air
  fAir = man->FindOrBuildMaterial(name = "G4_AIR");
  // Copper
  fCopper = man->FindOrBuildMaterial(name = "G4_Cu");
  // Titanium
  fTitanium = man->FindOrBuildMaterial(name = "G4_Ti");
  //
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  //

  //--------- Material definition ---------
  G4double a, z;
  G4double density;
  G4int nel;

  G4double air_temperature = ( 273.15 - 30.0 )*kelvin; //( 273.15 + 25.0 )*kelvin;
  G4double air_pressure    = 0.8687*atmosphere; //1.0 * atmosphere;

  // Element
  G4Element* H   = new G4Element("Hydrogen", "H", z=1., a=1.00794*g/mole );   

  G4Element* Be  = new G4Element("Beryrllium","Be", z=4., a=9.012182*g/mole );
  G4Element* B   = new G4Element("Boron", "B",      z=5., a=10.811*g/mole   ); 
  G4Element* C   = new G4Element("Carbon", "C",     z=6., a=12.0107*g/mole  );   

  G4Element* N   = new G4Element("Nitrogen", "N", z=7., a= 14.01*g/mole  );
  G4Element* O   = new G4Element("Oxygen"  , "O", z=8., a= 16.00*g/mole  );
  G4Element* F   = new G4Element("Fluorine", "F", z=9., a= 18.9984032*g/mole  );

  G4Element* Na  = new G4Element("Natrium",    "Na", z=11., a= 22.98976928*g/mole  );  
  G4Element* Mg  = new G4Element("Magnesium",  "Mg", z=12., a= 24.305*g/mole       );
  G4Element* Al  = new G4Element("Alminium",   "Al", z=13., a= 26.98*g/mole        );
  G4Element* Si  = new G4Element("Silicon" ,   "Si", z=14., a= 28.0855*g/mole      );
  G4Element* P   = new G4Element("Phosphorus", "P",  z=15., a= 30.973762*g/mole    );
  G4Element* S   = new G4Element("Sulfur" ,    "S",  z=16., a= 53.092*g/mole       );
  G4Element* Ar  = new G4Element("Argon" ,     "Ar", z=18., a= 39.948*g/mole       );

  G4Element* K   = new G4Element("Potassium", "K",   z=19., a= 39.0983*g/mole  );
  G4Element* Ca  = new G4Element("Calcium",   "Ca",  z=20., a= 40.078*g/mole   );

  G4Element* Ti  = new G4Element("Titanium",  "Ti",  z=22,  a=47.867*g/mole    );

  G4Element* Cr  = new G4Element("Chromium",  "Cr", z=24., a= 51.9961*g/mole   );
  G4Element* Mn  = new G4Element("Manganum" , "Mn", z=25., a= 54.938049*g/mole );
  G4Element* Fe  = new G4Element("Iron" ,     "Fe", z=26., a= 55.845*g/mole    ); 
  G4Element* Ni  = new G4Element("Nickel" ,   "Ni", z=28., a= 58.093*g/mole    );
  G4Element* Cu  = new G4Element("Copper",    "Cu", z=29., a= 63.546*g/mole    );

  G4Element* Ta  = new G4Element("Tantalum",  "Ta", z=73,  a=180.94788*g/mole  );
  G4Element* W   = new G4Element("Tungten",   "W",  z=74,  a=183.84*g/mole     );

  G4Element* Pb  = new G4Element("Lead",      "Pb", z=82., a= 207.19*g/mole    );

  G4Material* CO2 =  new G4Material("CO2", density=1.977*mg/cm3, nel=2 );
  CO2->AddElement(C, 1 );
  CO2->AddElement(O, 2 );

  // Water vapor
  G4Material* H2O =  new G4Material("H2O", density=1.000*g/cm3, nel=2 );
  H2O->AddElement(H, 2 );
  H2O->AddElement(O, 1 );
  
  // ICE Properties
  G4Material* ICE =  new G4Material("ICE", density=0.920*g/cm3, nel=2,
				    kStateSolid, 216.15*kelvin);

  ICE->AddElement(H, 2 );
  ICE->AddElement(O, 1 );
  
  // Give one characteristic property for now.  
  // This is *wrong* however, the chart for ICE properties doesn't
  // extend into radio regime: http://refractiveindex.info/?group=CRYSTALS&material=H2O-ice  
  //const G4int ICE_ENTRIES = 1;
  //G4double ppckov[ICE_ENTRIES]     = {4.0 * eV};
  //G4double rindex[ICE_ENTRIES]     = {1.3015};
  //G4double absorption[ICE_ENTRIES] = {1450.0*cm};
  
  //G4MaterialPropertiesTable *ICE_PROP = new G4MaterialPropertiesTable();
  //ICE_PROP->AddProperty("RINDEX",ppckov,rindex,ICE_ENTRIES);
  //ICE_PROP->AddProperty("ABSLENGTH",ppckov,absorption,ICE_ENTRIES);
  
  //ICE->SetMaterialPropertiesTable(ICE_PROP);
  
  // Aluminium Oxide
  G4Material* Al2O3 = new G4Material("Al2O3",  density=3.97*g/cm3, nel=2 );
  Al2O3->AddElement(Al, 2 );
  Al2O3->AddElement(O, 3  );

  G4Material* Cr2O3 = new G4Material("Cr2O3",  density=5.22*g/cm3, nel=2 );
  Cr2O3->AddElement(Cr, 2 );
  Cr2O3->AddElement(O,  3 );

  // G4Material* CH2 = new G4Material("CH2", density=0.94*g/cm3, nel=2 );
  // G4Material* CH2 = new G4Material("CH2", density=0.953*g/cm3, nel=2 ); // modified in 2012.04.12
  G4Material* C2H4 = new G4Material("C2H4", density=0.953*g/cm3, nel=2 ); // modified in 2012.06.22
  C2H4->AddElement(C, 2 );
  C2H4->AddElement(H, 4 );

  // Silicon Dioxide ( = Silica )
  G4Material* SiO2 = new G4Material("SiO2", density=2.2*g/cm3,  nel=2 );
  SiO2->AddElement(Si, 1);
  SiO2->AddElement(O, 2);
    
  // Magnesium Oxide 
  G4Material* MgO = new G4Material("MgO", density=3.65*g/cm3,  nel=2 );
  MgO->AddElement(Mg, 1);  
  MgO->AddElement(O, 1);

  // Potassium Oxide
  G4Material* K2O = new G4Material("K2O", density=2.35*g/cm3,  nel=2 );
  K2O->AddElement(K, 2);
  K2O->AddElement(O, 1);

  // Boron Trioxide
  G4Material* B2O3 = new G4Material("B2O3", density=1.85*g/cm3,  nel=2 );
  B2O3->AddElement(B, 2);
  B2O3->AddElement(O, 3);

  // Boron Teflon 2200 kg/m3  2.2 *10E+6/1Ecm
  G4Material* C2F4= new G4Material("C2F4_Teflon", density=2.2*g/cm3,  nel=2 );
  C2F4->AddElement(C, 2);
  C2F4->AddElement(F, 4);

  // Air ===============================================================
   // Cacluration of saturated_vapor_pressure 
   // input  = deg-C
   // output = hPa
   //saturated_vapor_pressure = atmform->Saturated_vapor_pressure_weler_hyland( Air_temperature[0] );
  
   OutputFileAir << " " << endl;
   OutputFileAir << " Temperature (deg-C) " << Air_temperature[0]       << endl;
   OutputFileAir << " Pressure    (hPa)   " << Air_pressure[0]          << endl; 
   OutputFileAir << " Relative Humidity   " << Air_relative_humidity[0] << endl;
   OutputFileAir << " " << endl;

   saturated_vapor_pressure = atmform->Saturated_vapor_pressure_tetens( Air_temperature[0] );
   OutputFileAir << " saturated_vapor_pressure (hPa) = " << saturated_vapor_pressure << endl;
  
   // Cacluration of vapor pressure (hPa)
   // input  = % 
   // output = hPa
   vapor_pressure = atmform->Vapor_pressure(Air_relative_humidity[0], saturated_vapor_pressure);
   OutputFileAir << " vapor_pressure (hPa) = " << vapor_pressure << endl;
   vapor_pressure_rate = atmform->Vapor_pressure_rate(Air_pressure[0], vapor_pressure );

   // Air composition ratio
   // output = %
   air_composition_nitrogen = Air_composition_nitrogen[0]*vapor_pressure_rate;
   air_composition_oxygen = Air_composition_oxygen[0]*vapor_pressure_rate;
   air_composition_argon = Air_composition_argon[0]*vapor_pressure_rate;
   air_composition_carbon_dioxide = Air_composition_carbon_dioxide[0]*vapor_pressure_rate;
   air_composition_water = vapor_pressure/Air_pressure[0]*100.0;

   // Air mass density 
   // output = g/cm3
   air_mass_density = atmform->Air_mass_density( Air_temperature[0], Air_pressure[0], vapor_pressure );
   OutputFileAir << " " << endl;
   OutputFileAir << " air_mass_density(g/cm^3) = " << air_mass_density  << endl;
   OutputFileAir << " " << endl;

   G4Material* Air = new G4Material("Air", density= air_mass_density*g/cm3, nel=5, 
	 			           kStateGas, atmform->teme_Kelvin(Air_temperature[0]), Air_pressure[0] );

   Air->AddElement(N, air_composition_nitrogen*perCent);
   Air->AddElement(O, air_composition_oxygen*perCent);

   Air->AddElement(Ar, air_composition_argon*perCent);
   Air->AddMaterial( CO2, air_composition_carbon_dioxide*perCent);
   Air->AddMaterial( H2O, air_composition_water*perCent);

   OutputFileAir << " " << endl;
   OutputFileAir << " N2 (%)  : " << air_composition_nitrogen       << endl;
   OutputFileAir << " O2 (%)  : " << air_composition_oxygen         << endl;
   OutputFileAir << " Ar (%)  : " << air_composition_argon          << endl;
   OutputFileAir << " CO2 (%) : " << air_composition_carbon_dioxide << endl;
   OutputFileAir << " H2O (%) : " << air_composition_water          << endl;
   OutputFileAir << " " << endl;

   OutputFileAir << "Air Condition, Output of geant4 " 
                 << setprecision(10) << " Density = " << Air->GetDensity()         << " "
                                 	   	              << Air->GetDensity()/(g/cm3) << " " 
                 << " Temperature " << Air->GetTemperature() << " (K) " 
	 	                            << atmform->teme_DegC( Air->GetTemperature() )   << " (degC) " 
                 << " Pressure =  " << Air->GetPressure()                            << " (hPa) " 
                 << endl;
   OutputFileAir << " " << endl;

   OutputFileAir <<  "Number of Atom in unit Vomule " << *Air->GetVecNbOfAtomsPerVolume() << "(/mm^3)" <<  " "
                                                      << *Air->GetVecNbOfAtomsPerVolume()*1000.0 << "(/cm^3)"<< endl;  

   // Comment out by T.Shibata in 2012.04.26
   /*
   G4double AirMassConst = 0.01*(*Air->GetVecNbOfAtomsPerVolume())/Avogadro/(g/mole)*1.E3; 
   
   G4double MassFraction_N = air_composition_nitrogen*N->GetA()*AirMassConst; 
   G4double MassFraction_O = air_composition_oxygen*O->GetA()*AirMassConst
	               	         + air_composition_carbon_dioxide*2.0*O->GetA()*AirMassConst
                   		     + air_composition_water*O->GetA()*AirMassConst;
   G4double MassFraction_Ar = air_composition_argon*Ar->GetA()*AirMassConst;
   G4double MassFraction_C  = air_composition_carbon_dioxide*C->GetA()*AirMassConst;
   G4double MassFraction_H  = air_composition_water*2.0*H->GetA()*AirMassConst;  // bugfixed in 2012.04.26 

   //OutputFileAir  <<  "mass density of Air is  "  << 100.0*AirMassConst*28.966*(g/mole)*(Air->GetDensity()/(g/cm3))/0.001293 << endl; 
   OutputFileAir  <<  "mass of Air " << AirMassConst << endl;
   OutputFileAir  <<  "mass density of N "  << MassFraction_N  << endl;  
   OutputFileAir  <<  "mass density of O "  << MassFraction_O  << endl;
   OutputFileAir  <<  "mass density of Ar " << MassFraction_Ar << endl;
   OutputFileAir  <<  "mass density of C "  << MassFraction_C  << endl;
   OutputFileAir  <<  "mass density of H "  << MassFraction_H  << endl; 
   OutputFileAir <<  " Total = " << MassFraction_N + MassFraction_O + MassFraction_Ar + MassFraction_C + MassFraction_H << endl; 
   OutputFileAir << " " << endl;
   */

  // N2 Air 
   vc_mass_density = atmform->N2Air_mass_density( vc_temperature, vc_pressure );
   G4Material* N2Air = new G4Material("N2Air", density= vc_mass_density*g/cm3, nel=1,
				      kStateGas, atmform->teme_Kelvin(vc_temperature), vc_pressure );
   N2Air->AddElement(N, 100.0*perCent);

  //=============================================================================================-
  // --- Air Properties ( Rayleigh and refraction index ) --- added in 2012.09.24
   const G4int NumAir=500;  // Not more than 10000!!!
   const int intNumAir(NumAir);

   //----------------------------                                            
   // SetUp Rayleigh Scattering
   //---------------------------- 
   /*
   lambda_min = 200.0; // unit = m 
   lambda_max = 700.0; // unit = m        
   */
   lambda_min = 250.0; // unit = m 
   lambda_max = 450.0; // unit = m                                          

   rayleigh->attenuation_length( intNumAir,
				                 lambda_min, lambda_max,
 		   	                     Air->GetTemperature(), 
                                 Air->GetPressure()-vapor_pressure,  
                                 vapor_pressure, 
                                 Air->GetDensity()/(g/cm3),
				                 Mlambda, 
                                 RefractionIndexAir, 
                                 RaylieghAttenuation );

   mie->attenuation_length( intNumAir,
							lambda_min, lambda_max,
							Air->GetTemperature(),
							Air->GetPressure(),
							vapor_pressure,
							Air->GetDensity()/(g/cm3),
							Mlambda,
                            MieAttenuation ); 

   Mie_forward_g   = mie->Forward_g();  
   Mie_barckward_g = mie->Backward_g();
   Mie_r           = mie->ForwardBackward_r();

   G4double EphotonAir[NumAir];
   G4double RefractiveAir[NumAir];
   G4double RayleighSactteringLength[NumAir];
   G4double MieSactteringLength[NumAir];

   int iNumAir(Mlambda.size());
   for( vector<double>::iterator it = Mlambda.begin();
	it!=Mlambda.end(); it++ ){
     double it_lambda = *it;
     EphotonAir[iNumAir-1]=PhotonLambdatoEnergy(it_lambda)*eV;
     iNumAir--;
   }

   // Air Refraction Index 
   iNumAir=RefractionIndexAir.size();
   for( vector<double>::iterator it = RefractionIndexAir.begin();
	it!=RefractionIndexAir.end(); it++ ){
     double it_refraction = *it;
     RefractiveAir[iNumAir-1]=it_refraction;
     iNumAir--;
   }

   // Rayliegh Scattering Attenuation length (m)
   iNumAir=RaylieghAttenuation.size();
   for( vector<double>::iterator it = RaylieghAttenuation.begin();
   	    it!=RaylieghAttenuation.end(); it++ ){
        double it_attenuation = *it;
        RayleighSactteringLength[iNumAir-1]=it_attenuation*m;  // attenuation length ( m )
        iNumAir--;
   }
   
   // Mie Scattering Attenuation length (m)
   iNumAir=MieAttenuation.size();
   for( vector<double>::iterator it = MieAttenuation.begin();
   	    it!=MieAttenuation.end(); it++ ){
        double it_attenuation = *it;
        MieSactteringLength[iNumAir-1]=it_attenuation*m;  // attenuation length ( m )
        iNumAir--;
   }

   G4MaterialPropertiesTable* AIRMPT = new G4MaterialPropertiesTable();
   AIRMPT->AddProperty("RINDEX", EphotonAir, RefractiveAir, NumAir );
   AIRMPT->AddProperty("RAYLEIGH", EphotonAir, RayleighSactteringLength, NumAir);
   AIRMPT->AddProperty("MIEHG", EphotonAir, MieSactteringLength, NumAir);
   AIRMPT->AddConstProperty("MIEHG_FORWARD",  Mie_forward_g   );
   AIRMPT->AddConstProperty("MIEHG_BACKWARD", Mie_barckward_g );
   AIRMPT->AddConstProperty("MIEHG_FORWARD_RATIO", Mie_r );

   /*
   G4double AirScint[NumAir];
   for( int i=0; i<NumAir; i++ ){
   	    AirScint[i] = 1.0;
   }
   AIRMPT->AddProperty("FASTCOMPONENT", EphotonAir, AirScint, NumAir 
  ->SetSpline(true);

   AIRMPT->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
   AIRMPT->AddConstProperty("RESOLUTIONSCALE",1.0);
   AIRMPT->AddConstProperty("FASTTIMECONSTANT",1.*ns);
   AIRMPT->AddConstProperty("YIELDRATIO",0.8);
  */

   Air->SetMaterialPropertiesTable(AIRMPT);   
   //-------------------------------------------------------------------------------------- 

  // Vaccum ============================================================
  G4Material* Vacuum = new G4Material("Vacuum", density=universe_mean_density, nel=1,
				      //kStateGas, 2.73*kelvin, 3.e-18*pascal );
                                        kStateGas, 2.73*kelvin, 1.e-25*pascal );
  Vacuum->AddElement( N, 100.000*perCent );
 
  // ELS Container ====================================================
  G4Material* MaterialContainer = new G4Material("Container", density=7.874*g/cm3, nel=1 );
  MaterialContainer->AddElement( Fe, 100.000*perCent );  

  // Cover Box ========================================================
  G4Material* MaterialCoverBox = new G4Material("Container", density=7.874*g/cm3, nel=1 );
  MaterialCoverBox->AddElement( Fe, 100.000*perCent );

 //  Base A & B in 40FT Container ======================================
  G4Material* MaterialBase = new G4Material("Base", density=7.874*g/cm3, nel=1 );
  MaterialBase->AddElement( Fe, 100.000*perCent );  

  // Cu Collimeter =====================================================
  G4Material* MaterialCuCollimeter = new G4Material("CuCollimeter", density=8.920*g/cm3,  nel=1 );
  MaterialCuCollimeter->AddElement( Cu, 100.000*perCent );

  // Beam Duct SUS304 ==================================================
  G4Material* SUS304 = new G4Material("SUS304", density=7.90*g/cm3, nel=8 );
  SUS304->AddElement( C,   0.080*perCent );
  SUS304->AddElement( Si,  1.000*perCent );
  SUS304->AddElement( Mn,  2.000*perCent );
  SUS304->AddElement( P,   0.045*perCent );
  SUS304->AddElement( S,   0.030*perCent );
  SUS304->AddElement( Ni, 10.500*perCent );
  SUS304->AddElement( Cr, 20.000*perCent );
  SUS304->AddElement( Fe, 66.345*perCent );
  
  // York ====================================================
  G4Material* MaterialYork = new G4Material("York", density=7.874*g/cm3, nel=1 );
  MaterialYork->AddElement( Fe, 100.000*perCent );  

  // Coil ==============================================
  G4Material* MaterialCoil = new G4Material("Coil", density=8.920*g/cm3,  nel=1 );
  MaterialCoil->AddElement( Cu, 100.000*perCent );
  
  // Collimeter  ==============================================
  G4Material* MaterialCollimeter = new G4Material("TaCollimeter", density=16.665*g/cm3,  nel=1 );
  MaterialCollimeter->AddElement( Ta, 100.000*perCent );
     
  // OutputWindow ==============================================
  G4Material*  MaterialtTiwindow = new G4Material("Tiwindow", density=4.507*g/cm3,  nel=1 );
  MaterialtTiwindow->AddElement( Ti,  100.000*perCent );

  // Faraday Cup ==============================================
  // AL Cylinder
  G4Material*  MaterialtFDAlCylinder = new G4Material("FDAlCylinder", density=2.700*g/cm3,  nel=1 );
  MaterialtFDAlCylinder->AddElement( Al,  100.000*perCent );
  // FD Body
  G4Material* MaterialtFDPbBody =  new G4Material("FCPbBody",  density=11.35*g/cm3,  nel=1 );
  MaterialtFDPbBody->AddElement( Pb, 100.000*perCent );

  G4Material* MaterialtFDCBody =  new G4Material("FCCBody",  density=1.60*g/cm3,  nel=1 );
  MaterialtFDCBody->AddElement( C, 100.000*perCent );

  // FC2 Front/Body
  // --> FC4 
  // Creat Faraday Cup 4 in 2012.10.02
  // Geometry was fixed except position         

  // Tough-Pitch Copper(TPC)
  //G4Material* MaterialToughPitchCopper =  new G4Material("TPC",  density=8.89*g/cm3,  nel=1 );
  G4Material* MaterialToughPitchCopper =  new G4Material("TPC",  density=8.960*g/cm3,  nel=1 );
  MaterialToughPitchCopper->AddElement( Cu, 100.000*perCent );

  // Oxgen-Free Copper
  //G4Material* MaterialOxgenFreeCopper =  new G4Material("OFC",  density=8.89*g/cm3,  nel=1 );
  G4Material* MaterialOxgenFreeCopper =  new G4Material("OFC",  density=8.960*g/cm3,  nel=1 );
  MaterialOxgenFreeCopper->AddElement( Cu, 100.000*perCent );

  // MACOR = Machinable Glass Ceramics
  G4Material* MaterialMacor =  new G4Material("MACOR",  density=2.52*g/cm3,  nel=6 );
  MaterialMacor->AddMaterial( SiO2,  46.0*perCent );
  MaterialMacor->AddMaterial( MgO,   17.0*perCent );
  MaterialMacor->AddMaterial( Al2O3, 16.0*perCent );
  MaterialMacor->AddMaterial( K2O,   10.0*perCent );
  MaterialMacor->AddMaterial( B2O3,   7.0*perCent );
  MaterialMacor->AddElement( F,       4.0*perCent );
  
  // Pure Al
  G4Material*  MaterialtPureAl = new G4Material("PureAl", density=2.70*g/cm3,  nel=1 );
  MaterialtPureAl->AddElement( Al, 100.0*perCent );
  
  // AL5052
  G4Material*  MaterialtA5052 = new G4Material("A5052", density=2.68*g/cm3,  nel=2 );
  MaterialtA5052->AddElement( Mg,  2.5*perCent );
  MaterialtA5052->AddElement( Al, 97.5*perCent );

  // Material for Faraday Cup 5
  G4Material*  MaterialTi = new G4Material("Ti", density=4.507*g/cm3,  nel=1 );
  MaterialTi->AddElement( Ti,  100.000*perCent );
  G4Material*  MaterialAl = new G4Material("Al", density=2.700*g/cm3,  nel=1 );
  MaterialAl->AddElement( Al,  100.000*perCent );
  G4Material* MaterialBeCu =  new G4Material("BeCu",  density=8.36*g/cm3,  nel=2 );
  MaterialBeCu->AddElement( Cu, 98.000*perCent );
  MaterialBeCu->AddElement( Be, 2.000*perCent );

 // Screen Monitor 3 ( modified in 2012.01.26, added Cr2O3  )
  G4Material* MaterialSM3 = new G4Material("SM3", density=3.97625*g/cm3,  nel=2 ); 
  MaterialSM3->AddMaterial( Al2O3, 99.5*perCent );
  MaterialSM3->AddMaterial( Cr2O3,  0.5*perCent );

  // Beam Attenuator 
  G4Material* MaterialBeamAttenuator;
  if( BeamAttenuator_material == 1 ){
    MaterialBeamAttenuator = new G4Material( "BeamAttenuator", density=8.920*g/cm3,  nel=1 );
    MaterialBeamAttenuator->AddElement( Cu, 100.000*perCent );
  }else if( BeamAttenuator_material == 2 ){
    MaterialBeamAttenuator = new G4Material( "BeamAttenuator", density=2.70*g/cm3,  nel=1 );
    MaterialBeamAttenuator->AddElement( Al, 100.000*perCent );
  }else if( BeamAttenuator_material == 3 ){
    MaterialBeamAttenuator = new G4Material( "BeamAttenuator", density=0.953*g/cm3,  nel=1 );
    MaterialBeamAttenuator->AddMaterial( C2H4, 100.000*perCent );
  }else if( BeamAttenuator_material == 4 ){
    MaterialBeamAttenuator = new G4Material( "BeamAttenuator", density=16.665*g/cm3,  nel=1 );
    MaterialBeamAttenuator->AddElement( Ta, 100.000*perCent );
  }else if( BeamAttenuator_material == 5 ){
    MaterialBeamAttenuator = new G4Material( "BeamAttenuator", density=11.35*g/cm3,  nel=1 );
    MaterialBeamAttenuator->AddElement( Pb, 100.000*perCent );
  }

  // Pb Collimator added in 2012.09.25
  //G4Material *MaterialPbCollimator = new G4Material("PbCollimator", density=11.35*g/cm3,  nel=1 );
  //MaterialPbCollimator->AddElement( Pb, 100.000*perCent );
   G4Material *MaterialPbCollimator = new G4Material("PbCollimator", density=16.665*g/cm3,  nel=1 );
   MaterialPbCollimator->AddElement( Ta, 100.000*perCent );

  // Pb Block ==============================================
  G4Material* MaterialPbBlock = new G4Material("PbBlock", density=11.35*g/cm3,  nel=1 );
  MaterialPbBlock->AddElement( Pb, 100.000*perCent );

  // Concrete Block ========================================
  G4Material* MaterialConcrete = new G4Material("Concrete", density=2.3*g/cm3, nel=10 );
  MaterialConcrete->AddElement( H,   0.56*perCent );
  MaterialConcrete->AddElement( O,  49.83*perCent );
  MaterialConcrete->AddElement( Si, 31.58*perCent );
  MaterialConcrete->AddElement( Al,  4.56*perCent );
  MaterialConcrete->AddElement( Ca,  8.26*perCent );
  MaterialConcrete->AddElement( S,   0.12*perCent );
  MaterialConcrete->AddElement( Fe,  1.22*perCent );
  MaterialConcrete->AddElement( Mg,  0.24*perCent );
  MaterialConcrete->AddElement( Na,  1.71*perCent );
  MaterialConcrete->AddElement( K,   1.92*perCent );

  //--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------   

    G4ThreeVector G4PositionTemp = G4ThreeVector( 0.0*mm, 0.0*mm, 0.0*mm );

    G4RotationMatrix* G4RotationTemp = new G4RotationMatrix;
     G4RotationTemp->rotateX( 0.0*deg );
     G4RotationTemp->rotateY( 0.0*deg );
     G4RotationTemp->rotateZ( 0.0*deg );
    
  //------------------------------------------------
  // World Volume = Experimental Hole ( Air Space )
  //------------------------------------------------   
  G4ThreeVector PExperimentalHole = G4ThreeVector(0,0,0);
  //solidWorld = new G4Box("world", SizeOfWorld[0]*0.5*m, SizeOfWorld[1]*0.5*m, SizeOfWorld[2]*0.5*m );
  solidWorld = new G4Box("world", SizeOfWorld[0]*m, SizeOfWorld[1]*m, SizeOfWorld[2]*0.5*m );
  logicWorld = new G4LogicalVolume( solidWorld, Air, "World", 0, 0, 0);

  G4VisAttributes* WorldVisAtt = new G4VisAttributes(true, G4Colour(1.0,1.0,1.0));
  WorldVisAtt->SetForceWireframe(true);
  logicWorld->SetVisAttributes(WorldVisAtt);
  physiWorld = new G4PVPlacement(0, PExperimentalHole, logicWorld, "World", 0, false, 0 );       
  
  //------------------------------------------------
  // Ground Plane 
  //------------------------------------------------   
  PositionOfGroundPlane[0] = 0.0;
  PositionOfGroundPlane[1] = 0.0;
  //PositionOfGroundPlane[2] = (-1)*(SizeOfGroundPlane[2]*0.5);

  RotationOfGroundPlane[0] = 0.0;
  RotationOfGroundPlane[1] = 0.0;
  RotationOfGroundPlane[2] = 0.0;

  solidGroundPlane = new G4Box("Ground_Plane", 
			       SizeOfGroundPlane[0]*0.5*m ,
                               SizeOfGroundPlane[1]*0.5*m ,
                               SizeOfGroundPlane[2]*0.5*mm );
                                                          
  G4ThreeVector G4PositionOfGroundPlane
    = G4ThreeVector( PositionOfGroundPlane[0]*mm,
                     PositionOfGroundPlane[1]*mm,
                     PositionOfGroundPlane[2]*mm );

  G4RotationMatrix* G4RotationOfGroundPlane = new G4RotationMatrix;
  G4RotationOfGroundPlane->rotateX( RotationOfGroundPlane[0]*deg );
  G4RotationOfGroundPlane->rotateY( RotationOfGroundPlane[1]*deg );
  G4RotationOfGroundPlane->rotateZ( RotationOfGroundPlane[2]*deg );

  G4Transform3D Trans3DGroundPlane( *G4RotationOfGroundPlane, G4PositionOfGroundPlane );

  logicGroundPlane = new G4LogicalVolume( solidGroundPlane,
					  Air, 
					  "Ground",
					  0,0,0 );

   physiGroundPlane = new G4PVPlacement( Trans3DGroundPlane,
                                         logicGroundPlane,
                                         "Ground",
                                         logicWorld,
                                         false, 0 );

  //------------------------------------------------
  // Concrete Pad  
  //------------------------------------------------   
   solidConcretepad = new G4Box( "Concrete Pad",
                                 elsparam->sizeConcretepad(0)*0.5*mm, 
                                 elsparam->sizeConcretepad(1)*0.5*mm,
                                 elsparam->sizeConcretepad(2)*0.5*mm );

   G4ThreeVector G4PositionOfConcretepad
     = G4ThreeVector( elsparam->positionConcretepad(0)*mm,
                      elsparam->positionConcretepad(1)*mm,
		      elsparam->positionConcretepad(2)*mm );

  G4RotationMatrix* G4RotationOfConcretepad = new G4RotationMatrix;
    G4RotationOfConcretepad->rotateX( elsparam->rotationConcretepad(0)*deg );
    G4RotationOfConcretepad->rotateY( elsparam->rotationConcretepad(1)*deg );
    G4RotationOfConcretepad->rotateZ( elsparam->rotationConcretepad(2)*deg );

  G4Transform3D Trans3DConcretepad( *G4RotationOfConcretepad, 
	  			     G4PositionOfConcretepad );

  logicConcretepad = new G4LogicalVolume( solidConcretepad,
                                          MaterialConcrete,
                                          "Concrete Pad",
                                          0,0,0 );

  G4VisAttributes* ConcretePadVisAtt = new G4VisAttributes(G4Colour(0.2,0.2,0.2));  
  ConcretePadVisAtt->SetForceWireframe(true);
  ConcretePadVisAtt->SetForceSolid(true);
  logicConcretepad->SetVisAttributes(ConcretePadVisAtt);

  physiConcretepad = new G4PVPlacement( Trans3DConcretepad, 
                                        logicConcretepad,
                                        "Concrete Pad",
                                        logicWorld,
                                        false, 0 );  

  //-------------------------------------------------
  // ELS Container 
  //-------------------------------------------------
  G4VisAttributes* ELSVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.0,1.0));
  ELSVisAtt->SetForceWireframe(true);
  //ELSVisAtt->SetForceSolid(true);
  
  solidELSContainer_outer = new G4Box("Container_Outer", 
				      elsparam->sizeOuterELSContainer(0)*0.5*mm,
                                      elsparam->sizeOuterELSContainer(1)*0.5*mm,
				      elsparam->sizeOuterELSContainer(2)*0.5*mm );

  solidELSContainer_inner = new G4Box("Container_Inner",
				      elsparam->sizeInnerELSContainer(0)*0.5*mm,
				      elsparam->sizeInnerELSContainer(1)*0.5*mm,
				      elsparam->sizeInnerELSContainer(2)*0.5*mm );

  solidELSContainer_hole = new G4Box("Container_hole",
				     elsparam->sizeELSinjectionhole(0)*0.5*mm,
                                     elsparam->sizeELSinjectionhole(1)*0.5*mm,
				     elsparam->sizeELSinjectionhole(2)*0.5*mm );

  G4ThreeVector G4PositionOfELSContainer_Outer
    = G4ThreeVector( elsparam->positionOuterELSContainer(0)*mm,
                     elsparam->positionOuterELSContainer(1)*mm,
		     elsparam->positionOuterELSContainer(2)*mm );
  G4ThreeVector G4PositionOfELSContainer_Inner
    = G4ThreeVector( elsparam->positionInnerELSContainer(0)*mm,
                     elsparam->positionInnerELSContainer(1)*mm,
		     elsparam->positionInnerELSContainer(2)*mm );

  G4ThreeVector G4PositionOfELSContainer_Hole
    = G4ThreeVector( elsparam->positionELSinjectionhole(0)*mm,
                     elsparam->positionELSinjectionhole(1)*mm,
                     elsparam->positionELSinjectionhole(2)*mm ); 

  solidELSContainer = new G4SubtractionSolid("ELS Container",
                                              solidELSContainer_outer, 
                                              solidELSContainer_inner,
                                              0,
                                              G4PositionOfELSContainer_Inner -
                                              G4PositionOfELSContainer_Outer );  

  G4ThreeVector G4PositionOfELSContainer 
     = G4ThreeVector( elsparam->positionELSContainer(0)*mm,
                      elsparam->positionELSContainer(1)*mm,
		      elsparam->positionELSContainer(2)*mm );  

  G4RotationMatrix* G4RotationOfELSContainer = new G4RotationMatrix;
  G4RotationOfELSContainer->rotateX( elsparam->rotationELSContainer(0)*deg );
  G4RotationOfELSContainer->rotateY( elsparam->rotationELSContainer(1)*deg );
  G4RotationOfELSContainer->rotateZ( elsparam->rotationELSContainer(2)*deg );

  G4Transform3D Trans3DELSContainer( *G4RotationOfELSContainer,
                                      G4PositionOfELSContainer );

   solidELSContainer_withhole = new G4SubtractionSolid("ELS Container",
						       solidELSContainer,
                                                       solidELSContainer_hole,
                                                       0,
                                                       G4PositionOfELSContainer_Hole -
						       G4PositionOfELSContainer );

   logicELSContainer = new G4LogicalVolume( solidELSContainer_withhole, //solidELSContainer,
                                            MaterialContainer, 
                                            "ELS Container",
                                            0,0,0 );

   logicELSContainer->SetVisAttributes(ELSVisAtt);

   physiELSContainer = new G4PVPlacement( Trans3DELSContainer, 
                                          logicELSContainer,
                                          "ELS Container",
                                          logicWorld,
                                          false, 0 );  

  //-------------------------------------------------
  // ELS Cover Box
  //-------------------------------------------------
  G4VisAttributes* CoverBoxVisAtt = new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));
  CoverBoxVisAtt->SetForceWireframe(true);

  solidCoverBox_outer =  new G4Box("CoverBox Outer",
			  	    elsparam->sizeOuterCoverbox(0)*0.5*mm,
				    elsparam->sizeOuterCoverbox(1)*0.5*mm,
				    elsparam->sizeOuterCoverbox(2)*0.5*mm );

  solidCoverBox_inner =  new G4Box("CoverBox Inner",
				   elsparam->sizeInnerCoverbox(0)*0.5*mm,
				   elsparam->sizeInnerCoverbox(1)*0.5*mm,
				   elsparam->sizeInnerCoverbox(2)*0.5*mm );

  G4ThreeVector G4PositionOfCoverBox_outer 
    =  G4ThreeVector( elsparam->positionOuterCoverbox(0)*mm,
		      elsparam->positionOuterCoverbox(1)*mm,
		      elsparam->positionOuterCoverbox(2)*mm );

   G4ThreeVector G4PositionOfCoverBox_inner
     =  G4ThreeVector( elsparam->positionInnerCoverbox(0)*mm,
		       elsparam->positionInnerCoverbox(1)*mm,
		       elsparam->positionInnerCoverbox(2)*mm );

   solidsolidCoverBox = new G4SubtractionSolid("ELS CoverBox",
					        solidCoverBox_outer,
		 			        solidCoverBox_inner,
                                                0,
                                                G4PositionOfCoverBox_inner -
                                                G4PositionOfCoverBox_outer );


   solidCoverBox_board =  new G4Box("CoverBox Board", 
 				     elsparam->sizeCoverboxBoard(0)*0.5*mm,
                                     elsparam->sizeCoverboxBoard(1)*0.5*mm,
				     elsparam->sizeCoverboxBoard(2)*0.5*mm );

   solidCoverBox_boardhole = new G4Tubs( "CoverBox Boardhole",
					 elsparam->rminCoverboxBoardHole()*mm, 
                                         elsparam->rmaxCoverboxBoardHole()*mm,
                                         elsparam->lCoverboxBoardHole()*0.5*mm,
					 0.*deg, 360.0*deg );

    G4ThreeVector G4PositionOfCoverBox_board
      = G4ThreeVector( elsparam->positionCoverboxBoard(0)*mm,
		       elsparam->positionCoverboxBoard(1)*mm,
		       elsparam->positionCoverboxBoard(2)*mm );

    G4ThreeVector G4PositionOfCoverBox_boardhole
      = G4ThreeVector( elsparam->positionCoverboxBoardHole(0)*mm,
		       elsparam->positionCoverboxBoardHole(1)*mm,
		       elsparam->positionCoverboxBoardHole(2)*mm );
      
    solidCoverBox_board_withhole  = new G4SubtractionSolid("ELS CoverBoxBoard",
							   solidCoverBox_board,
							   solidCoverBox_boardhole,
							   0,
                                                           G4PositionOfCoverBox_boardhole -
							   G4PositionOfCoverBox_board );


    G4RotationMatrix* G4RotationOfCoverBox = new G4RotationMatrix;
    G4RotationOfCoverBox->rotateX( 0.0*deg );
    G4RotationOfCoverBox->rotateY( 0.0*deg );
    G4RotationOfCoverBox->rotateZ( 0.0*deg );

    G4Transform3D Trans3DCoverBox( *G4RotationOfCoverBox, G4PositionOfCoverBox_outer );
    G4Transform3D Trans3DCoverBoxBoard( *G4RotationOfCoverBox, G4PositionOfCoverBox_board );
 
    logicCoverBox = new G4LogicalVolume( solidsolidCoverBox,
                                         MaterialCoverBox,
                                         "ELS CoverBox",
					 0,0,0 );    

    logicCoverBox->SetVisAttributes(CoverBoxVisAtt);

    physiCoverBox = new G4PVPlacement(  Trans3DCoverBox,
					logicCoverBox, 
                                        "ELS CoverBox",
					logicWorld,
					false, 0 );

    logicCoverBoxBoard = new G4LogicalVolume( solidCoverBox_board_withhole,
					      MaterialCoverBox,
					      "ELS CoverBox",
					      0,0,0 );
    
    logicCoverBoxBoard->SetVisAttributes(CoverBoxVisAtt);
     
    physiCoverBoxBoard = new G4PVPlacement(  Trans3DCoverBoxBoard,
					     logicCoverBoxBoard,
					     "ELS CoverBox",
					     logicWorld,
					     false, 0 );

  //-------------------------------------------------
  // Base A & B & C in ELS Container
  //-------------------------------------------------
  G4VisAttributes* BaseVisAtt = new G4VisAttributes(true,G4Colour(0.0,1.0,0.0));
    BaseVisAtt->SetForceWireframe(true);
    //BaseVisAtt->SetForceSolid(true);

  for( int i=0; i<3; i++ ){

  solidBase[i] = new G4Box( "Base", 
			    elsparam->sizeBase(i,0)*0.5*mm,
			    elsparam->sizeBase(i,1)*0.5*mm,
			    elsparam->sizeBase(i,2)*0.5*mm );

   G4ThreeVector G4PositionOfBase 
     = G4ThreeVector( elsparam->positionBase(i,0)*mm,
                      elsparam->positionBase(i,1)*mm,
		      elsparam->positionBase(i,2)*mm );   

    G4RotationMatrix* G4RotationOfBase = new G4RotationMatrix;
     G4RotationOfBase->rotateX( elsparam->rotationBase(i,0)*deg );
     G4RotationOfBase->rotateY( elsparam->rotationBase(i,1)*deg );
     G4RotationOfBase->rotateZ( elsparam->rotationBase(i,2)*deg );

    G4Transform3D Trans3DBase( *G4RotationOfBase, G4PositionOfBase );

    logicBase[i] = new G4LogicalVolume( solidBase[i], MaterialBase, "Base", 0, 0, 0 );
    logicBase[i]->SetVisAttributes(BaseVisAtt);

    physiBase[i] = new G4PVPlacement( Trans3DBase, logicBase[i], "Base", logicWorld, false, 0 );

  }

  //-------------------------------------------------
  // Beam Line Components
  //-------------------------------------------------
  //-------------------------------------------------
  // Tube Component
  //-------------------------------------------------
  G4VisAttributes* BLVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,1.0));  
     BLVisAtt->SetForceWireframe(true);
     BLVisAtt->SetForceSolid(true);

  G4VisAttributes* BLFFVisAtt = new G4VisAttributes(true,G4Colour(1.0,0.5,0.5));  
     BLFFVisAtt->SetForceWireframe(true);
     BLFFVisAtt->SetForceSolid(true);

  G4VisAttributes* BLBFVisAtt = new G4VisAttributes(true,G4Colour(0.5,1.0,0.5));  
     BLBFVisAtt->SetForceWireframe(true);
     BLBFVisAtt->SetForceSolid(true);

  for( int i=0; i < elsparam->numBLcomp(); i++ ){
   
    // Beam Duct & Vacuum
    solidBL_tubs[i] = new G4Tubs( "BL", elsparam->rminBL(i)*mm, elsparam->rmaxBL(i)*mm,
                                        elsparam->lBL(i)*0.5*mm, 0.*deg, 360.0*deg );  
    solidBL_vac_tubs[i] = new G4Tubs( "BL VAC", 0.*mm, elsparam->rminBL(i)*mm,
                                                elsparam->lBL(i)*0.5*mm, 0.*deg, 360.0*deg );  
    // Beam Fowward Frange & Vacuum
    solidBL_FF_tubs[i] = new G4Tubs( "BLFF", elsparam->rminBLff(i)*mm, elsparam->rmaxBLff(i)*mm,
                                             elsparam->lBLff(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidBL_FF_vac_tubs[i] = new G4Tubs( "BLFF VAC", 0.*mm, elsparam->rminBLff(i)*mm,
				                     elsparam->lBLff(i)*0.5*mm, 0.*deg, 360.0*deg );
    // Beam Backward Frange & Vacuum
    solidBL_BF_tubs[i] = new G4Tubs( "BLBF", elsparam->rminBLbf(i)*mm, elsparam->rmaxBLbf(i)*mm,
                                             elsparam->lBLbf(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidBL_BF_vac_tubs[i] = new G4Tubs( "BLBF VAC", 0.*mm, elsparam->rminBLbf(i)*mm,
 				                     elsparam->lBLbf(i)*0.5*mm, 0.*deg, 360.0*deg );

   G4ThreeVector G4PositionOfBL
     = G4ThreeVector( elsparam->positionBL(i,0)*mm, 
                      elsparam->positionBL(i,1)*mm, 
                      elsparam->positionBL(i,2)*mm ); 

   G4ThreeVector G4PositionOfBLFF
     = G4ThreeVector( elsparam->positionBLff(i,0)*mm,
		      elsparam->positionBLff(i,1)*mm,
		      elsparam->positionBLff(i,2)*mm );

   G4ThreeVector G4PositionOfBLBF
     = G4ThreeVector( elsparam->positionBLbf(i,0)*mm,
	  	      elsparam->positionBLbf(i,1)*mm,
		      elsparam->positionBLbf(i,2)*mm );

   G4RotationMatrix* G4RotationOfBL = new G4RotationMatrix;
    G4RotationOfBL->rotateX( elsparam->rotationBL(i,0)*deg );
    G4RotationOfBL->rotateY( elsparam->rotationBL(i,1)*deg );
    G4RotationOfBL->rotateZ( elsparam->rotationBL(i,2)*deg );

   G4Transform3D Trans3DBL(   *G4RotationOfBL, G4PositionOfBL   );
   G4Transform3D Trans3DBLFF( *G4RotationOfBL, G4PositionOfBLFF );
   G4Transform3D Trans3DBLBF( *G4RotationOfBL, G4PositionOfBLBF );

   logicBL_tubs[i] = new G4LogicalVolume( solidBL_tubs[i], 
                                          SUS304, "BL", 0, 0, 0 );
   logicBL_tubs[i]->SetVisAttributes( BLVisAtt );
   physiBL_tubs[i] = new G4PVPlacement( Trans3DBL, logicBL_tubs[i], 
                                        "BL", logicWorld, false, 0 );

   logicBL_vac_tubs[i] = new G4LogicalVolume( solidBL_vac_tubs[i], 
                                              Vacuum, "BL VAC", 0, 0, 0 );
   //logicBL_vac_tubs[i]->SetVisAttributes( BLVisAtt ); 
   physiBL_vac_tubs[i] = new G4PVPlacement( Trans3DBL, logicBL_vac_tubs[i], 
                                            "BL VAC", logicWorld, false, 0 );

   logicBL_FF_tubs[i] = new G4LogicalVolume( solidBL_FF_tubs[i], 
                                             SUS304, "BL FF", 0, 0, 0 );
   logicBL_FF_tubs[i]->SetVisAttributes( BLFFVisAtt );
   physiBL_FF_tubs[i] = new G4PVPlacement( Trans3DBLFF, logicBL_FF_tubs[i], 
                                           "BL FF", logicWorld, false, 0 );

   logicBL_FF_vac_tubs[i] = new G4LogicalVolume( solidBL_FF_vac_tubs[i], 
                                                 Vacuum, "BL FF VAC", 0, 0, 0 );
   //logicBL_FF_vac_tubs[i]->SetVisAttributes( BLFFVisAtt );
   physiBL_FF_vac_tubs[i] = new G4PVPlacement( Trans3DBLFF, logicBL_FF_vac_tubs[i],
                                               "BL FF VAC", logicWorld, false, 0 );

   logicBL_BF_tubs[i] = new G4LogicalVolume( solidBL_BF_tubs[i],
                                             SUS304, "BL BF", 0, 0, 0 );
   logicBL_BF_tubs[i]->SetVisAttributes( BLBFVisAtt );
   physiBL_BF_tubs[i] = new G4PVPlacement( Trans3DBLBF, logicBL_BF_tubs[i], 
                                           "BL BF", logicWorld, false, 0 );
   logicBL_BF_vac_tubs[i] = new G4LogicalVolume( solidBL_BF_vac_tubs[i], 
                                                 Vacuum, "BL BF VAC", 0, 0, 0 );
   //logicBL_BF_vac_tubs[i]->SetVisAttributes( BLBFVisAtt );
   physiBL_BF_vac_tubs[i] = new G4PVPlacement( Trans3DBLBF, logicBL_BF_vac_tubs[i],
                                               "BL BF VAC", logicWorld, false, 0 );
  }

    //--------------------------------------------------
    // Blank of EGUN
    //--------------------------------------------------
     G4VisAttributes* EGUNBLANKVisAtt = new G4VisAttributes(true,G4Colour(1.0,0.2,0.2));  
      EGUNBLANKVisAtt->SetForceWireframe(true);
      EGUNBLANKVisAtt->SetForceSolid(true);

      solidEGUNBlankFrange = new G4Tubs("EGUNBlank",
                                    elsparam->rminEGUNBlank()*mm, elsparam->rmaxEGUNBlank()*mm, 
                                    elsparam->lEGUNBlank()*0.5*mm, 0.*deg, 360.*deg ); 
      logicEGUNBlankFrange = new G4LogicalVolume( solidEGUNBlankFrange, SUS304, "EGUNBlank",  0, 0, 0 );
      logicEGUNBlankFrange->SetVisAttributes(EGUNBLANKVisAtt);
       
      G4ThreeVector G4PositionOfEGUNBlank
        = G4ThreeVector( elsparam->positionEGUNBlank(0)*mm,
                         elsparam->positionEGUNBlank(1)*mm,
			 elsparam->positionEGUNBlank(2)*mm );
      G4RotationMatrix* G4RotationOfEGUNBlank = new G4RotationMatrix;
          G4RotationOfEGUNBlank->rotateX( elsparam->rotationEGUNBlank(0)*deg );
          G4RotationOfEGUNBlank->rotateY( elsparam->rotationEGUNBlank(1)*deg );
          G4RotationOfEGUNBlank->rotateZ( elsparam->rotationEGUNBlank(2)*deg );
      G4Transform3D Trans3DEGUNBlank( *G4RotationOfEGUNBlank, G4PositionOfEGUNBlank );
      
      physiEGUNBlankFrange = new G4PVPlacement( Trans3DEGUNBlank, logicEGUNBlankFrange ,
					    "EGUNBlank", logicWorld , false, 0 );

    //--------------------------------------------------
    // Anode of EGUN
    //--------------------------------------------------
     G4VisAttributes* ANODEVisAtt = new G4VisAttributes(true,G4Colour(1.0,0.2,0.2));  
      ANODEVisAtt->SetForceWireframe(true);
      ANODEVisAtt->SetForceSolid(true);

      solidAnode = new G4Tubs("Anode",
                                    elsparam->rminAnode()*mm, elsparam->rmaxAnode()*mm, 
                                    elsparam->lAnode()*0.5*mm, 0.*deg, 360.*deg ); 
      logicAnode = new G4LogicalVolume( solidAnode, SUS304, "Anode",  0, 0, 0 );
      logicAnode->SetVisAttributes(ANODEVisAtt);
       
      G4ThreeVector G4PositionOfAnode
        = G4ThreeVector( elsparam->positionAnode(0)*mm,
                         elsparam->positionAnode(1)*mm,
			 elsparam->positionAnode(2)*mm );
      G4RotationMatrix* G4RotationOfAnode = new G4RotationMatrix;
          G4RotationOfAnode->rotateX( elsparam->rotationAnode(0)*deg );
          G4RotationOfAnode->rotateY( elsparam->rotationAnode(1)*deg );
          G4RotationOfAnode->rotateZ( elsparam->rotationAnode(2)*deg );
      G4Transform3D Trans3DAnode( *G4RotationOfAnode, G4PositionOfAnode );
      
      physiAnode = new G4PVPlacement( Trans3DAnode, logicAnode ,				     
                                     "Anode", logicWorld , false, 0 );
  
  //--------------------------------------------------
  // Screen Monitor
  //--------------------------------------------------
    G4VisAttributes* SMFVisAtt = new G4VisAttributes(true,G4Colour(0.3,1.0,0.3));  
      SMFVisAtt->SetForceWireframe(true);
      SMFVisAtt->SetForceSolid(true);
    G4VisAttributes* SMBVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,0.2));  
      SMBVisAtt->SetForceWireframe(true);
      SMBVisAtt->SetForceSolid(true);
    G4VisAttributes* SMDVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,0.4));  
      SMDVisAtt->SetForceWireframe(true);
      SMDVisAtt->SetForceSolid(true);

      for( int i=0; i<2; i++ ){

     // Frange 
     solidSM_FF_tubs[i] = new G4Tubs( "SM FF",
                                   elsparam->rminSMff(i)*mm, elsparam->rmaxSMff(i)*mm,
                                   elsparam->lSMff(i)*0.5*mm, 0.*deg, 360.0*deg );
     solidSM_FF_vac_tubs[i] = new G4Tubs( "SM FF Vac",
                                   0.0*mm, elsparam->rminSMff(i)*mm, 
                                   elsparam->lSMff(i)*0.5*mm, 0.*deg, 360.0*deg );
     solidSM_BF_tubs[i] = new G4Tubs( "SM BF",
                                   elsparam->rminSMbf(i)*mm, elsparam->rmaxSMbf(i)*mm,
                                   elsparam->lSMbf(i)*0.5*mm, 0.*deg, 360.0*deg );     
     solidSM_BF_vac_tubs[i] = new G4Tubs( "SM BF Vac",
                                   0.0*mm, elsparam->rminSMbf(i)*mm, 
                                   elsparam->lSMbf(i)*0.5*mm, 0.*deg, 360.0*deg );     

     G4ThreeVector G4PositionOfSMff
	 = G4ThreeVector( elsparam->positionSMff( i, 0 )*mm,
			  elsparam->positionSMff( i, 1 )*mm,
			  elsparam->positionSMff( i, 2 )*mm );
     G4ThreeVector G4PositionOfSMbf 
	 = G4ThreeVector( elsparam->positionSMbf( i, 0 )*mm,
		 	  elsparam->positionSMbf( i, 1 )*mm,
			  elsparam->positionSMbf( i, 2 )*mm );				

     G4RotationMatrix* G4RotationOfSMf = new G4RotationMatrix;     
	 G4RotationOfSMf->rotateX( elsparam->rotationSMf( i, 0 )*deg );	 
	 G4RotationOfSMf->rotateY( elsparam->rotationSMf( i, 1 )*deg );	 
	 G4RotationOfSMf->rotateZ( elsparam->rotationSMf( i, 2 )*deg );	 

      G4Transform3D Trans3DSMFF( *G4RotationOfSMf, G4PositionOfSMff );
      G4Transform3D Trans3DSMBF( *G4RotationOfSMf, G4PositionOfSMbf );

    logicSM_FF_tubs[i] = new G4LogicalVolume( solidSM_FF_tubs[i], SUS304, "SM FF", 0, 0, 0 );
     logicSM_FF_tubs[i]->SetVisAttributes( SMFVisAtt );
    logicSM_FF_vac_tubs[i] = new G4LogicalVolume( solidSM_FF_vac_tubs[i], Vacuum, "SM FF Vac", 0, 0, 0 );
     //logicSM_FF_vac_tubs[i]->SetVisAttributes( SMFVisAtt );

    logicSM_BF_tubs[i] = new G4LogicalVolume( solidSM_BF_tubs[i], SUS304, "SM BF", 0, 0, 0 );
     logicSM_BF_tubs[i]->SetVisAttributes( SMFVisAtt );
    logicSM_BF_vac_tubs[i] = new G4LogicalVolume( solidSM_BF_vac_tubs[i], Vacuum, "SM BF Vac", 0, 0, 0 );    
     //logicSM_BF_vac_tubs[i]->SetVisAttributes( SMFVisAtt );

    physiSM_FF_tubs[i] = new G4PVPlacement( Trans3DSMFF, logicSM_FF_tubs[i], 
                                            "SM FF", logicWorld, false, 0 );
    physiSM_FF_vac_tubs[i] = new G4PVPlacement( Trans3DSMFF, logicSM_FF_vac_tubs[i], 
                                            "SM FF Vac", logicWorld, false, 0 );

    physiSM_BF_tubs[i] = new G4PVPlacement( Trans3DSMBF, logicSM_BF_tubs[i], 
                                            "SM BF", logicWorld, false, 0 );
    physiSM_BF_vac_tubs[i] = new G4PVPlacement( Trans3DSMBF, logicSM_BF_vac_tubs[i], 

                                            "SM BF Vac", logicWorld, false, 0 );
    // Body & Duct
    solidSM_parts1_tubs[i] =  new G4Tubs( "SM parts1",
                                          elsparam->rminSMparts1(i)*mm, elsparam->rmaxSMparts1(i)*mm,
                                          elsparam->lSMparts1(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidSM_parts12_tubs[i] =  new G4Tubs( "SM parts12",
                                          0.0*mm, elsparam->rmaxSMparts1(i)*mm,
                                          elsparam->lSMparts1(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidSM_parts13_tubs[i] =  new G4Tubs( "SM parts13",
                                          0.0*mm, elsparam->rminSMparts1(i)*mm,
                                          elsparam->lSMparts1(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidSM_parts2_tubs[i] =  new G4Tubs( "SM parts2",
                                          elsparam->rminSMparts2(i)*mm, elsparam->rmaxSMparts2(i)*mm,
                                          elsparam->lSMparts2(i)*0.5*mm, 0.*deg, 360.0*deg );
    solidSM_parts3_tubs[i] =  new G4Tubs( "SM parts3",
                                          elsparam->rminSMparts3(i)*mm, elsparam->rmaxSMparts3(i)*mm,
                                          elsparam->lSMparts3(i)*0.5*mm, 0.*deg, 360.0*deg );

    G4ThreeVector G4PositionOfSMparts1
      = G4ThreeVector( elsparam->positionSMparts1( i, 0 )*mm,
		       elsparam->positionSMparts1( i, 1 )*mm,
		       elsparam->positionSMparts1( i, 2 )*mm );
    G4ThreeVector G4PositionOfSMparts2
      = G4ThreeVector( elsparam->positionSMparts2( i, 0 )*mm,
		       elsparam->positionSMparts2( i, 1 )*mm,
		       elsparam->positionSMparts2( i, 2 )*mm );
    G4ThreeVector G4PositionOfSMparts3
      = G4ThreeVector( elsparam->positionSMparts3( i, 0 )*mm,
                       elsparam->positionSMparts3( i, 1 )*mm,
                       elsparam->positionSMparts3( i, 2 )*mm );
   
    G4RotationMatrix* G4RotationOfSMcomp = new G4RotationMatrix;
      G4RotationOfSMcomp->rotateX( elsparam->rotationSMcomp( i, 0 )*deg );
      G4RotationOfSMcomp->rotateY( elsparam->rotationSMcomp( i, 1 )*deg );
      G4RotationOfSMcomp->rotateZ( elsparam->rotationSMcomp( i, 2 )*deg );

    G4RotationMatrix* G4RotationOfSMduct = new G4RotationMatrix;
      G4RotationOfSMduct->rotateX( elsparam->rotationSMduct( i, 0 )*deg );
      G4RotationOfSMduct->rotateY( elsparam->rotationSMduct( i, 1 )*deg );
      G4RotationOfSMduct->rotateZ( elsparam->rotationSMduct( i, 2 )*deg );

    G4RotationMatrix* G4RotationOfSMbody = new G4RotationMatrix;
      G4RotationOfSMbody->rotateX( elsparam->rotationSMbody( i, 0 )*deg );
      G4RotationOfSMbody->rotateY( elsparam->rotationSMbody( i, 1 )*deg );
      G4RotationOfSMbody->rotateZ( elsparam->rotationSMbody( i, 2 )*deg );

    G4Transform3D Trans3DSMDUCT( *G4RotationOfSMduct, G4PositionOfSMparts1 );    
    G4Transform3D Trans3DSMBODY( *G4RotationOfSMbody, G4PositionOfSMparts3 );

    solidSM_duct[i] = new G4SubtractionSolid( "SM Duct",
					      solidSM_parts1_tubs[i],
					      solidSM_parts2_tubs[i],
					      G4RotationOfSMcomp,
					      G4PositionOfSMparts2 - G4PositionOfSMparts1 );
 
    solidSM_duct_vac[i]= new G4SubtractionSolid( "SM Duct Vacuum",
						 solidSM_parts13_tubs[i],
						 solidSM_parts2_tubs[i],
						 G4RotationOfSMcomp,
						 G4PositionOfSMparts2 - G4PositionOfSMparts1 );

    solidSM_body[i] = new G4SubtractionSolid( "SM Body",
                                              solidSM_parts3_tubs[i],
                                              solidSM_parts12_tubs[i],
                                             G4RotationOfSMcomp,
                                              G4PositionOfSMparts1 - G4PositionOfSMparts3 );
    solidSM_body_vac[i] = new G4Tubs( "SM Body Vac",
                                      elsparam->rminSMparts2(i)*mm, elsparam->rmaxSMparts2(i)*mm,
                                      elsparam->lSMparts2(i)*0.5*mm, 0.*deg, 360.0*deg );
    
    logicSM_duct[i] = new G4LogicalVolume( solidSM_duct[i],
                                           SUS304, "SM Duct", 0, 0, 0 );
    logicSM_duct[i]->SetVisAttributes( SMDVisAtt );
    logicSM_duct_vac[i] = new G4LogicalVolume( solidSM_duct_vac[i],
	 				       Vacuum, "SM Duct Vac", 0, 0, 0 );
    //logicSM_duct_vac[i]->SetVisAttributes( SMDVisAtt );

    physiSM_duct[i] = new G4PVPlacement( Trans3DSMDUCT, logicSM_duct[i],
                                         "SM Duct", logicWorld, false, 0 );
    physiSM_duct_vac[i] = new G4PVPlacement(  Trans3DSMDUCT,logicSM_duct_vac[i],
					      "SM Duct Vac", logicWorld, false, 0 );     
    
    logicSM_body[i] = new G4LogicalVolume( solidSM_body[i],
                                           SUS304, "SM Body", 0, 0, 0 );    
    logicSM_body[i]->SetVisAttributes( SMBVisAtt );  
    logicSM_body_vac[i] =  new G4LogicalVolume( solidSM_body_vac[i], 
						Vacuum, "SM Body Vac", 0, 0, 0 );
    //logicSM_body_vac[i]->SetVisAttributes( SMBVisAtt );  
    physiSM_body[i] = new G4PVPlacement( Trans3DSMBODY, logicSM_body[i],
                                         "SM Body", logicWorld, false, 0 );
    physiSM_body_vac[i]= new G4PVPlacement( Trans3DSMBODY, logicSM_body_vac[i],
                                         "SM Body Vac", logicWorld, false, 0 );
   }

  //--------------------------------------------------
  // BM Duct
  //--------------------------------------------------

    //--------------------------------------------------
    //  Cu Collimeter 
    //--------------------------------------------------
      G4VisAttributes* CuCollimeterVisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.5, 0.5));
       CuCollimeterVisAtt->SetForceWireframe(true);
       CuCollimeterVisAtt->SetForceSolid(true);

       solidCuCollimeter_body = new G4Box( "CuCollimeter box",
					   elsparam->sizeCuCollimeter(0)*0.5*mm,
					   elsparam->sizeCuCollimeter(1)*0.5*mm,
					   elsparam->sizeCuCollimeter(2)*0.5*mm );
       solidCuCollimeter_hole = new G4Tubs( "CuCollimeter hole",
					   elsparam->rminCuCollimeterhole()*mm, 
					   CuCollimeter_Diameter/2.0*mm,
					   elsparam->lCuCollimeterhole()*0.5*mm,
                                            0.*deg, 360.0*deg );

       solidCuCollimeter_vacuum = new G4Tubs( "CuCollimeter vac",
					      elsparam->rminCuCollimeterhole()*mm,
					      CuCollimeter_Diameter/2.0*mm,
					      elsparam->lCuCollimeterhole()*0.5*mm,
                                              0.*deg, 360.0*deg );
        
    G4ThreeVector G4PositionOfCuCollimeter_body
      = G4ThreeVector( elsparam->positionCuCollimeter(0)*mm,
                       elsparam->positionCuCollimeter(1)*mm,
		       elsparam->positionCuCollimeter(2)*mm );
     G4ThreeVector G4PositionOfCuCollimeter_hole
       = G4ThreeVector( elsparam->positionCuCollimeterhole(0)*mm,
			elsparam->positionCuCollimeterhole(1)*mm,
			elsparam->positionCuCollimeterhole(2)*mm );

    G4RotationMatrix* G4RotationOfCuCollimeter_body = new G4RotationMatrix;
       G4RotationOfCuCollimeter_body->rotateX( elsparam->rotationCuCollimeter(0)*deg );
       G4RotationOfCuCollimeter_body->rotateY( elsparam->rotationCuCollimeter(1)*deg );
       G4RotationOfCuCollimeter_body->rotateZ( elsparam->rotationCuCollimeter(2)*deg );
       
    G4RotationMatrix* G4RotationOfCuCollimeter_hole = new G4RotationMatrix;
       G4RotationOfCuCollimeter_hole->rotateX( elsparam->rotationCuCollimeterhole(0)*deg );
       G4RotationOfCuCollimeter_hole->rotateY( elsparam->rotationCuCollimeterhole(1)*deg );
       G4RotationOfCuCollimeter_hole->rotateZ( elsparam->rotationCuCollimeterhole(2)*deg );

    G4Transform3D Trans3DCuCollimeter( *G4RotationOfCuCollimeter_body, G4PositionOfCuCollimeter_body );
    G4Transform3D Trans3DCuCollimeter_vacuum( *G4RotationOfCuCollimeter_hole, G4PositionOfCuCollimeter_hole );
    solidCuCollimeter = new G4SubtractionSolid( "CuCollimeter",
					        solidCuCollimeter_body,
					        solidCuCollimeter_hole,
					        G4RotationOfCuCollimeter_hole,
					        G4PositionOfCuCollimeter_hole - G4PositionOfCuCollimeter_body ); 

    logicCuCollimeter 
        = new G4LogicalVolume( solidCuCollimeter, MaterialCuCollimeter, "CuCollimeter", 0, 0, 0 );
    logicCuCollimeter->SetVisAttributes( CuCollimeterVisAtt );

    logicCuCollimeter_vacuum 
        = new G4LogicalVolume( solidCuCollimeter_vacuum, Vacuum, "CuCollimeter vac", 0, 0, 0 );
  
    physiCuCollimeter = new G4PVPlacement( Trans3DCuCollimeter, logicCuCollimeter, 
                                           "CuCollimeter", logicWorld, false, 0 );
    physiCuCollimeter_vacuum =  new G4PVPlacement( Trans3DCuCollimeter_vacuum, logicCuCollimeter_vacuum,
						   "CuCollimeter vac", logicWorld, false, 0 );
       
  //--------------------------------------------------
  //  BM Duct Frange
  //--------------------------------------------------
    G4VisAttributes* BMFVisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.1, 0.1));  
       BMFVisAtt->SetForceWireframe(true);
       BMFVisAtt->SetForceSolid(true);

    solidBMF_outer = new G4Tubs( "BMFOuter", 
                                  elsparam->rminBMf()*mm,  elsparam->rmaxBMf()*mm,
                                  elsparam->lBMf()*0.5*mm, 0.*deg, 360.0*deg );     
    solidBMF_inner  = new G4Box( "BMFInner",
				 elsparam->sizeInnerBMf(0)*0.5*mm,
                                 elsparam->sizeInnerBMf(1)*0.5*mm,
				 elsparam->sizeInnerBMf(2)*0.5*mm );

    solidBMF = new G4SubtractionSolid("BMF",
				      solidBMF_outer, solidBMF_inner, 
                                      G4RotationTemp, G4PositionTemp );

    logicBMF = new G4LogicalVolume( solidBMF, SUS304, "BMF", 0, 0, 0 );
    logicBMF->SetVisAttributes( BMFVisAtt );

    solidBMF_vac  = new G4Box( "BMF vac",
                                 elsparam->sizeInnerBMf(0)*0.5*mm,
                                 elsparam->sizeInnerBMf(1)*0.5*mm,
                                 elsparam->sizeInnerBMf(2)*0.5*mm );
    logicBMF_vac = new G4LogicalVolume( solidBMF_vac, Vacuum, "BMF vac", 0, 0, 0 );
    //logicBMF_vac->SetVisAttributes( BMFVisAtt );

    for( int i=0; i<3; i++ ){
       G4ThreeVector G4PositionOfBMF
	 = G4ThreeVector( elsparam->positionBMf( i, 0 )*mm,
                          elsparam->positionBMf( i, 1 )*mm,
			  elsparam->positionBMf( i, 2 )*mm );
       G4RotationMatrix* G4RotationOfBMF = new G4RotationMatrix;     
	 G4RotationOfBMF->rotateX( elsparam->rotationBMf( i, 0 )*deg );	 
         G4RotationOfBMF->rotateY( elsparam->rotationBMf( i, 1 )*deg );
	 G4RotationOfBMF->rotateZ( elsparam->rotationBMf( i, 2 )*deg );
       G4Transform3D Trans3DBMF( *G4RotationOfBMF, G4PositionOfBMF );

       physiBMF[i] = new G4PVPlacement( Trans3DBMF, logicBMF, "BMF", logicWorld, false, 0 );

       // Change by T.Shibata 2010.02.11
       // i=0 was omitted, because of Cu Collimeter was install into BM duct.
       if(i>0){   
          physiBMF_vac[i] = new G4PVPlacement( Trans3DBMF, logicBMF_vac, "BMF vac", logicWorld, false, 0 );
       }
    }

  //--------------------------------------------------
  //  BM Duct part-1-3
  //--------------------------------------------------
    G4VisAttributes* BMLVisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.2, 0.2));  
       BMLVisAtt->SetForceWireframe(true);
       BMLVisAtt->SetForceSolid(true);

  for( int i=0; i<3; i++ ){
    solidBMLduct_outer[i] = new G4Box("BMLouter",
                                      elsparam->sizeOuterBML(i, 0)*0.5*mm,
                                      elsparam->sizeOuterBML(i, 1)*0.5*mm,
                                      elsparam->sizeOuterBML(i, 2)*0.5*mm );
    solidBMLduct_inner[i] = new G4Box("BMLinner",
                                      elsparam->sizeInnerBML(i, 0)*0.5*mm,
                                      elsparam->sizeInnerBML(i, 1)*0.5*mm,
                                      elsparam->sizeInnerBML(i, 2)*0.5*mm );
    solidBMLduct_vac[i] = new G4Box("BMLinner vac",
                                      elsparam->sizeBMLvac(i, 0)*0.5*mm,
                                      elsparam->sizeBMLvac(i, 1)*0.5*mm,
                                      elsparam->sizeBMLvac(i, 2)*0.5*mm );
     G4ThreeVector G4PositionOfBML
       =  G4ThreeVector( elsparam->positionOuterBML(i, 0)*mm,
			 elsparam->positionOuterBML(i, 1)*mm,
			 elsparam->positionOuterBML(i, 2)*mm ); 
     G4RotationMatrix* G4RotationOfBML = new G4RotationMatrix;
     G4RotationOfBML->rotateX( 0. );
     G4RotationOfBML->rotateY( 0. );
     G4RotationOfBML->rotateZ( 0. );
     G4Transform3D Trans3DBML( *G4RotationOfBML, G4PositionOfBML );

     G4ThreeVector G4PositionOfBML_vac
       =  G4ThreeVector( elsparam->positionBMLvac(i, 0)*mm,
			 elsparam->positionBMLvac(i, 1)*mm,
			 elsparam->positionBMLvac(i, 2)*mm ); 
     G4RotationMatrix* G4RotationOfBML_vac = new G4RotationMatrix;
     G4RotationOfBML_vac->rotateX( 0. );
     G4RotationOfBML_vac->rotateY( 0. );
     G4RotationOfBML_vac->rotateZ( 0. );
     G4Transform3D Trans3DBML_vac( *G4RotationOfBML_vac, G4PositionOfBML_vac );

    solidBMLduct[i] = new G4SubtractionSolid("BML",
			  solidBMLduct_outer[i],
                          solidBMLduct_inner[i],
 	                  G4RotationTemp, G4PositionTemp );

    logicBMLduct[i] = new G4LogicalVolume( solidBMLduct[i], SUS304, "BML", 0, 0, 0 );
    logicBMLduct[i]->SetVisAttributes( BMLVisAtt );
 
    logicBMLduct_vac[i] = new G4LogicalVolume( solidBMLduct_vac[i], Vacuum, "BML vac", 0, 0, 0 );    
    //logicBMLduct_vac[i]->SetVisAttributes( BMLVisAtt );

    physiBMLduct[i] = new G4PVPlacement( Trans3DBML, 
                                         logicBMLduct[i],
                                         "BML", logicWorld, false, 0 );
    physiBMLduct_vac[i] = new G4PVPlacement( Trans3DBML_vac,
                                             logicBMLduct_vac[i],
                                             "BML vac", logicWorld, false, 0 );
  }
  //--------------------------------------------------
  // BM Duct parts 4 (1)-(5)
  //--------------------------------------------------
    G4VisAttributes* BMp41VisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.4, 0.4));  
       BMp41VisAtt->SetForceWireframe(true);
       BMp41VisAtt->SetForceSolid(true);
    G4VisAttributes* BMp42VisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.6, 0.6));  
       BMp42VisAtt->SetForceWireframe(true);
       BMp42VisAtt->SetForceSolid(true);

  solidBMmainSub1 = new G4Tubs( "BMmainSub1", 
                                elsparam->rminBMparts43()*mm, elsparam->rmaxBMparts43()*mm,
                                elsparam->lBMparts43()*0.5*mm, 0.0*deg, 360.0*deg ); 
  solidBMmainSub2 = new G4Box( "BMmainSub2", 
			       elsparam->sizeBMparts44(0)*0.5*mm,
                               elsparam->sizeBMparts44(1)*0.5*mm,
                               elsparam->sizeBMparts44(2)*0.5*mm );
  solidBMmainSub3 = new G4Tubs( "BMmainSub3",
                                elsparam->rminBMparts45()*mm, elsparam->rmaxBMparts45()*mm,
                                elsparam->lBMparts45()*0.5*mm, 0.0*deg, 360.0*deg );

   G4ThreeVector G4PositionOfBMmainSub1
     =  G4ThreeVector( elsparam->positionBMparts43(0)*mm,
                       elsparam->positionBMparts43(1)*mm,
                       elsparam->positionBMparts43(2)*mm );
   G4ThreeVector G4PositionOfBMmainSub2
     =  G4ThreeVector( elsparam->positionBMparts44(0)*mm,
         	       elsparam->positionBMparts44(1)*mm,
		       elsparam->positionBMparts44(2)*mm );
   G4ThreeVector G4PositionOfBMmainSub3
     =  G4ThreeVector( elsparam->positionBMparts45(0)*mm,
                       elsparam->positionBMparts45(1)*mm,
                       elsparam->positionBMparts45(2)*mm );

   G4RotationMatrix* G4RotationOfBMmainSub1 = new G4RotationMatrix;
     G4RotationOfBMmainSub1->rotateX( elsparam->rotationBMparts43(0)*deg );
     G4RotationOfBMmainSub1->rotateY( elsparam->rotationBMparts43(1)*deg );
     G4RotationOfBMmainSub1->rotateZ( elsparam->rotationBMparts43(2)*deg );
   G4RotationMatrix* G4RotationOfBMmainSub2 = new G4RotationMatrix;
     G4RotationOfBMmainSub2->rotateX( 0. );
     G4RotationOfBMmainSub2->rotateY( 0. );
     G4RotationOfBMmainSub2->rotateZ( 0. );
   G4RotationMatrix* G4RotationOfBMmainSub3 = new G4RotationMatrix;
     G4RotationOfBMmainSub3->rotateX( elsparam->rotationBMparts45(0)*deg );
     G4RotationOfBMmainSub3->rotateY( elsparam->rotationBMparts45(1)*deg );
     G4RotationOfBMmainSub3->rotateZ( elsparam->rotationBMparts45(2)*deg );

   for( int i=0; i<2; i++ ){
     solidBMmainBody[i] = new G4Box("BMmain",
                                    elsparam->sizeBMparts4m( i, 0 )*0.5*mm,
                                    elsparam->sizeBMparts4m( i, 1 )*0.5*mm,
                                    elsparam->sizeBMparts4m( i, 2 )*0.5*mm );
     G4ThreeVector G4PositionOfBMmain
       = G4ThreeVector( elsparam->positionBMparts4m( i, 0 )*mm,
                        elsparam->positionBMparts4m( i, 1 )*mm,
			elsparam->positionBMparts4m( i, 2 )*mm );
     G4RotationMatrix* G4RotationOfBMmain = new G4RotationMatrix;
     G4RotationOfBMmain->rotateX( 0. );
     G4RotationOfBMmain->rotateY( 0. );
     G4RotationOfBMmain->rotateZ( 0. );

   G4Transform3D Trans3DBMmain( *G4RotationOfBMmain, G4PositionOfBMmain );

     solidBMmain1[i] = new G4SubtractionSolid( "BMmain",
                                            solidBMmainBody[i], 
                                            solidBMmainSub1,
                                            G4RotationOfBMmainSub1,
                                            G4PositionOfBMmainSub1 - G4PositionOfBMmain );
     solidBMmain2[i] = new G4SubtractionSolid( "BMmain",
                                            solidBMmain1[i],
                                            solidBMmainSub2,
                                            G4RotationOfBMmainSub2,
                                            G4PositionOfBMmainSub2 - G4PositionOfBMmain );
     solidBMmain[i] = new G4SubtractionSolid( "BMmain",
                                             solidBMmain2[i],
                                             solidBMmainSub3,
                                             G4RotationOfBMmainSub3,
                                             G4PositionOfBMmainSub3 - G4PositionOfBMmain );

     logicBMmain[i] = new G4LogicalVolume( solidBMmain[i], SUS304, "BMmain", 0, 0, 0, true );
     if( i==0 ) logicBMmain[0]->SetVisAttributes( BMp41VisAtt );
     if( i==1 ) logicBMmain[1]->SetVisAttributes( BMp42VisAtt );
     physiBMmain[i] = new G4PVPlacement( Trans3DBMmain, 
                                         logicBMmain[i],
                                         "BMmain", logicWorld, false, 0 );
   }
  //--------------------------------------------------
  // BM Duct parts 4 (6)-(9)
  //--------------------------------------------------
   G4VisAttributes* BMSubVisAtt = new G4VisAttributes(true,G4Colour(0.3, 0.3, 1.0));  
       BMSubVisAtt->SetForceWireframe(true);
       BMSubVisAtt->SetForceSolid(true);

   solidBMSub1 = new G4Tubs( "BMSub1", 
                              elsparam->rminBMparts46()*mm, elsparam->rmaxBMparts46()*mm,
                              elsparam->lBMparts46()*0.5*mm, 
                              elsparam->phistartBMparts46()*deg,
                              elsparam->phistopBMparts46()*deg ); 
   solidBMSub2 = new G4Box( "BMSub2",
			    elsparam->sizeBMparts47(0)*0.5*mm, 
                            elsparam->sizeBMparts47(1)*0.5*mm,
			    elsparam->sizeBMparts47(2)*0.5*mm );
   solidBMSub3 = new G4Tubs( "BMSub3",
			      elsparam->rminBMparts48()*mm, elsparam->rmaxBMparts48()*mm,
			      elsparam->lBMparts48()*0.5*mm,
			      elsparam->phistartBMparts48()*deg,
			      elsparam->phistopBMparts48()*deg );
   solidBMSub4 = new G4Box( "BMSub4",
			    elsparam->sizeBMparts49(0)*0.5*mm, 
                            elsparam->sizeBMparts49(1)*0.5*mm,
			    elsparam->sizeBMparts49(2)*0.5*mm );

   G4ThreeVector G4PositionOfBMSub1
     = G4ThreeVector( elsparam->positionBMparts46(0)*mm, 
                      elsparam->positionBMparts46(1)*mm,
		      elsparam->positionBMparts46(2)*mm );                               
   G4ThreeVector G4PositionOfBMSub2
     = G4ThreeVector( elsparam->positionBMparts47(0)*mm,
		      elsparam->positionBMparts47(1)*mm,
		      elsparam->positionBMparts47(2)*mm );
   G4ThreeVector G4PositionOfBMSub3
     = G4ThreeVector( elsparam->positionBMparts48(0)*mm,
		      elsparam->positionBMparts48(1)*mm,
		      elsparam->positionBMparts48(2)*mm );
   G4ThreeVector G4PositionOfBMSub4
     = G4ThreeVector( elsparam->positionBMparts49(0)*mm,
	 	      elsparam->positionBMparts49(1)*mm,
		      elsparam->positionBMparts49(2)*mm );
   
   G4RotationMatrix* G4RotationOfBMSub1 = new G4RotationMatrix;
      G4RotationOfBMSub1->rotateX( elsparam->rotationBMparts46(0)*deg );                 
      G4RotationOfBMSub1->rotateY( elsparam->rotationBMparts46(1)*deg );
      G4RotationOfBMSub1->rotateZ( elsparam->rotationBMparts46(2)*deg );

   G4RotationMatrix* G4RotationOfBMSub2 = new G4RotationMatrix;
      G4RotationOfBMSub2->rotateX( 0. ); 
      G4RotationOfBMSub2->rotateY( 0. ); 
      G4RotationOfBMSub2->rotateZ( 0. ); 

   G4RotationMatrix* G4RotationOfBMSub3 = new G4RotationMatrix;
      G4RotationOfBMSub3->rotateX( elsparam->rotationBMparts48(0)*deg );
      G4RotationOfBMSub3->rotateY( elsparam->rotationBMparts48(1)*deg );
      G4RotationOfBMSub3->rotateZ( elsparam->rotationBMparts48(2)*deg );

   G4RotationMatrix* G4RotationOfBMSub4  = new G4RotationMatrix;
      G4RotationOfBMSub4->rotateX( 0. ); 
      G4RotationOfBMSub4->rotateY( 0. ); 
      G4RotationOfBMSub4->rotateZ( 0. ); 

   G4Transform3D Trans3DBMSub1( *G4RotationOfBMSub1, G4PositionOfBMSub1 );
   G4Transform3D Trans3DBMSub2( *G4RotationOfBMSub2, G4PositionOfBMSub2 );
   G4Transform3D Trans3DBMSub3( *G4RotationOfBMSub3, G4PositionOfBMSub3 );
   G4Transform3D Trans3DBMSub4( *G4RotationOfBMSub4 , G4PositionOfBMSub4 );
   
   logicBMSub1 = new G4LogicalVolume( solidBMSub1, SUS304, "BMsub1", 0, 0, 0 );
   logicBMSub2 = new G4LogicalVolume( solidBMSub2, SUS304, "BMsub2", 0, 0, 0 );
   logicBMSub3 = new G4LogicalVolume( solidBMSub3, SUS304, "BMsub3", 0, 0, 0 );
   logicBMSub4 = new G4LogicalVolume( solidBMSub4, SUS304, "BMsub4", 0, 0, 0 );

   logicBMSub1->SetVisAttributes( BMSubVisAtt );
   logicBMSub2->SetVisAttributes( BMSubVisAtt );
   logicBMSub3->SetVisAttributes( BMSubVisAtt );
   logicBMSub4->SetVisAttributes( BMSubVisAtt );

   physiBMSub1 = new G4PVPlacement( Trans3DBMSub1, logicBMSub1,"BMSub1", logicWorld, false, 0 ); 
   physiBMSub2 = new G4PVPlacement( Trans3DBMSub2, logicBMSub2,"BMSub2", logicWorld, false, 0 );
   physiBMSub3 = new G4PVPlacement( Trans3DBMSub3, logicBMSub3,"BMSub3", logicWorld, false, 0 );
   physiBMSub4 = new G4PVPlacement( Trans3DBMSub4, logicBMSub4,"BMSub4", logicWorld, false, 0 );

   //-------------------------------------------------   
   // BM Duct vacuum ( parts5-(1)-(4)
   //-------------------------------------------------
  solidBMmainSub1_vac = new G4Tubs( "BMmainSub1 vac", 
                                  elsparam->rminBMparts52()*mm, elsparam->rmaxBMparts52()*mm,
                                  elsparam->lBMparts52()*0.5*mm, 0.0*deg, 360.0*deg ); 
  solidBMmainSub2_vac = new G4Box( "BMmainSub2 vac", 
			       elsparam->sizeBMparts53(0)*0.5*mm,
                               elsparam->sizeBMparts53(1)*0.5*mm,
                               elsparam->sizeBMparts53(2)*0.5*mm );
  solidBMmainSub3_vac = new G4Tubs( "BMmainSub3_vac",
                                elsparam->rminBMparts54()*mm, elsparam->rmaxBMparts54()*mm,
                                elsparam->lBMparts54()*0.5*mm, 0.0*deg, 360.0*deg );

   G4ThreeVector G4PositionOfBMmainSub1_vac
     =  G4ThreeVector( elsparam->positionBMparts52(0)*mm,
                       elsparam->positionBMparts52(1)*mm,
                       elsparam->positionBMparts52(2)*mm );
   G4ThreeVector G4PositionOfBMmainSub2_vac
     =  G4ThreeVector( elsparam->positionBMparts53(0)*mm,
         	       elsparam->positionBMparts53(1)*mm,
		       elsparam->positionBMparts53(2)*mm );
   G4ThreeVector G4PositionOfBMmainSub3_vac
     =  G4ThreeVector( elsparam->positionBMparts54(0)*mm,
                       elsparam->positionBMparts54(1)*mm,
                       elsparam->positionBMparts54(2)*mm );

   G4RotationMatrix* G4RotationOfBMmainSub1_vac = new G4RotationMatrix;
     G4RotationOfBMmainSub1_vac->rotateX( elsparam->rotationBMparts52(0)*deg );
     G4RotationOfBMmainSub1_vac->rotateY( elsparam->rotationBMparts52(1)*deg );
     G4RotationOfBMmainSub1_vac->rotateZ( elsparam->rotationBMparts52(2)*deg );
   G4RotationMatrix* G4RotationOfBMmainSub2_vac = new G4RotationMatrix;
     G4RotationOfBMmainSub2_vac->rotateX( 0. );
     G4RotationOfBMmainSub2_vac->rotateY( 0. );
     G4RotationOfBMmainSub2_vac->rotateZ( 0. );
   G4RotationMatrix* G4RotationOfBMmainSub3_vac = new G4RotationMatrix;
     G4RotationOfBMmainSub3_vac->rotateX( elsparam->rotationBMparts54(0)*deg );
     G4RotationOfBMmainSub3_vac->rotateY( elsparam->rotationBMparts54(1)*deg );
     G4RotationOfBMmainSub3_vac->rotateZ( elsparam->rotationBMparts54(2)*deg );

     solidBMmainBody_vac = new G4Box("BMmain vac",
                                    elsparam->sizeBMparts51(0)*0.5*mm,
                                    elsparam->sizeBMparts51(1)*0.5*mm,
                                    elsparam->sizeBMparts51(2)*0.5*mm );
     G4ThreeVector G4PositionOfBMmain_vac
       = G4ThreeVector( elsparam->positionBMparts51(0)*mm,
                        elsparam->positionBMparts51(1)*mm,
			elsparam->positionBMparts51(2)*mm );
     G4RotationMatrix* G4RotationOfBMmain_vac = new G4RotationMatrix;
     G4RotationOfBMmain_vac->rotateX( 0. );
     G4RotationOfBMmain_vac->rotateY( 0. );
     G4RotationOfBMmain_vac->rotateZ( 0. );

   G4Transform3D Trans3DBMmain_vac( *G4RotationOfBMmain_vac, G4PositionOfBMmain_vac );

     solidBMmain1_vac = new G4SubtractionSolid( "BMmain vac",
                                            solidBMmainBody_vac, 
                                            solidBMmainSub1_vac,
                                            G4RotationOfBMmainSub1_vac,
                                            G4PositionOfBMmainSub1_vac - G4PositionOfBMmain_vac );
     solidBMmain2_vac = new G4SubtractionSolid( "BMmain vac",
                                            solidBMmain1_vac,
                                            solidBMmainSub2_vac,
                                            G4RotationOfBMmainSub2_vac,
                                            G4PositionOfBMmainSub2_vac - G4PositionOfBMmain_vac );
     solidBMmain_vac = new G4SubtractionSolid( "BMmain vac",
                                             solidBMmain2_vac,
                                             solidBMmainSub3_vac,
                                             G4RotationOfBMmainSub3_vac,
                                             G4PositionOfBMmainSub3_vac - G4PositionOfBMmain_vac );

     logicBMmain_vac = new G4LogicalVolume( solidBMmain_vac, Vacuum, "BMmain vac", 0, 0, 0, true );
     //logicBMmain_vac->SetVisAttributes( BMp41VisAtt );

     physiBMmain_vac = new G4PVPlacement( Trans3DBMmain_vac, 
                                          logicBMmain_vac,
                                          "BMmain vac", logicWorld, false, 0 );

  //--------------------------------------------------
  // BM York & Coil
  //--------------------------------------------------
  // York
     G4VisAttributes* BMYorkVisAtt = new G4VisAttributes(true,G4Colour(0.1, 1.0, 0.1));
     BMYorkVisAtt->SetForceWireframe(true);
     BMYorkVisAtt->SetForceSolid(true);

     for( int i=0; i<5; i++ ){
      solidBMYork[i] = new G4Box( "BM York",
 		                   elsparam->sizeBMYork(i,0)*0.5*mm,
				   elsparam->sizeBMYork(i,1)*0.5*mm,
  				   elsparam->sizeBMYork(i,2)*0.5*mm );
      logicBMYork[i] = new G4LogicalVolume( solidBMYork[i], MaterialYork, 
                                            "BM York", 0, 0, 0 );
      logicBMYork[i]->SetVisAttributes( BMYorkVisAtt );

       G4ThreeVector G4PositionOfBMYork
	 = G4ThreeVector( elsparam->positionBMYork(i,0)*mm,
			  elsparam->positionBMYork(i,1)*mm,
			  elsparam->positionBMYork(i,2)*mm );
       G4RotationMatrix* G4RotationOfBMYork = new G4RotationMatrix;
       G4RotationOfBMYork->rotateX( 0. );
       G4RotationOfBMYork->rotateY( 0. );
       G4RotationOfBMYork->rotateZ( 0. );

       G4Transform3D Trans3DBMYork( *G4RotationOfBMYork, G4PositionOfBMYork );
 
       physiBMYork[i]= new G4PVPlacement( Trans3DBMYork,
					  logicBMYork[i],
					    "BM York", logicWorld, false, 0 );
    }

   // Coil  
     G4VisAttributes* BMCoilVisAtt = new G4VisAttributes(true,G4Colour(1.0, 1.0, 0.1));
     BMCoilVisAtt->SetForceWireframe(true);
     BMCoilVisAtt->SetForceSolid(true);
   
     solidBMCoilOuter  = new G4Box( "BM Coil Outer",
				    elsparam->sizeBMCoilOuter(0)*0.5*mm,
				    elsparam->sizeBMCoilOuter(1)*0.5*mm,
                                    elsparam->sizeBMCoilOuter(2)*0.5*mm );
     solidBMCoilInner  = new G4Box( "BM Coil Inner",
				    elsparam->sizeBMCoilInner(0)*0.5*mm,
				    elsparam->sizeBMCoilInner(1)*0.5*mm,
                                    elsparam->sizeBMCoilInner(2)*0.5*mm );

      solidBMCoil =  new G4SubtractionSolid( "BM Coil",
					     solidBMCoilOuter,
                                             solidBMCoilInner,
                                             G4RotationTemp, G4PositionTemp );
       
      logicBMCoil = new G4LogicalVolume( solidBMCoil,  MaterialCoil,
					      "BM Coil", 0, 0, 0 );
      logicBMCoil->SetVisAttributes( BMCoilVisAtt );

      for( int i=0; i<2; i++ ){

       G4ThreeVector G4PositionOfBMCoil
         = G4ThreeVector( elsparam->positionBMCoil(i,0)*mm,
                          elsparam->positionBMCoil(i,1)*mm,
                          elsparam->positionBMCoil(i,2)*mm );
       G4RotationMatrix* G4RotationOfBMCoil = new G4RotationMatrix;
       G4RotationOfBMCoil->rotateX( 0. );
       G4RotationOfBMCoil->rotateY( 0. );
       G4RotationOfBMCoil->rotateZ( 0. );

       G4Transform3D Trans3DBMCoil( *G4RotationOfBMCoil, G4PositionOfBMCoil );
   
      physiBMCoil[i]= new G4PVPlacement( Trans3DBMCoil,
				         logicBMCoil,
				         "BM Coil", logicWorld, false, 0 );
      }

  //--------------------------------------------------
  // SLIT
  //--------------------------------------------------
    G4VisAttributes* SLITFVisAtt = new G4VisAttributes(true,G4Colour(0.3,1.0,0.3));  
      SLITFVisAtt->SetForceWireframe(true);
      SLITFVisAtt->SetForceSolid(true);

    G4VisAttributes* SLITBVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,0.2));  
      SLITBVisAtt->SetForceWireframe(true);
      //SLITBVisAtt->SetForceSolid(true);

    G4VisAttributes* SLITDVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,0.4));  
      SLITDVisAtt->SetForceWireframe(true);
      SLITDVisAtt->SetForceSolid(true);

    G4VisAttributes* COLLIVisAtt = new G4VisAttributes(true,G4Colour(0.0,1.0,1.0));
      COLLIVisAtt->SetForceWireframe(true);
      COLLIVisAtt->SetForceSolid(true);

     // Frange 
     solidSLIT_FF_tubs = new G4Tubs( "SLIT FF",
                                   elsparam->rminSLITff()*mm, elsparam->rmaxSLITff()*mm,
                                   elsparam->lSLITff()*0.5*mm, 0.*deg, 360.0*deg );
     solidSLIT_FF_vac_tubs = new G4Tubs( "SLIT FF Vac",
                                   0.0*mm, elsparam->rminSLITff()*mm, 
                                   elsparam->lSLITff()*0.5*mm, 0.*deg, 360.0*deg );

     solidSLIT_BF_tubs = new G4Tubs( "SLIT BF",
                                   elsparam->rminSLITbf()*mm, elsparam->rmaxSLITbf()*mm,
                                   elsparam->lSLITbf()*0.5*mm, 0.*deg, 360.0*deg );     
     solidSLIT_BF_vac_tubs = new G4Tubs( "SLIT BF Vac",
                                   0.0*mm, elsparam->rminSLITbf()*mm, 
                                   elsparam->lSLITbf()*0.5*mm, 0.*deg, 360.0*deg );

     G4ThreeVector G4PositionOfSLITff
	 = G4ThreeVector( elsparam->positionSLITff(0)*mm,
			  elsparam->positionSLITff(1)*mm,
			  elsparam->positionSLITff(2)*mm );
     G4ThreeVector G4PositionOfSLITbf 
	 = G4ThreeVector( elsparam->positionSLITbf(0)*mm,
		 	  elsparam->positionSLITbf(1)*mm,
			  elsparam->positionSLITbf(2)*mm );				

     G4RotationMatrix* G4RotationOfSLITf = new G4RotationMatrix;     
	 G4RotationOfSLITf->rotateX( elsparam->rotationSLITf(0)*deg );	 
	 G4RotationOfSLITf->rotateY( elsparam->rotationSLITf(1)*deg );	 
	 G4RotationOfSLITf->rotateZ( elsparam->rotationSLITf(2)*deg );	 

      G4Transform3D Trans3DSLITFF( *G4RotationOfSLITf, G4PositionOfSLITff );
      G4Transform3D Trans3DSLITBF( *G4RotationOfSLITf, G4PositionOfSLITbf );

    logicSLIT_FF_tubs = new G4LogicalVolume( solidSLIT_FF_tubs, SUS304, "SLIT FF", 0, 0, 0 );
    logicSLIT_FF_tubs->SetVisAttributes( SLITFVisAtt );
    logicSLIT_FF_vac_tubs = new G4LogicalVolume( solidSLIT_FF_vac_tubs, Vacuum, "SLIT FF Vac", 0, 0, 0 );
    //logicSLIT_FF_vac_tubs->SetVisAttributes( SLITFVisAtt );

    logicSLIT_BF_tubs = new G4LogicalVolume( solidSLIT_BF_tubs, SUS304, "SLIT BF", 0, 0, 0 );
    logicSLIT_BF_tubs->SetVisAttributes( SLITFVisAtt );
    logicSLIT_BF_vac_tubs = new G4LogicalVolume( solidSLIT_BF_vac_tubs, Vacuum, "SLIT BF Vac", 0, 0, 0 );
    //logicSLIT_BF_vac_tubs->SetVisAttributes( SLITFVisAtt );

    physiSLIT_FF_tubs = new G4PVPlacement( Trans3DSLITFF, logicSLIT_FF_tubs, 
                                            "SLIT FF", logicWorld, false, 0 );
    physiSLIT_FF_vac_tubs = new G4PVPlacement( Trans3DSLITFF, logicSLIT_FF_vac_tubs,
					   "SLIT FF Vac", logicWorld, false, 0 );

    physiSLIT_BF_tubs = new G4PVPlacement( Trans3DSLITBF, logicSLIT_BF_tubs, 
                                            "SLIT BF", logicWorld, false, 0 );
    physiSLIT_BF_vac_tubs = new G4PVPlacement( Trans3DSLITBF, logicSLIT_BF_vac_tubs,
					   "SLIT BF Vac", logicWorld, false, 0 );

    // Body & Duct
    solidSLIT_parts1_tubs =  new G4Tubs( "SLIT parts1",
                                          elsparam->rminSLITparts1()*mm, elsparam->rmaxSLITparts1()*mm,
                                          elsparam->lSLITparts1()*0.5*mm, 0.*deg, 360.0*deg );
    solidSLIT_parts12_tubs =  new G4Tubs( "SLIT parts12",
                                          0.0*mm, elsparam->rmaxSLITparts1()*mm,
                                          elsparam->lSLITparts1()*0.5*mm, 0.*deg, 360.0*deg );
    solidSLIT_parts13_tubs =  new G4Tubs( "SLIT parts13",
                                          0.0*mm, elsparam->rminSLITparts1()*mm,
                                          elsparam->lSLITparts1()*0.5*mm, 0.*deg, 360.0*deg );
    solidSLIT_parts2_tubs =  new G4Tubs( "SLIT parts2",
                                          elsparam->rminSLITparts2()*mm, elsparam->rmaxSLITparts2()*mm,
                                          elsparam->lSLITparts2()*0.5*mm, 0.*deg, 360.0*deg );
    solidSLIT_parts3_tubs =  new G4Tubs( "SLIT parts3",
                                          elsparam->rminSLITparts3()*mm, elsparam->rmaxSLITparts3()*mm,
                                          elsparam->lSLITparts3()*0.5*mm, 0.*deg, 360.0*deg );

    G4ThreeVector G4PositionOfSLITparts1
      = G4ThreeVector( elsparam->positionSLITparts1(0)*mm,
		       elsparam->positionSLITparts1(1)*mm,
		       elsparam->positionSLITparts1(2)*mm );
    G4ThreeVector G4PositionOfSLITparts2
      = G4ThreeVector( elsparam->positionSLITparts2(0)*mm,
		       elsparam->positionSLITparts2(1)*mm,
		       elsparam->positionSLITparts2(2)*mm );
    G4ThreeVector G4PositionOfSLITparts3
      = G4ThreeVector( elsparam->positionSLITparts3(0)*mm,
                       elsparam->positionSLITparts3(1)*mm,
                       elsparam->positionSLITparts3(2)*mm );
   
    G4RotationMatrix* G4RotationOfSLITcomp = new G4RotationMatrix;
      G4RotationOfSLITcomp->rotateX( elsparam->rotationSLITcomp(0)*deg );
      G4RotationOfSLITcomp->rotateY( elsparam->rotationSLITcomp(1)*deg );
      G4RotationOfSLITcomp->rotateZ( elsparam->rotationSLITcomp(2)*deg );

    G4RotationMatrix* G4RotationOfSLITduct = new G4RotationMatrix;
      G4RotationOfSLITduct->rotateX( elsparam->rotationSLITduct(0)*deg );
      G4RotationOfSLITduct->rotateY( elsparam->rotationSLITduct(1)*deg );
      G4RotationOfSLITduct->rotateZ( elsparam->rotationSLITduct(2)*deg );

    G4RotationMatrix* G4RotationOfSLITbody = new G4RotationMatrix;
      G4RotationOfSLITbody->rotateX( elsparam->rotationSLITbody(0)*deg );
      G4RotationOfSLITbody->rotateY( elsparam->rotationSLITbody(1)*deg );
      G4RotationOfSLITbody->rotateZ( elsparam->rotationSLITbody(2)*deg );

    G4Transform3D Trans3DSLITDUCT( *G4RotationOfSLITduct, G4PositionOfSLITparts1 );    
    G4Transform3D Trans3DSLITBODY( *G4RotationOfSLITbody, G4PositionOfSLITparts3 );

    solidSLIT_duct = new G4SubtractionSolid( "SLIT Duct",
					      solidSLIT_parts1_tubs,
					      solidSLIT_parts2_tubs,
					      G4RotationOfSLITcomp,
					      G4PositionOfSLITparts2 - G4PositionOfSLITparts1 );
    solidSLIT_duct_vac = new G4SubtractionSolid( "SLIT Duct Vacuum",
						 solidSLIT_parts13_tubs,
						 solidSLIT_parts2_tubs,
						 G4RotationOfSLITcomp,
						 G4PositionOfSLITparts2 - G4PositionOfSLITparts1 );

    solidSLIT_body = new G4SubtractionSolid( "SLIT Body",
                                              solidSLIT_parts3_tubs,
                                              solidSLIT_parts12_tubs,
                                              G4RotationOfSLITcomp,
                                              G4PositionOfSLITparts1 - G4PositionOfSLITparts3 );
    
    logicSLIT_duct = new G4LogicalVolume( solidSLIT_duct,
                                           SUS304, "SLIT Duct", 0, 0, 0 );
    logicSLIT_duct->SetVisAttributes( SLITDVisAtt );
    logicSLIT_duct_vac = new G4LogicalVolume( solidSLIT_duct_vac,
	 				       Vacuum, "SLIT Duct Vac", 0, 0, 0 );
    //logicSLIT_duct_vac->SetVisAttributes( SLITDVisAtt );

    physiSLIT_duct = new G4PVPlacement( Trans3DSLITDUCT, logicSLIT_duct,
                                         "SLIT Duct", logicWorld, false, 0 );
    physiSLIT_duct_vac = new G4PVPlacement(  Trans3DSLITDUCT,logicSLIT_duct_vac,
					      "SLIT Duct Vac", logicWorld, false, 0 );     

    logicSLIT_body = new G4LogicalVolume( solidSLIT_body,
                                           SUS304, "SLIT Body", 0, 0, 0 );
    logicSLIT_body->SetVisAttributes( SLITBVisAtt );  

    physiSLIT_body = new G4PVPlacement( Trans3DSLITBODY, logicSLIT_body,
                                         "SLIT Body", logicWorld, false, 0 );

    // Collimeter

    //SLIT BODY VAC 
    solidCollimeter_vac = new G4Box("Collimeter Vac",
                             	   elsparam->sizeCollimeter_Vac(0)*0.5*mm,
	  			   elsparam->sizeCollimeter_Vac(1)*0.5*mm,
				   elsparam->sizeCollimeter_Vac(2)*0.5*mm );
      G4ThreeVector G4PositionOfCollimeter_vac_0
 	= G4ThreeVector( elsparam->positionCollimeter_Vac(0,0)*mm,
 	                 elsparam->positionCollimeter_Vac(0,1)*mm,
                         elsparam->positionCollimeter_Vac(0,2)*mm + Collimeter_Width/2.0 );
      G4ThreeVector G4PositionOfCollimeter_vac_1
  	= G4ThreeVector( elsparam->positionCollimeter_Vac(1,0)*mm,
 	                 elsparam->positionCollimeter_Vac(1,1)*mm,
                         elsparam->positionCollimeter_Vac(1,2)*mm - Collimeter_Width/2.0 );

    solidSLIT_body_vac_0 = new G4Tubs( "SLIT Body Vac 0",
                                      elsparam->rminSLITparts2()*mm, elsparam->rmaxSLITparts2()*mm,
                                      elsparam->lSLITparts2()*0.5*mm, 0.*deg, 360.0*deg );

      solidSLIT_body_vac_1 = new G4SubtractionSolid( "SLIT BODY Vac 1",
					                           solidSLIT_parts2_tubs,
					                           solidCollimeter_vac,
                                               0,
                                               G4PositionOfCollimeter_vac_0 - G4PositionOfSLITparts3 );
       
      solidSLIT_body_vac = new G4SubtractionSolid( "SLIT BODY Vac",
                                                   solidSLIT_body_vac_1,
                                                   solidCollimeter_vac,
						   0,
						   G4PositionOfCollimeter_vac_1 - G4PositionOfSLITparts3 );

      logicSLIT_body_vac  =  new G4LogicalVolume( solidSLIT_body_vac, 
				 		  Vacuum, "SLIT Body Vac", 0, 0, 0 );
      //logicSLIT_body_vac->SetVisAttributes( SLITBVisAtt );  
      physiSLIT_body_vac = new G4PVPlacement( Trans3DSLITBODY, logicSLIT_body_vac,
       	  				      "SLIT Body Vac", logicWorld, false, 0 );

    G4double cw(1.);
    for( int i=0; i < 2; i++ ){

    solidCollimeter[i] = new G4Box("Collimeter",
	      		  	   elsparam->sizeCollimeter(0)*0.5*mm,
	  			   elsparam->sizeCollimeter(1)*0.5*mm,
				   elsparam->sizeCollimeter(2)*0.5*mm );    
      logicCollimeter[i] = new G4LogicalVolume( solidCollimeter[i], MaterialCollimeter,
                                             "Collimeter",  0, 0, 0 );
      logicCollimeter[i]->SetVisAttributes( COLLIVisAtt );

      if( i==0 ) cw=1.;
      if( i==1 ) cw=-1.;

      cout << " Collimeter_Width = " << Collimeter_Width << G4endl;

      G4ThreeVector G4PositionOfCollimeter
          = G4ThreeVector( elsparam->positionCollimeter(i,0)*mm + cw*Collimeter_Width/2.0,
	   	           elsparam->positionCollimeter(i,1)*mm,
		           elsparam->positionCollimeter(i,2)*mm );

      G4RotationMatrix* G4RotationOfCollimeter = new G4RotationMatrix;
        G4RotationOfCollimeter->rotateX( 0.0*deg );
        G4RotationOfCollimeter->rotateY( 0.0*deg );
        G4RotationOfCollimeter->rotateZ( 0.0*deg );
      G4Transform3D Trans3DCollimeter( *G4RotationOfCollimeter, G4PositionOfCollimeter );
      physiCollimeter[i] = new G4PVPlacement( Trans3DCollimeter, logicCollimeter[i] ,
      				      "Collimeter", logicWorld, false, 0 );
      }

    //--------------------------------------------------
    // Blank Frange
    //--------------------------------------------------
     G4VisAttributes* BLANKVisAtt = new G4VisAttributes(true,G4Colour(1.0,0.2,0.2));  
      BLANKVisAtt->SetForceWireframe(true);
      BLANKVisAtt->SetForceSolid(true);

      solidBlankFrange = new G4Tubs("Blank",
                                    elsparam->rminBlank()*mm, elsparam->rmaxBlank()*mm, 
                                    elsparam->lBlank()*0.5*mm, 0.*deg, 360.*deg ); 
      logicBlankFrange = new G4LogicalVolume( solidBlankFrange, SUS304, "Blank",  0, 0, 0 );
      logicBlankFrange->SetVisAttributes(BLANKVisAtt);
       
      G4ThreeVector G4PositionOfBlank
        = G4ThreeVector( elsparam->positionBlank(0)*mm,
                         elsparam->positionBlank(1)*mm,
			 elsparam->positionBlank(2)*mm );
      G4RotationMatrix* G4RotationOfBlank = new G4RotationMatrix;
          G4RotationOfBlank->rotateX( elsparam->rotationBlank(0)*deg );
          G4RotationOfBlank->rotateY( elsparam->rotationBlank(1)*deg );
          G4RotationOfBlank->rotateZ( elsparam->rotationBlank(2)*deg );
      G4Transform3D Trans3DBlank( *G4RotationOfBlank, G4PositionOfBlank );
      
      physiBlankFrange = new G4PVPlacement( Trans3DBlank, logicBlankFrange ,
					    "Blank", logicWorld, false, 0 );

    //-- ------------------------------------------------
    // Output Frange & Window
    //--------------------------------------------------

     G4VisAttributes* WINDOWVisAtt = new G4VisAttributes(true,G4Colour(0.6,0.6,0.6));  
      WINDOWVisAtt->SetForceWireframe(true);
      WINDOWVisAtt->SetForceSolid(true);

      solidWindow = new G4Tubs("Window",
                                    elsparam->rminWindow()*mm, elsparam->rmaxWindow()*mm, 
                                    elsparam->lWindow()*0.5*mm, 0.*deg, 360.*deg ); 
      logicWindow = new G4LogicalVolume( solidWindow, SUS304, "Window",  0, 0, 0 );
      logicWindow->SetVisAttributes(WINDOWVisAtt);

      solidWindow_vac = new G4Tubs("Window Vac",
                                    0.0*mm, elsparam->rminWindow()*mm,
                                    elsparam->lWindow_vac()*0.5*mm, 
                                    0.*deg, 360.*deg ); 
      logicWindow_vac = new G4LogicalVolume( solidWindow_vac, Vacuum, "Window Vac",  0, 0, 0 );

      G4ThreeVector G4PositionOfWindow
        = G4ThreeVector( elsparam->positionWindow(0)*mm,
                         elsparam->positionWindow(1)*mm,
			 elsparam->positionWindow(2)*mm );
      G4ThreeVector G4PositionOfWindow_vac
        = G4ThreeVector( elsparam->positionWindow(0)*mm,
                         elsparam->positionWindow(1)*mm,
			 elsparam->position_z_Window_vac()*mm );
      G4RotationMatrix* G4RotationOfWindow = new G4RotationMatrix;
          G4RotationOfWindow->rotateX( elsparam->rotationWindow(0)*deg );
          G4RotationOfWindow->rotateY( elsparam->rotationWindow(1)*deg );
          G4RotationOfWindow->rotateZ( elsparam->rotationWindow(2)*deg );
      G4Transform3D Trans3DWindow( *G4RotationOfWindow, G4PositionOfWindow );
      G4Transform3D Trans3DWindow_vac( *G4RotationOfWindow, G4PositionOfWindow_vac );
      
      physiWindow = new G4PVPlacement( Trans3DWindow, logicWindow ,
					    "Window", logicWorld, false, 0 );
      physiWindow_vac = new G4PVPlacement( Trans3DWindow_vac, logicWindow_vac ,
					    "Window Vac", logicWorld, false, 0 );

      // TiWindow
     G4VisAttributes* TIWINDOWVisAtt = new G4VisAttributes(true,G4Colour(0.3,0.3,0.3));  
      TIWINDOWVisAtt->SetForceWireframe(true);
      TIWINDOWVisAtt->SetForceSolid(true);

      solidTiWindow = new G4Tubs("TiWindow",
                                    elsparam->rminTiWindow()*mm, elsparam->rmaxTiWindow()*mm, 
                                    elsparam->lTiWindow()*0.5*mm, 0.*deg, 360.*deg ); 

      logicTiWindow = new G4LogicalVolume( solidTiWindow, MaterialtTiwindow,
      //  Changed(6/7) for FC study 2011.04.14      
      //logicTiWindow = new G4LogicalVolume( solidTiWindow, SUS304,
					     "TiWindow",  0, 0, 0 );
      logicTiWindow->SetVisAttributes(TIWINDOWVisAtt);

       
      G4ThreeVector G4PositionOfTiWindow
        = G4ThreeVector( elsparam->positionTiWindow(0)*mm,
                         elsparam->positionTiWindow(1)*mm,
			 elsparam->positionTiWindow(2)*mm );
      G4RotationMatrix* G4RotationOfTiWindow = new G4RotationMatrix;
          G4RotationOfTiWindow->rotateX( elsparam->rotationTiWindow(0)*deg );
          G4RotationOfTiWindow->rotateY( elsparam->rotationTiWindow(1)*deg );
          G4RotationOfTiWindow->rotateZ( elsparam->rotationTiWindow(2)*deg );
      G4Transform3D Trans3DTiWindow( *G4RotationOfTiWindow, G4PositionOfTiWindow );
      
      physiTiWindow = new G4PVPlacement( Trans3DTiWindow, logicTiWindow ,
					    "TiWindow", logicWorld, false, 0 );
  

    //-- ------------------------------------------------
    // Top Plate
    //---------------------------------------------------
     G4VisAttributes* TopPlateVisAtt = new G4VisAttributes(true,G4Colour(0.6,0.6,0.6));  
      TopPlateVisAtt->SetForceWireframe(true);
      TopPlateVisAtt->SetForceSolid(true);
    
      solidTopPlate_body = new G4Box("TopPlate",
                                      elsparam->sizeTopPlate(0)*0.5*mm, 
                                      elsparam->sizeTopPlate(1)*0.5*mm, 
                                      elsparam->sizeTopPlate(2)*0.5*mm );

      G4ThreeVector G4PositionOfTopPlate
                    = G4ThreeVector( elsparam->positionTopPlate(0)*mm,
                                     elsparam->positionTopPlate(1)*mm,
			             elsparam->positionTopPlate(2)*mm );      

      G4RotationMatrix* G4RotationOfTopPlate = new G4RotationMatrix;
      G4RotationOfTopPlate->rotateX( elsparam->rotationTopPlate(0)*deg );
      G4RotationOfTopPlate->rotateY( elsparam->rotationTopPlate(1)*deg );
      G4RotationOfTopPlate->rotateZ( elsparam->rotationTopPlate(2)*deg );

      G4Transform3D Trans3DTopPlate( *G4RotationOfTopPlate, G4PositionOfTopPlate );

      solidTopPlate_hole  = new G4Tubs("TopPlateHole",
				       elsparam->rmin_TopPlate_hole()*mm,
				       elsparam->rmax_TopPlate_hole()*mm,
				       elsparam->l_TopPlate_hole()*0.5*mm,
				       0.*deg, 360.*deg );
      G4ThreeVector G4PositionOfTopPlateHole
	  = G4ThreeVector( elsparam->positionTopPlate_hole(0)*mm,
			   elsparam->positionTopPlate_hole(1)*mm,
			   elsparam->positionTopPlate_hole(2)*mm );
      G4RotationMatrix* G4RotationOfTopPlateHole = new G4RotationMatrix;
      G4RotationOfTopPlateHole->rotateX( elsparam->rotationTopPlate_hole(0)*deg );
      G4RotationOfTopPlateHole->rotateY( elsparam->rotationTopPlate_hole(1)*deg );
      G4RotationOfTopPlateHole->rotateZ( elsparam->rotationTopPlate_hole(2)*deg );

      solidTopPlate = new G4SubtractionSolid( "TopPlate",
					      solidTopPlate_body,
					      solidTopPlate_hole,
                                              0,
                                              G4PositionOfTopPlateHole - G4PositionOfTopPlate );  

      logicTopPlate = new G4LogicalVolume( solidTopPlate,
                                           SUS304,
					   "TopPlate", 0, 0, 0 );

      // Comment out by T.Shibata in 2012.10.02, bacause we already removed this plate in Mar.2012.
      /*  
      physiTopPlate = new G4PVPlacement( Trans3DTopPlate, logicTopPlate,
					  "TopPlate", logicWorld, false, 0 );
      */

    //--------------------------------------------------
    // Faraday Cup
    //--------------------------------------------------
    if( Faradaycup_flag == 1 && Faradaycup4_flag == 0 && ScreenMonitor3_flag == 0){      

     // Al Cylinder
     G4VisAttributes* FDAlCylinderVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,1.0));  
      FDAlCylinderVisAtt->SetForceWireframe(true);
      FDAlCylinderVisAtt->SetForceSolid(true);

      for( int i=0; i<3; i++ ){

      solidFDAlCylinder[i] = new G4Tubs("FDAlCylinder",
                                        elsparam->rminFDAlCylinder(i)*mm, 
                                        elsparam->rmaxFDAlCylinder(i)*mm,
                                        elsparam->lFDAlCylinder(i)*0.5*mm,
                                        0.*deg, 360.*deg ); 
      logicFDAlCylinder[i]= new G4LogicalVolume( solidFDAlCylinder[i], 
                                                 MaterialtFDAlCylinder, 
                                                 "FDAlCylinder", 0, 0, 0 );
      logicFDAlCylinder[i]->SetVisAttributes( FDAlCylinderVisAtt );
       
      G4ThreeVector G4PositionOfFDAlCylinder
        = G4ThreeVector( elsparam->positionFDAlCylinder(i,0)*mm,
                         elsparam->positionFDAlCylinder(i,1)*mm,
			 elsparam->positionFDAlCylinder(i,2)*mm );
      G4RotationMatrix* G4RotationOfFDAlCylinder = new G4RotationMatrix;
          G4RotationOfFDAlCylinder->rotateX( elsparam->rotationFDAlCylinder(i,0)*deg );
          G4RotationOfFDAlCylinder->rotateY( elsparam->rotationFDAlCylinder(i,1)*deg );
          G4RotationOfFDAlCylinder->rotateZ( elsparam->rotationFDAlCylinder(i,2)*deg );
      G4Transform3D Trans3DFDAlCylinder( *G4RotationOfFDAlCylinder, G4PositionOfFDAlCylinder ) ;

     
      physiFDAlCylinder[i] = new G4PVPlacement( Trans3DFDAlCylinder, logicFDAlCylinder[i],
	   			 		    "FDALCylinder", logicWorld, false, 0 );
        // Changed(7/7) for FC study 2011.04.14           
      }

     // FC Body
     G4VisAttributes* FDBodyVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));  
      FDBodyVisAtt->SetForceWireframe(true);
      FDBodyVisAtt->SetForceSolid(true);

      for( int i=0; i<4; i++ ){

      solidFDBody[i] = new G4Tubs("FCBody",
                                        elsparam->rminFDBody(i)*mm, 
                                        elsparam->rmaxFDBody(i)*mm,
                                        elsparam->lFDBody(i)*0.5*mm,
                                        0.*deg, 360.*deg ); 
      if( i!= 3 ){
      logicFDBody[i]= new G4LogicalVolume( solidFDBody[i], 
					   MaterialtFDPbBody,
                                             "FCBody", 0, 0, 0 );
      }else if( i == 3 ){
      logicFDBody[i]= new G4LogicalVolume( solidFDBody[i],
					   MaterialtFDCBody,
					   "FCBody", 0, 0, 0 );
      }

      logicFDBody[i]->SetVisAttributes( FDBodyVisAtt );
       
      G4ThreeVector G4PositionOfFDBody
        = G4ThreeVector( elsparam->positionFDBody(i,0)*mm,
                         elsparam->positionFDBody(i,1)*mm,
			 elsparam->positionFDBody(i,2)*mm );

      G4cout << "FCBody " << G4PositionOfFDBody << endl;

      G4RotationMatrix* G4RotationOfFDBody = new G4RotationMatrix;
          G4RotationOfFDBody->rotateX( elsparam->rotationFDBody(i,0)*deg );
          G4RotationOfFDBody->rotateY( elsparam->rotationFDBody(i,1)*deg );
          G4RotationOfFDBody->rotateZ( elsparam->rotationFDBody(i,2)*deg );
      G4Transform3D Trans3DFDBody( *G4RotationOfFDBody, G4PositionOfFDBody ) ;

      physiFDBody[i] = new G4PVPlacement( Trans3DFDBody, logicFDBody[i],
    	           	                  "FCBody", logicWorld, false, 0 );

      }

    }

   //--------------------------------------------------
   // Faraday Cup 4
   // FC2-> FC4 in 2012.05.01                                                                                    
   // Front : Carbon                                                                                             
   //       --> Front1/2 thin shield ( 1 and 2 ) Cu or Al or Any other material ?  in 2012.05.01         
   // Front(Carbon)/Body(Copper)
   // Creat Faraday Cup 4 in 2012.10.02                
   // Geometry was fixed except position                                  
   //----- ---------------------------------------------
    if( ( Faradaycup_flag == 4 || Faradaycup4_flag == 1 ) && ScreenMonitor3_flag == 0 ){

   // Ground Layer (G4Polycone)                

      const G4int numZPlanesGNDLayerFC4 = elsparam->numZPlanesGNDLayerFC4();
      for( int i=0; i<numZPlanesGNDLayerFC4; i++ ){
	irGNDLayerFC4[i]=elsparam->irGNDLayerFC4(i);
        orGNDLayerFC4[i]=elsparam->orGNDLayerFC4(i);
        zGNDLayerFC4[i]=elsparam->zGNDLayerFC4(i);
      }

     G4VisAttributes* FC4GNDLayerVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4GNDLayerVisAtt->SetForceWireframe(true);
     FC4GNDLayerVisAtt->SetForceSolid(true);

     solidGNDLayerFC4 = new G4Polycone("FC4GNDLayer",
				        0.0*deg, 360.0*deg, 
				        numZPlanesGNDLayerFC4, zGNDLayerFC4, irGNDLayerFC4, orGNDLayerFC4 );

      logicGNDLayerFC4 = new G4LogicalVolume( solidGNDLayerFC4, MaterialToughPitchCopper, 
                                              "FC4GNDLayer", 0, 0, 0 );
      logicGNDLayerFC4->SetVisAttributes(FC4GNDLayerVisAtt);

      G4ThreeVector G4PositionOfGNDLayerFC4
	   = G4ThreeVector( elsparam->positionGNDLayerFC4(0)*mm,
			    elsparam->positionGNDLayerFC4(1)*mm,
			    elsparam->positionGNDLayerFC4(2)*mm );

      G4RotationMatrix* G4RotationOfGNDLayerFC4 = new G4RotationMatrix;
	 G4RotationOfGNDLayerFC4->rotateX( elsparam->rotationGNDLayerFC4(0)*deg );
	 G4RotationOfGNDLayerFC4->rotateY( elsparam->rotationGNDLayerFC4(1)*deg );
	 G4RotationOfGNDLayerFC4->rotateZ( elsparam->rotationGNDLayerFC4(2)*deg );
	 G4Transform3D Trans3DGNDLayerFC4( *G4RotationOfGNDLayerFC4, G4PositionOfGNDLayerFC4 ) ;

         physiGNDLayerFC4 = new G4PVPlacement( Trans3DGNDLayerFC4, logicGNDLayerFC4,
       			         	      "FC4GNDLayer", logicWorld, false, 0 );

   // Shield Layer (G4Polycone)    
      numZPlanesShieldLayerFC4=elsparam->numZPlanesShieldLayerFC4();
      for( int i=0; i<numZPlanesShieldLayerFC4; i++ ){
	   irShieldLayerFC4[i]=elsparam->irShieldLayerFC4(i);
           orShieldLayerFC4[i]=elsparam->orShieldLayerFC4(i);
           zShieldLayerFC4[i]=elsparam->zShieldLayerFC4(i);
      }

     G4VisAttributes* FC4ShieldLayerVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4ShieldLayerVisAtt->SetForceWireframe(true);
     FC4ShieldLayerVisAtt->SetForceSolid(true);

     solidShieldLayerFC4 = new G4Polycone("FC4ShieldLayer",
				       0.0, 360.0*deg, 
				       numZPlanesShieldLayerFC4, zShieldLayerFC4, irShieldLayerFC4, orShieldLayerFC4 );

      logicShieldLayerFC4 = new G4LogicalVolume( solidShieldLayerFC4, MaterialToughPitchCopper, 
                                              "FC4ShieldLayer", 0, 0, 0 );
      logicShieldLayerFC4->SetVisAttributes(FC4ShieldLayerVisAtt);

      G4ThreeVector G4PositionOfShieldLayerFC4
	   = G4ThreeVector( elsparam->positionShieldLayerFC4(0)*mm,
			    elsparam->positionShieldLayerFC4(1)*mm,
			    elsparam->positionShieldLayerFC4(2)*mm );

      G4RotationMatrix* G4RotationOfShieldLayerFC4 = new G4RotationMatrix;
	 G4RotationOfShieldLayerFC4->rotateX( elsparam->rotationShieldLayerFC4(0)*deg );
	 G4RotationOfShieldLayerFC4->rotateY( elsparam->rotationShieldLayerFC4(1)*deg );
	 G4RotationOfShieldLayerFC4->rotateZ( elsparam->rotationShieldLayerFC4(2)*deg );
	 G4Transform3D Trans3DShieldLayerFC4( *G4RotationOfShieldLayerFC4, G4PositionOfShieldLayerFC4 ) ;

         physiShieldLayerFC4 = new G4PVPlacement( Trans3DShieldLayerFC4, logicShieldLayerFC4,
        			         	          "FC4ShieldLayer", logicWorld, false, 0 );

   // Body : Copper (G4Polycone)
      numZPlanesBodyFC4=elsparam->numZPlanesBodyFC4();
      for( int i=0; i<numZPlanesBodyFC4; i++ ){
   	   irBodyFC4[i]=elsparam->irBodyFC4(i);
           orBodyFC4[i]=elsparam->orBodyFC4(i);
           zBodyFC4[i]=elsparam->zBodyFC4(i);
      }

     G4VisAttributes* FC4BodyVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4BodyVisAtt->SetForceWireframe(true);
     FC4BodyVisAtt->SetForceSolid(true);

     solidBodyFC4 = new G4Polycone("FC4Body",
				       0.0, 360.0*deg, 
				       numZPlanesBodyFC4, zBodyFC4, irBodyFC4, orBodyFC4 );

      logicBodyFC4 = new G4LogicalVolume( solidBodyFC4, MaterialToughPitchCopper, 
                                              "FC4Body", 0, 0, 0 );
      logicBodyFC4->SetVisAttributes(FC4BodyVisAtt);

      G4ThreeVector G4PositionOfBodyFC4
	   = G4ThreeVector( elsparam->positionBodyFC4(0)*mm,
			    elsparam->positionBodyFC4(1)*mm,
			    elsparam->positionBodyFC4(2)*mm );

      G4RotationMatrix* G4RotationOfBodyFC4 = new G4RotationMatrix;
	 G4RotationOfBodyFC4->rotateX( elsparam->rotationBodyFC4(0)*deg );
	 G4RotationOfBodyFC4->rotateY( elsparam->rotationBodyFC4(1)*deg );
	 G4RotationOfBodyFC4->rotateZ( elsparam->rotationBodyFC4(2)*deg );
	 G4Transform3D Trans3DBodyFC4( *G4RotationOfBodyFC4, G4PositionOfBodyFC4 ) ;

         physiBodyFC4 = new G4PVPlacement( Trans3DBodyFC4, logicBodyFC4,
      			        	   "FC4Body", logicWorld, false, 0 );

   // Isolator-1 : Macor                
     G4VisAttributes* FC4Iso1VisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4Iso1VisAtt->SetForceWireframe(true);
     FC4Iso1VisAtt->SetForceSolid(true);

     solidIso1FC4 = new G4Tubs("FC4Iso1",
			        elsparam->rminIso1FC4()*mm,
			        elsparam->rmaxIso1FC4()*mm,
				elsparam->lIso1FC4()*0.5*mm,
				0.*deg, 360.*deg );

      logicIso1FC4 = new G4LogicalVolume( solidIso1FC4, MaterialMacor, 
                                              "FC4Iso1", 0, 0, 0 );
      logicIso1FC4->SetVisAttributes(FC4Iso1VisAtt);

      G4ThreeVector G4PositionOfIso1FC4
	   = G4ThreeVector( elsparam->positionIso1FC4(0)*mm,
			    elsparam->positionIso1FC4(1)*mm,
			    elsparam->positionIso1FC4(2)*mm );

      G4RotationMatrix* G4RotationOfIso1FC4 = new G4RotationMatrix;
	 G4RotationOfIso1FC4->rotateX( elsparam->rotationIso1FC4(0)*deg );
	 G4RotationOfIso1FC4->rotateY( elsparam->rotationIso1FC4(1)*deg );
	 G4RotationOfIso1FC4->rotateZ( elsparam->rotationIso1FC4(2)*deg );
	 G4Transform3D Trans3DIso1FC4( *G4RotationOfIso1FC4, G4PositionOfIso1FC4 ) ;

         physiIso1FC4 = new G4PVPlacement( Trans3DIso1FC4, logicIso1FC4,
       			         	      "FC4Iso1", logicWorld, false, 0 );

   // Isolator-2 : Macor                
     G4VisAttributes* FC4Iso2VisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4Iso2VisAtt->SetForceWireframe(true);
     FC4Iso2VisAtt->SetForceSolid(true);

     solidIso2FC4 = new G4Tubs("FC4Iso2",
			        elsparam->rminIso2FC4()*mm,
			        elsparam->rmaxIso2FC4()*mm,
				elsparam->lIso2FC4()*0.5*mm,
				0.*deg, 360.*deg );

      logicIso2FC4 = new G4LogicalVolume( solidIso2FC4, MaterialMacor, 
                                              "FC4Iso2", 0, 0, 0 );
      logicIso2FC4->SetVisAttributes(FC4Iso2VisAtt);

      G4ThreeVector G4PositionOfIso2FC4
	   = G4ThreeVector( elsparam->positionIso2FC4(0)*mm,
			    elsparam->positionIso2FC4(1)*mm,
			    elsparam->positionIso2FC4(2)*mm );

      G4RotationMatrix* G4RotationOfIso2FC4 = new G4RotationMatrix;
	 G4RotationOfIso2FC4->rotateX( elsparam->rotationIso2FC4(0)*deg );
	 G4RotationOfIso2FC4->rotateY( elsparam->rotationIso2FC4(1)*deg );
	 G4RotationOfIso2FC4->rotateZ( elsparam->rotationIso2FC4(2)*deg );
	 G4Transform3D Trans3DIso2FC4( *G4RotationOfIso2FC4, G4PositionOfIso2FC4 ) ;

        physiIso2FC4 = new G4PVPlacement( Trans3DIso2FC4, logicIso2FC4,
        	                	      "FC4Iso2", logicWorld, false, 0 );

   // Top Plate : A5052   
     G4VisAttributes* FC4TopPlateVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));      
     FC4TopPlateVisAtt->SetForceWireframe(true);
     FC4TopPlateVisAtt->SetForceSolid(true);

     solidTopPlateFC4 = new G4Box("FC4TopPlate",
				  elsparam->sizeofTopPlateFC4(0)*0.5*mm,
				  elsparam->sizeofTopPlateFC4(1)*0.5*mm,
				  elsparam->sizeofTopPlateFC4(2)*0.5*mm );
			           
      logicTopPlateFC4 = new G4LogicalVolume( solidTopPlateFC4, MaterialtA5052, 
                                              "FC4TopPlate", 0, 0, 0 );
      logicTopPlateFC4->SetVisAttributes(FC4TopPlateVisAtt);

      G4ThreeVector G4PositionOfTopPlateFC4
	   = G4ThreeVector( elsparam->positionTopPlateFC4(0)*mm,
			    elsparam->positionTopPlateFC4(1)*mm,
			    elsparam->positionTopPlateFC4(2)*mm );

      G4RotationMatrix* G4RotationOfTopPlateFC4 = new G4RotationMatrix;
	 G4RotationOfTopPlateFC4->rotateX( elsparam->rotationTopPlateFC4(0)*deg );
	 G4RotationOfTopPlateFC4->rotateY( elsparam->rotationTopPlateFC4(1)*deg );
	 G4RotationOfTopPlateFC4->rotateZ( elsparam->rotationTopPlateFC4(2)*deg );
	 G4Transform3D Trans3DTopPlateFC4( *G4RotationOfTopPlateFC4, G4PositionOfTopPlateFC4 ) ;

            physiTopPlateFC4 = new G4PVPlacement( Trans3DTopPlateFC4, logicTopPlateFC4,
        	             	      "FC4TopPlate", logicWorld, false, 0 );

    }

     //============================================================
     /* FC 5  added by 130405*/
     // made by B.K.Shin-san  
     //============================================================
    if( Faradaycup_flag == 5 && Faradaycup4_flag == 0 && ScreenMonitor3_flag == 0 ){

     // Body : Copper beam dump (G4Polycone) No 05            
     numZPlanesBodyFC5=elsparam->numZPlanesBodyFC5();
     for( int i=0; i<numZPlanesBodyFC5; i++ ){
       irBodyFC5[i]=elsparam->irBodyFC5(i);
       orBodyFC5[i]=elsparam->orBodyFC5(i);
       zBodyFC5[i]=elsparam->zBodyFC5(i);
     }

     G4VisAttributes* FC5BodyVisAtt = new G4VisAttributes(true,G4Colour(0.9,0.2,0.2));
     FC5BodyVisAtt->SetForceWireframe(true);
     FC5BodyVisAtt->SetForceSolid(true);

     solidBodyFC5 = new G4Polycone("FC5Body",
				   00.0,360.0*deg,
				   numZPlanesBodyFC5, zBodyFC5, irBodyFC5, orBodyFC5 );

     fCopper = man->FindOrBuildMaterial(name = "G4_Cu");
     fTitanium = man->FindOrBuildMaterial(name = "G4_Ti");

     logicBodyFC5 = new G4LogicalVolume( solidBodyFC5, MaterialOxgenFreeCopper, 
                                          //  MaterialBeCu --> must be "MaterialOxgenFreeCopper"
					 "FC5Body", 0, 0, 0 );
     logicBodyFC5->SetVisAttributes(FC5BodyVisAtt);

        G4ThreeVector G4PositionOfBodyFC5
          = G4ThreeVector( elsparam->positionBodyFC5(0)*mm,
                           elsparam->positionBodyFC5(1)*mm,
                           elsparam->positionBodyFC5(2)*mm );

        G4RotationMatrix* G4RotationOfBodyFC5 = new G4RotationMatrix;
        G4RotationOfBodyFC5->rotateX( elsparam->rotationBodyFC5(0)*deg );
        G4RotationOfBodyFC5->rotateY( elsparam->rotationBodyFC5(1)*deg );
        G4RotationOfBodyFC5->rotateZ( elsparam->rotationBodyFC5(2)*deg );
        G4Transform3D Trans3DBodyFC5( *G4RotationOfBodyFC5, G4PositionOfBodyFC5 ) ;

	// BS : berrel supporter Teflon (G4Polycone) No 04                                                                                               
        numZPlanesBSFC5=elsparam->numZPlanesBSFC5();
        for( int i=0; i<numZPlanesBSFC5; i++ ){
          irBSFC5[i]=elsparam->irBSFC5(i);
          orBSFC5[i]=elsparam->orBSFC5(i);
          zBSFC5[i]=elsparam->zBSFC5(i);
        }

        G4VisAttributes* FC5BSVisAtt = new G4VisAttributes(true,G4Colour(0.80,0.8,0.8));
        FC5BSVisAtt->SetForceWireframe(true);
        FC5BSVisAtt->SetForceSolid(true);

        solidBSFC5 = new G4Polycone("FC5BS",
                                    0.*deg,360.*deg,
                                    numZPlanesBSFC5, zBSFC5, irBSFC5, orBSFC5 );

	logicBSFC5 = new G4LogicalVolume( solidBSFC5,C2F4,
					  "FC5BS", 0, 0, 0 );
        logicBSFC5->SetVisAttributes(FC5BSVisAtt);

	G4ThreeVector G4PositionOfBSFC5
          = G4ThreeVector( elsparam->positionBSFC5(0)*mm,
                           elsparam->positionBSFC5(1)*mm,
                           elsparam->positionBSFC5(2)*mm );

        G4RotationMatrix* G4RotationOfBSFC5 = new G4RotationMatrix;
        G4RotationOfBSFC5->rotateX( elsparam->rotationBSFC5(0)*deg );
        G4RotationOfBSFC5->rotateY( elsparam->rotationBSFC5(1)*deg );
        G4RotationOfBSFC5->rotateZ( elsparam->rotationBSFC5(2)*deg );
	G4Transform3D Trans3DBSFC5( *G4RotationOfBSFC5, G4PositionOfBSFC5 ) ;

	// BSV : berrel supporter Vaccum  (G4Tubs)   
        G4VisAttributes* FC5BSVVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.,0.));
        FC5BSVVisAtt->SetForceWireframe(true);
        FC5BSVVisAtt->SetForceSolid(true);

        solidBSVFC5 = new G4Tubs("FC5BSV",
                                 elsparam->rminBSVFC5(),elsparam->rmaxBSVFC5(),elsparam->hBSVFC5(),
                                 0.*deg,360.*deg);

        logicBSVFC5 = new G4LogicalVolume( solidBSVFC5, Vacuum,
					   "FC5BSV", 0, 0, 0 );
        //      logicBSVFC5->SetVisAttributes(FC5BSVVisAtt);                                                                                             

	G4ThreeVector G4PositionOfBSVFC5
          = G4ThreeVector( elsparam->positionBSVFC5(0)*mm,
                           elsparam->positionBSVFC5(1)*mm,
                           elsparam->positionBSVFC5(2)*mm );

        G4RotationMatrix* G4RotationOfBSVFC5 = new G4RotationMatrix;
	G4RotationOfBSVFC5->rotateX( elsparam->rotationBSVFC5(0)*deg );
	G4RotationOfBSVFC5->rotateY( elsparam->rotationBSVFC5(1)*deg );
	G4RotationOfBSVFC5->rotateZ( elsparam->rotationBSVFC5(2)*deg );
	G4Transform3D Trans3DBSVFC5( *G4RotationOfBSVFC5, G4PositionOfBSVFC5 ) ;
	G4VisAttributes* FC5BSVBVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.,0.));
        FC5BSVBVisAtt->SetForceWireframe(true);
        FC5BSVBVisAtt->SetForceSolid(true);

        solidBSVBFC5 = new G4Tubs("FC5BSVB",
				  elsparam->rminBSVBFC5(),elsparam->rmaxBSVBFC5(),elsparam->hBSVBFC5(),
				  0.*deg,360.*deg);

        logicBSVBFC5 = new G4LogicalVolume( solidBSVBFC5, Vacuum,
                                            "FC5BSVB", 0, 0, 0 );


        //      logicBSVBFC5->SetVisAttributes(FC5BSVBVisAtt);                                                                                           

        G4ThreeVector G4PositionOfBSVBFC5
          = G4ThreeVector( elsparam->positionBSVBFC5(0)*mm,
                           elsparam->positionBSVBFC5(1)*mm,
                           elsparam->positionBSVBFC5(2)*mm );

        G4RotationMatrix* G4RotationOfBSVBFC5 = new G4RotationMatrix;
        G4RotationOfBSVBFC5->rotateX( elsparam->rotationBSVBFC5(0)*deg );
        G4RotationOfBSVBFC5->rotateY( elsparam->rotationBSVBFC5(1)*deg );
        G4RotationOfBSVBFC5->rotateZ( elsparam->rotationBSVBFC5(2)*deg );
        G4Transform3D Trans3DBSVBFC5( *G4RotationOfBSVBFC5, G4PositionOfBSVBFC5 ) ;

	// TS : Titan sheild No 07  (G4Polycone)   
        numZPlanesTSFC5=elsparam->numZPlanesTSFC5();
        for( int i=0; i<numZPlanesTSFC5; i++ ){
          irTSFC5[i]=elsparam->irTSFC5(i);
          orTSFC5[i]=elsparam->orTSFC5(i);
          zTSFC5[i]=elsparam->zTSFC5(i);
        }

        G4VisAttributes* FC5TSVisAtt = new G4VisAttributes(true,G4Colour(0.00,0.0,0.8));
        FC5TSVisAtt->SetForceWireframe(true);
        FC5TSVisAtt->SetForceSolid(true);

        solidTSFC5 = new G4Polycone("FC5TS",
                                    0.*deg,360.*deg,
                                    numZPlanesTSFC5, zTSFC5, irTSFC5, orTSFC5 );

        logicTSFC5 = new G4LogicalVolume( solidTSFC5, MaterialTi,
					  "FC5TS", 0, 0, 0 );
        logicTSFC5->SetVisAttributes(FC5TSVisAtt);

        G4ThreeVector G4PositionOfTSFC5
          = G4ThreeVector( elsparam->positionTSFC5(0)*mm,
                           elsparam->positionTSFC5(1)*mm,
                           elsparam->positionTSFC5(2)*mm );

        G4RotationMatrix* G4RotationOfTSFC5 = new G4RotationMatrix;
        G4RotationOfTSFC5->rotateX( elsparam->rotationTSFC5(0)*deg );
        G4RotationOfTSFC5->rotateY( elsparam->rotationTSFC5(1)*deg );
        G4RotationOfTSFC5->rotateZ( elsparam->rotationTSFC5(2)*deg );
        G4Transform3D Trans3DTSFC5( *G4RotationOfTSFC5, G4PositionOfTSFC5 ) ;


	// TSV : Titan sheild Vaccum  (G4Cons) 
        G4VisAttributes* FC5TSVVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.,0.));
        FC5TSVVisAtt->SetForceWireframe(true);
        FC5TSVVisAtt->SetForceSolid(true);

        solidTSVFC5 = new G4Cons("FC5TSV",
                                 elsparam->brminTSVFC5(),elsparam->brmaxTSVFC5(),
                                 elsparam->trminTSVFC5(),elsparam->trmaxTSVFC5(),elsparam->hTSVFC5(),
                                 0.*deg,360.*deg);

        logicTSVFC5 = new G4LogicalVolume( solidTSVFC5, Vacuum,
					   "FC5TSV", 0, 0, 0 );
        //      logicTSVFC5->SetVisAttributes(FC5TSVVisAtt);  

        G4ThreeVector G4PositionOfTSVFC5
          = G4ThreeVector( elsparam->positionTSVFC5(0)*mm,
                           elsparam->positionTSVFC5(1)*mm,
                           elsparam->positionTSVFC5(2)*mm );

        G4RotationMatrix* G4RotationOfTSVFC5 = new G4RotationMatrix;
        G4RotationOfTSVFC5->rotateX( elsparam->rotationTSVFC5(0)*deg );
        G4RotationOfTSVFC5->rotateY( elsparam->rotationTSVFC5(1)*deg );
        G4RotationOfTSVFC5->rotateZ( elsparam->rotationTSVFC5(2)*deg );
        G4Transform3D Trans3DTSVFC5( *G4RotationOfTSVFC5, G4PositionOfTSVFC5 ) ;

	// BotS : Bottom Supporter No 06  (G4Polycone)  
        numZPlanesBotSFC5=elsparam->numZPlanesBotSFC5();
        for( int i=0; i<numZPlanesBotSFC5; i++ ){
          irBotSFC5[i]=elsparam->irBotSFC5(i);
          orBotSFC5[i]=elsparam->orBotSFC5(i);
          zBotSFC5[i]=elsparam->zBotSFC5(i);
        }

        G4VisAttributes* FC5BotSVisAtt = new G4VisAttributes(true,G4Colour(.90,0.9,0.8));
        FC5BotSVisAtt->SetForceWireframe(true);
        FC5BotSVisAtt->SetForceSolid(true);

        solidBotSFC5 = new G4Polycone("FC5BotS",
				      0.*deg,360.*deg,
				      numZPlanesBotSFC5, zBotSFC5, irBotSFC5, orBotSFC5 );

        logicBotSFC5 = new G4LogicalVolume( solidBotSFC5, C2F4,
                                            "FC5BotS", 0, 0, 0 );
        logicBotSFC5->SetVisAttributes(FC5BotSVisAtt);

        G4ThreeVector G4PositionOfBotSFC5
          = G4ThreeVector( elsparam->positionBotSFC5(0)*mm,
                           elsparam->positionBotSFC5(1)*mm,
                           elsparam->positionBotSFC5(2)*mm );

        G4RotationMatrix* G4RotationOfBotSFC5 = new G4RotationMatrix;
        G4RotationOfBotSFC5->rotateX( elsparam->rotationBotSFC5(0)*deg );
        G4RotationOfBotSFC5->rotateY( elsparam->rotationBotSFC5(1)*deg );
        G4RotationOfBotSFC5->rotateZ( elsparam->rotationBotSFC5(2)*deg );
        G4Transform3D Trans3DBotSFC5( *G4RotationOfBotSFC5, G4PositionOfBotSFC5 ) ;



	// TopS : Top Supporter No 15  (G4Polycone)
        numZPlanesTopSFC5=elsparam->numZPlanesTopSFC5();
        for( int i=0; i<numZPlanesTopSFC5; i++ ){
          irTopSFC5[i]=elsparam->irTopSFC5(i);
          orTopSFC5[i]=elsparam->orTopSFC5(i);
          zTopSFC5[i]=elsparam->zTopSFC5(i);
        }

        G4VisAttributes* FC5TopSVisAtt = new G4VisAttributes(true,G4Colour(.90,0.9,0.8));
        FC5TopSVisAtt->SetForceWireframe(true);
        FC5TopSVisAtt->SetForceSolid(true);

        solidTopSFC5 = new G4Polycone("FC5TopS",
				      0.*deg,360.*deg,
				      numZPlanesTopSFC5, zTopSFC5, irTopSFC5, orTopSFC5 );

        logicTopSFC5 = new G4LogicalVolume( solidTopSFC5, C2F4,
                                            "FC5TopS", 0, 0, 0 );
        logicTopSFC5->SetVisAttributes(FC5TopSVisAtt);

        G4ThreeVector G4PositionOfTopSFC5
          = G4ThreeVector( elsparam->positionTopSFC5(0)*mm,
                           elsparam->positionTopSFC5(1)*mm,
                           elsparam->positionTopSFC5(2)*mm );

        G4RotationMatrix* G4RotationOfTopSFC5 = new G4RotationMatrix;
        G4RotationOfTopSFC5->rotateX( elsparam->rotationTopSFC5(0)*deg );
        G4RotationOfTopSFC5->rotateY( elsparam->rotationTopSFC5(1)*deg );
        G4RotationOfTopSFC5->rotateZ( elsparam->rotationTopSFC5(2)*deg );
        G4Transform3D Trans3DTopSFC5( *G4RotationOfTopSFC5, G4PositionOfTopSFC5 ) ;


	// TSV : Top supporter Vaccum  Bottom(G4Tubs)  
        G4VisAttributes* FC5TopSVVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.,0.));
        FC5TopSVVisAtt->SetForceWireframe(true);
        FC5TopSVVisAtt->SetForceSolid(true);

        solidTopSVFC5 = new G4Tubs("FC5TopSV",
				   elsparam->rminTopSVFC5(),elsparam->rmaxTopSVFC5(),elsparam->hTopSVFC5(),
				   0.*deg,360.*deg);

        logicTopSVFC5 = new G4LogicalVolume( solidTopSVFC5, Vacuum,
					     "FC5TopSV", 0, 0, 0 );
        //      logicTopSVFC5->SetVisAttributes(FC5TopSVVisAtt);                                                                                         

        G4ThreeVector G4PositionOfTopSVFC5
          = G4ThreeVector( elsparam->positionTopSVFC5(0)*mm,
                           elsparam->positionTopSVFC5(1)*mm,
                           elsparam->positionTopSVFC5(2)*mm );

        G4RotationMatrix* G4RotationOfTopSVFC5 = new G4RotationMatrix;
        G4RotationOfTopSVFC5->rotateX( elsparam->rotationTopSVFC5(0)*deg );
        G4RotationOfTopSVFC5->rotateY( elsparam->rotationTopSVFC5(1)*deg );
        G4RotationOfTopSVFC5->rotateZ( elsparam->rotationTopSVFC5(2)*deg );
        G4Transform3D Trans3DTopSVFC5( *G4RotationOfTopSVFC5, G4PositionOfTopSVFC5 ) ;


	// VDV :  Vacuum  Duct Vacuum(G4Tubs)    
        G4VisAttributes* FC5VDVVisAtt = new G4VisAttributes(true,G4Colour(0.0,0.,0.));
        FC5VDVVisAtt->SetForceWireframe(true);
        FC5VDVVisAtt->SetForceSolid(true);

        solidVDVFC5 = new G4Tubs("FC5VDV",
                                 elsparam->rminVDVFC5(),elsparam->rmaxVDVFC5(),elsparam->hVDVFC5(),
                                 0.*deg,360.*deg);

        logicVDVFC5 = new G4LogicalVolume( solidVDVFC5, Vacuum, 
					   "FC5VDV", 0, 0, 0 );
        //      logicVDVFC5->SetVisAttributes(FC5VDVVisAtt);    

        G4ThreeVector G4PositionOfVDVFC5
          = G4ThreeVector( elsparam->positionVDVFC5(0)*mm,
                           elsparam->positionVDVFC5(1)*mm,
                           elsparam->positionVDVFC5(2)*mm );

        G4RotationMatrix* G4RotationOfVDVFC5 = new G4RotationMatrix;
        G4RotationOfVDVFC5->rotateX( elsparam->rotationVDVFC5(0)*deg );
        G4RotationOfVDVFC5->rotateY( elsparam->rotationVDVFC5(1)*deg );
        G4RotationOfVDVFC5->rotateZ( elsparam->rotationVDVFC5(2)*deg );
        G4Transform3D Trans3DVDVFC5( *G4RotationOfVDVFC5, G4PositionOfVDVFC5 ) ;

	// VD :  Vacuum  Pipe (G4Tubs) No ??
        G4VisAttributes* FC5VDVisAtt = new G4VisAttributes(true,G4Colour(0.7,0.7,0.7));
        FC5VDVisAtt->SetForceWireframe(true);
        FC5VDVisAtt->SetForceSolid(true);

        solidVDFC5 = new G4Tubs("FC5VD",
				elsparam->rminVDFC5(),elsparam->rmaxVDFC5(),elsparam->hVDFC5(),
				0.*deg,360.*deg);

        logicVDFC5 = new G4LogicalVolume( solidVDFC5, MaterialAl,
					  "FC5VD", 0, 0, 0 );
        logicVDFC5->SetVisAttributes(FC5VDVisAtt);

        G4ThreeVector G4PositionOfVDFC5
          = G4ThreeVector( elsparam->positionVDFC5(0)*mm,
                           elsparam->positionVDFC5(1)*mm,
                           elsparam->positionVDFC5(2)*mm );

        G4RotationMatrix* G4RotationOfVDFC5 = new G4RotationMatrix;
        G4RotationOfVDFC5->rotateX( elsparam->rotationVDFC5(0)*deg );
        G4RotationOfVDFC5->rotateY( elsparam->rotationVDFC5(1)*deg );
        G4RotationOfVDFC5->rotateZ( elsparam->rotationVDFC5(2)*deg );
        G4Transform3D Trans3DVDFC5( *G4RotationOfVDFC5, G4PositionOfVDFC5 ) ;


	// CF90 CF90 Frange No  02  (G4Polycone)
        numZPlanesCF90FC5=elsparam->numZPlanesCF90FC5();
        for( int i=0; i<numZPlanesCF90FC5; i++ ){
          irCF90FC5[i]=elsparam->irCF90FC5(i);
          orCF90FC5[i]=elsparam->orCF90FC5(i);
          zCF90FC5[i]=elsparam->zCF90FC5(i);
        }


        G4VisAttributes* FC5CF90VisAtt = new G4VisAttributes(true,G4Colour(.50,0.5,0.5));
        FC5CF90VisAtt->SetForceWireframe(true);
        FC5CF90VisAtt->SetForceSolid(true);

        solidCF90FC5 = new G4Polycone("FC5CF90",
				      0.*deg,360.*deg,
				      numZPlanesCF90FC5, zCF90FC5, irCF90FC5, orCF90FC5 );

        G4SubtractionSolid *fsolidCF90FC5 
           = new G4SubtractionSolid( "FC5CF90withVacDuct", 
                                      solidCF90FC5, solidVDVFC5,0,G4ThreeVector(0,-20,0));

        logicCF90FC5 = new G4LogicalVolume( fsolidCF90FC5, SUS304,
                                            "FC5CF90", 0, 0, 0 );
        logicCF90FC5->SetVisAttributes(FC5CF90VisAtt);
        G4ThreeVector G4PositionOfCF90FC5
          = G4ThreeVector( elsparam->positionCF90FC5(0)*mm,
                           elsparam->positionCF90FC5(1)*mm,
                           elsparam->positionCF90FC5(2)*mm );

        G4RotationMatrix* G4RotationOfCF90FC5 = new G4RotationMatrix;
        G4RotationOfCF90FC5->rotateX( elsparam->rotationCF90FC5(0)*deg );
        G4RotationOfCF90FC5->rotateY( elsparam->rotationCF90FC5(1)*deg );
        G4RotationOfCF90FC5->rotateZ( elsparam->rotationCF90FC5(2)*deg );
        G4Transform3D Trans3DCF90FC5( *G4RotationOfCF90FC5, G4PositionOfCF90FC5 ) ;

	// TopSupporter  No 06b    (G4Polycone)  
        numZPlanesbTopSFC5=elsparam->numZPlanesbTopSFC5();
        for( int i=0; i<numZPlanesbTopSFC5; i++ ){
          irbTopSFC5[i]=elsparam->irbTopSFC5(i);
          orbTopSFC5[i]=elsparam->orbTopSFC5(i);
          zbTopSFC5[i]=elsparam->zbTopSFC5(i);
        }

        G4VisAttributes* FC5bTopSVisAtt = new G4VisAttributes(true,G4Colour(.90,0.9,0.9));
        FC5bTopSVisAtt->SetForceWireframe(true);
        FC5bTopSVisAtt->SetForceSolid(true);

        solidbTopSFC5 = new G4Polycone("FC5bTopS",
				       0.*deg,360.*deg,
				       numZPlanesbTopSFC5, zbTopSFC5, irbTopSFC5, orbTopSFC5 );

        logicbTopSFC5 = new G4LogicalVolume( solidbTopSFC5,C2F4
                                             ,"FC5bTopS", 0, 0, 0 );
        logicbTopSFC5->SetVisAttributes(FC5bTopSVisAtt);

        G4ThreeVector G4PositionOfbTopSFC5
          = G4ThreeVector( elsparam->positionbTopSFC5(0)*mm,
                           elsparam->positionbTopSFC5(1)*mm,
                           elsparam->positionbTopSFC5(2)*mm );

        G4RotationMatrix* G4RotationOfbTopSFC5 = new G4RotationMatrix;
        G4RotationOfbTopSFC5->rotateX( elsparam->rotationbTopSFC5(0)*deg );
        G4RotationOfbTopSFC5->rotateY( elsparam->rotationbTopSFC5(1)*deg );
        G4RotationOfbTopSFC5->rotateZ( elsparam->rotationbTopSFC5(2)*deg );
        G4Transform3D Trans3DbTopSFC5( *G4RotationOfbTopSFC5, G4PositionOfbTopSFC5 ) ;


	// Copper Chamber No 01    (G4Polycone)   
        numZPlanesCPFC5=elsparam->numZPlanesCPFC5();
        for( int i=0; i<numZPlanesCPFC5; i++ ){
          irCPFC5[i]=elsparam->irCPFC5(i);
          orCPFC5[i]=elsparam->orCPFC5(i);
          zCPFC5[i]=elsparam->zCPFC5(i);
        }

        G4VisAttributes* FC5CPVisAtt = new G4VisAttributes(true,G4Colour(.90,0.0,0.0));
        FC5CPVisAtt->SetForceWireframe(true);
        FC5CPVisAtt->SetForceSolid(true);

        solidCPFC5 = new G4Polycone("FC5CP",
                                    0.*deg,360.*deg,
                                    numZPlanesCPFC5, zCPFC5, irCPFC5, orCPFC5 );

        logicCPFC5 = new G4LogicalVolume( solidCPFC5, MaterialOxgenFreeCopper, 
                                          // MaterialBeCu --> "MaterialOxgenFreeCopper"
					  "FC5CP", 0, 0, 0 );
        logicCPFC5->SetVisAttributes(FC5CPVisAtt);

        G4ThreeVector G4PositionOfCPFC5
          = G4ThreeVector( elsparam->positionCPFC5(0)*mm,
                           elsparam->positionCPFC5(1)*mm,
                           elsparam->positionCPFC5(2)*mm );

        G4RotationMatrix* G4RotationOfCPFC5 = new G4RotationMatrix;
        G4RotationOfCPFC5->rotateX( elsparam->rotationCPFC5(0)*deg );
        G4RotationOfCPFC5->rotateY( elsparam->rotationCPFC5(1)*deg );
        G4RotationOfCPFC5->rotateZ( elsparam->rotationCPFC5(2)*deg );
        G4Transform3D Trans3DCPFC5( *G4RotationOfCPFC5, G4PositionOfCPFC5 ) ;


	// Top Plate  No 12    (G4Box) 
        G4VisAttributes* FC5TPVisAtt = new G4VisAttributes(true,G4Colour(.50,0.5,0.5));
        FC5TPVisAtt->SetForceWireframe(true);
        FC5TPVisAtt->SetForceSolid(true);

        solidTPFC5 = new G4Box("FC5TP",
                               elsparam->xTPFC5(),  elsparam->yTPFC5(), elsparam->zTPFC5());
        G4Tubs *TPHole = new G4Tubs("TPHole",0,30.5/2.,5.,0,360.*deg);
        G4double holePosition= -(195.5)/2. +56 - 20;
        G4SubtractionSolid *fsolidTPFC5 
           = new G4SubtractionSolid( "TPFC5withHole", solidTPFC5, TPHole,0,G4ThreeVector(0,holePosition,0));
        TPHole = new G4Tubs("TPHole",0,20.0/2.,5.,0,360.*deg);
        holePosition= -(195.5)/2. +56 +6;
        fsolidTPFC5 = new  G4SubtractionSolid("TPFC5withHole", fsolidTPFC5, TPHole,0,G4ThreeVector(0,holePosition,0));
        logicTPFC5 = new G4LogicalVolume( fsolidTPFC5, MaterialAl,
					  "FC5TP", 0, 0, 0 );
        logicTPFC5->SetVisAttributes(FC5TPVisAtt);

        G4ThreeVector G4PositionOfTPFC5
          = G4ThreeVector( elsparam->positionTPFC5(0)*mm,
                           elsparam->positionTPFC5(1)*mm,
                           elsparam->positionTPFC5(2)*mm );

        G4RotationMatrix* G4RotationOfTPFC5 = new G4RotationMatrix;
        G4RotationOfTPFC5->rotateX( elsparam->rotationTPFC5(0)*deg );
        G4RotationOfTPFC5->rotateY( elsparam->rotationTPFC5(1)*deg );
        G4RotationOfTPFC5->rotateZ( elsparam->rotationTPFC5(2)*deg );
        G4Transform3D Trans3DTPFC5( *G4RotationOfTPFC5, G4PositionOfTPFC5 ) ;

	// FT :  FeedThrough(G4Tubs) No 11 
        G4VisAttributes* FC5FTVisAtt = new G4VisAttributes(true,G4Colour(1.0,0.,0.));
        FC5FTVisAtt->SetForceWireframe(true);
        FC5FTVisAtt->SetForceSolid(true);

        solidFTFC5 = new G4Tubs("FC5FT",
				elsparam->rminFTFC5(),elsparam->rmaxFTFC5(),elsparam->hFTFC5(),
				0.*deg,360.*deg);

        logicFTFC5 = new G4LogicalVolume( solidFTFC5,  MaterialBeCu,
					  "FC5FT", 0, 0, 0 );
        logicFTFC5->SetVisAttributes(FC5FTVisAtt);

        G4ThreeVector G4PositionOfFTFC5
          = G4ThreeVector( elsparam->positionFTFC5(0)*mm,
                           elsparam->positionFTFC5(1)*mm,
                           elsparam->positionFTFC5(2)*mm );

        G4RotationMatrix* G4RotationOfFTFC5 = new G4RotationMatrix;
        G4RotationOfFTFC5->rotateX( elsparam->rotationFTFC5(0)*deg );
        G4RotationOfFTFC5->rotateY( elsparam->rotationFTFC5(1)*deg );
        G4RotationOfFTFC5->rotateZ( elsparam->rotationFTFC5(2)*deg );
        G4Transform3D Trans3DFTFC5( *G4RotationOfFTFC5, G4PositionOfFTFC5 ) ;


	// FT2 :  FeedThrough 2(G4Tubs) No 11   
        G4VisAttributes* FC5FT2VisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));
        FC5FT2VisAtt->SetForceWireframe(true);
        FC5FT2VisAtt->SetForceSolid(true);

        solidFT2FC5 = new G4Tubs("FC5FT2",
                                 elsparam->rminFT2FC5(),elsparam->rmaxFT2FC5(),elsparam->hFT2FC5(),
                                 0.*deg,360.*deg);

        logicFT2FC5 = new G4LogicalVolume( solidFT2FC5, MaterialAl,
					   "FC5FT2", 0, 0, 0 );
        logicFT2FC5->SetVisAttributes(FC5FT2VisAtt);

        G4ThreeVector G4PositionOfFT2FC5
          = G4ThreeVector( elsparam->positionFT2FC5(0)*mm,
                           elsparam->positionFT2FC5(1)*mm,
                           elsparam->positionFT2FC5(2)*mm );

        G4RotationMatrix* G4RotationOfFT2FC5 = new G4RotationMatrix;
        G4RotationOfFT2FC5->rotateX( elsparam->rotationFT2FC5(0)*deg );
        G4RotationOfFT2FC5->rotateY( elsparam->rotationFT2FC5(1)*deg );
        G4RotationOfFT2FC5->rotateZ( elsparam->rotationFT2FC5(2)*deg );
        G4Transform3D Trans3DFT2FC5( *G4RotationOfFT2FC5, G4PositionOfFT2FC5 ) ;

       	physiTPFC5    = new G4PVPlacement( Trans3DTPFC5,    logicTPFC5,    "FC5TP",      logicWorld, false, 0 );
   	    physibTopSFC5 = new G4PVPlacement( Trans3DbTopSFC5, logicbTopSFC5, "FC5bTopS",   logicWorld, false, 0 );
	    physiCF90FC5  = new G4PVPlacement( Trans3DCF90FC5,  logicCF90FC5,  "FC5CF90",    logicWorld, false, 0 );
	    physiTopSFC5  = new G4PVPlacement( Trans3DTopSFC5,  logicTopSFC5,  "FC5TopS",    logicWorld, false, 0 );
	    physiBotSFC5  = new G4PVPlacement( Trans3DBotSFC5,  logicBotSFC5,  "FC5BotS",    logicWorld, false, 0 );
	    physiTSVFC5   = new G4PVPlacement( Trans3DTSVFC5,   logicTSVFC5,   "FC5TSVac",   logicWorld, false, 0 );
	    physiBSFC5    = new G4PVPlacement( Trans3DBSFC5,    logicBSFC5,    "FC5BS",      logicWorld, false, 0 );
	    physiBSVFC5   = new G4PVPlacement( Trans3DBSVFC5,   logicBSVFC5,   "FC5BSVac",   logicWorld, false, 0 );
     	physiBSVBFC5  = new G4PVPlacement( Trans3DBSVBFC5,  logicBSVBFC5,  "FC5BSBVac",  logicWorld, false, 0 );
	    physiTopSVFC5 = new G4PVPlacement( Trans3DTopSVFC5, logicTopSVFC5, "FC5TopSVac", logicWorld, false, 0 );
	    physiVDFC5    = new G4PVPlacement( Trans3DVDFC5,    logicVDFC5,    "FC5VD",      logicWorld, false, 0 );
	    physiVDVFC5   = new G4PVPlacement( Trans3DVDVFC5,   logicVDVFC5,   "FC5VDVVac",  logicWorld, false, 0 );
	    physiFT2FC5   = new G4PVPlacement( Trans3DFT2FC5,   logicFT2FC5,   "FC5FT2",     logicWorld, false, 0 );

   	    physiBodyFC5  = new G4PVPlacement( Trans3DBodyFC5,  logicBodyFC5,  "FC5Body",    logicWorld, false, 0 );  // Cu-Dump
	    physiFTFC5    = new G4PVPlacement( Trans3DFTFC5,    logicFTFC5,    "FC5Body",    logicWorld, false, 0 );  // FeedThrough

	    physiTSFC5    = new G4PVPlacement( Trans3DTSFC5,    logicTSFC5,    "FC5TS",      logicWorld, false, 0 );  // Ti-Chamber
	    physiCPFC5    = new G4PVPlacement( Trans3DCPFC5,    logicCPFC5,    "FC5CP",      logicWorld, false, 0 );  // Copper     
    }

     //--------------------------------------------------           
     // Screen Monitor 3 added in 2011.12.12
     //--------------------------------------------------          
    if( Faradaycup_flag != 1 && Faradaycup_flag!=4 && Faradaycup_flag!=5 &&
        Faradaycup4_flag == 0 && ScreenMonitor3_flag == 1 ){

     G4VisAttributes* ScreenMonitor3VisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));
     ScreenMonitor3VisAtt->SetForceWireframe(true);
     ScreenMonitor3VisAtt->SetForceSolid(true);

     solidSM3 = new G4Box("SM3", 
                          elsparam->sizeScreenMonitor3(0)*0.5*mm,
  	   	          elsparam->sizeScreenMonitor3(1)*0.5*mm,
			  elsparam->sizeScreenMonitor3(2)*0.5*mm );
     
     logicSM3 = new G4LogicalVolume( solidSM3,
				     MaterialSM3,
				     "SM3", 0, 0, 0 );

     logicSM3->SetVisAttributes( ScreenMonitor3VisAtt );

      G4ThreeVector G4PositionOfScreenMonitor3
	= G4ThreeVector( elsparam->positionScreenMonitor3(0)*mm,
			 elsparam->positionScreenMonitor3(1)*mm,
			 elsparam->positionScreenMonitor3(2)*mm );

      G4cout << "G4PositionOfScreenMonitor3 = " << G4PositionOfScreenMonitor3 << G4endl;
   
      G4RotationMatrix* G4RotationOfScreenMonitor3 = new G4RotationMatrix;
      G4RotationOfScreenMonitor3->rotateX( elsparam->rotationScreenMonitor3(0)*deg );
      G4RotationOfScreenMonitor3->rotateY( elsparam->rotationScreenMonitor3(1)*deg );
      G4RotationOfScreenMonitor3->rotateZ( elsparam->rotationScreenMonitor3(2)*deg );
      G4Transform3D Trans3DScreenMonitor3( *G4RotationOfScreenMonitor3, G4PositionOfScreenMonitor3 ) ;

      G4cout << " ScreenMonitor3_flag = " << ScreenMonitor3_flag << G4endl;

	physiSM3 = new G4PVPlacement( Trans3DScreenMonitor3, logicSM3, 
                                      "SM3",  logicWorld, false, 0 );

    }
                 
      //--------------------------------------------------
      // Beam Attenuator added in 2011.12.12       
      //--------------------------------------------------     
    if( Faradaycup_flag != 1 && Faradaycup_flag!=4 && Faradaycup_flag!=5 &&
        Faradaycup4_flag == 0 && ScreenMonitor3_flag == 0 && 
         BeamAttenuator_flag == 1 && PbCollimator_flag == 0 ) {

      G4VisAttributes*  BeamAttenuatorVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));  
      BeamAttenuatorVisAtt->SetForceWireframe(true);
      BeamAttenuatorVisAtt->SetForceSolid(true);

      solidBeamAttenuator = new G4Tubs("BeamAttenuator",
				elsparam->rminBeamAttenuator()*mm,
                                elsparam->rmaxBeamAttenuator()*mm,
                                elsparam->lBeamAttenuator()*0.5*mm,
                                0.*deg, 360.*deg ); 

      logicBeamAttenuator = new G4LogicalVolume( solidBeamAttenuator, 
					         MaterialBeamAttenuator,
                                                 "BeamAttenuator", 0, 0, 0 );

      logicBeamAttenuator->SetVisAttributes( BeamAttenuatorVisAtt );

      G4ThreeVector G4PositionOfBeamAttenuator
	   = G4ThreeVector( elsparam->positionBeamAttenuator(0)*mm,
			    elsparam->positionBeamAttenuator(1)*mm,
			    elsparam->positionBeamAttenuator(2)*mm );

      G4RotationMatrix* G4RotationOfBeamAttenuator = new G4RotationMatrix;
	 G4RotationOfBeamAttenuator->rotateX( elsparam->rotationBeamAttenuator(0)*deg );
	 G4RotationOfBeamAttenuator->rotateY( elsparam->rotationBeamAttenuator(1)*deg );
	 G4RotationOfBeamAttenuator->rotateZ( elsparam->rotationBeamAttenuator(2)*deg );
	 G4Transform3D Trans3DBeamAttenuator( *G4RotationOfBeamAttenuator, G4PositionOfBeamAttenuator ) ;

         physiBeamAttenuator = new G4PVPlacement( Trans3DBeamAttenuator, logicBeamAttenuator,
      			         	          "BeamAttenuator", logicWorld, false, 0 );

    }

      //--------------------------------------------------
      // Pb Collimator added in 2012.09.25       
      //--------------------------------------------------     

      G4VisAttributes*  PbCollimatorVisAtt = new G4VisAttributes(true,G4Colour(0.5,0.5,0.5));  
      PbCollimatorVisAtt->SetForceWireframe(true);
      PbCollimatorVisAtt->SetForceSolid(true);

      solidPbCollimator = new G4Tubs("PbCollimator",
 	                  			elsparam->rminPbCollimator()*mm,
                                elsparam->rmaxPbCollimator()*mm,
                                elsparam->lPbCollimator()*0.5*mm,
                                0.*deg, 360.*deg ); 

      logicPbCollimator = new G4LogicalVolume( solidPbCollimator, 
					       MaterialPbCollimator,
                                               "PbCollimator", 0, 0, 0 );

      logicPbCollimator->SetVisAttributes( PbCollimatorVisAtt );

      G4ThreeVector G4PositionOfPbCollimator
	   = G4ThreeVector( elsparam->positionPbCollimator(0)*mm,
			    elsparam->positionPbCollimator(1)*mm,
			    elsparam->positionPbCollimator(2)*mm );

      G4RotationMatrix* G4RotationOfPbCollimator = new G4RotationMatrix;
	 G4RotationOfPbCollimator->rotateX( elsparam->rotationPbCollimator(0)*deg );
	 G4RotationOfPbCollimator->rotateY( elsparam->rotationPbCollimator(1)*deg );
	 G4RotationOfPbCollimator->rotateZ( elsparam->rotationPbCollimator(2)*deg );
	 G4Transform3D Trans3DPbCollimator( *G4RotationOfPbCollimator, G4PositionOfPbCollimator ) ;

     if( BeamAttenuator_flag == 0 && PbCollimator_flag == 1 ) {
         physiPbCollimator = new G4PVPlacement( Trans3DPbCollimator, logicPbCollimator,
      			         	        "PbCollimator", logicWorld, false, 0 );
     }

     //-------------------------------------------------------------   
     // Virtual Test Chamber added in 2011.12.13
     //   modified in 2013.04.17, to check backscattering effect            
     //                     all of gemetory parameters are fixed            
     //-------------------------------------------------------------
     G4VisAttributes* VirtualChamberVisAtt = new G4VisAttributes(true,G4Colour(0.1,0.1,0.1));
     VirtualChamberVisAtt->SetForceWireframe(true);
     VirtualChamberVisAtt->SetForceSolid(true);

     // Cylinder( Material ) parts 1-5 ---> 1-3 in 2013.04.17 --> 1-4 in 2013.06.26 --> 1-5 : in 2013.06.28
     char VOLUME_NAME_VirtualChamberCylinder[5][256];
     sprintf(VOLUME_NAME_VirtualChamberCylinder[0], "VirtualChamberCylinder");
     sprintf(VOLUME_NAME_VirtualChamberCylinder[1], "VirtualChamberCylinderAl");
     sprintf(VOLUME_NAME_VirtualChamberCylinder[2], "VirtualChamberCylinderWall");
     sprintf(VOLUME_NAME_VirtualChamberCylinder[3], "VirtualChamberCylinderWall");
     sprintf(VOLUME_NAME_VirtualChamberCylinder[4], "VirtualChamberCylinderRear");
     for( int i=0; i<5; i++ ){
          solidVirtualChamberCylinder[i] = new G4Tubs(VOLUME_NAME_VirtualChamberCylinder[i],
						      elsparam->rminVirtualChamberCylinder(i)*mm,
						      elsparam->rmaxVirtualChamberCylinder(i)*mm,
						      elsparam->lVirtualChamberCylinder(i)*0.5*mm,
					      	      0.*deg, 360.*deg );

	  if( i!=1 ){
          logicVirtualChamberCylinder[i] = new G4LogicalVolume( solidVirtualChamberCylinder[i],
						                SUS304,
						                VOLUME_NAME_VirtualChamberCylinder[i], 0, 0, 0 );
	  }else if( i==1 ){
          logicVirtualChamberCylinder[i] = new G4LogicalVolume( solidVirtualChamberCylinder[i],
						                MaterialtPureAl,
						                VOLUME_NAME_VirtualChamberCylinder[i], 0, 0, 0 );
	  }
	  logicVirtualChamberCylinder[i]->SetVisAttributes( VirtualChamberVisAtt );
          
          G4ThreeVector G4PositionOfVirtualChamberCylinder
	    = G4ThreeVector( elsparam->positionVirtualChamberCylinder(i,0)*mm,
			     elsparam->positionVirtualChamberCylinder(i,1)*mm,
			     elsparam->positionVirtualChamberCylinder(i,2)*mm );

	  G4RotationMatrix* G4RotationOfVirtualChamberCylinder = new G4RotationMatrix;
	  G4RotationOfVirtualChamberCylinder->rotateX( elsparam->rotationVirtualChamberCylinder(i,0)*deg );
	  G4RotationOfVirtualChamberCylinder->rotateY( elsparam->rotationVirtualChamberCylinder(i,1)*deg );
	  G4RotationOfVirtualChamberCylinder->rotateZ( elsparam->rotationVirtualChamberCylinder(i,2)*deg );

	  G4Transform3D Trans3DVirtualChamberCylinder( *G4RotationOfVirtualChamberCylinder, 
                                                        G4PositionOfVirtualChamberCylinder );

        if( Virtual_chamber_flag == 1 ) {
	    physiVirtualChamberCylinder[i] = new G4PVPlacement( Trans3DVirtualChamberCylinder, 
                                                                logicVirtualChamberCylinder[i],
						                VOLUME_NAME_VirtualChamberCylinder[i],
                                                                logicWorld, false, 0 );
	}

      }  // O.K,.


     // Target  : added in 2013.04.17, to check backscattering effect
      solidVirtualChamberTarget = new G4Tubs("VirtualChamberTarget",
                                             elsparam->rminVirtualChamberTarget()*mm,
                                             elsparam->rmaxVirtualChamberTarget()*mm,
                                             elsparam->lVirtualChamberTarget()*0.5*mm,
                                             0.*deg, 360.*deg );

      logicVirtualChamberTarget = new G4LogicalVolume( solidVirtualChamberTarget, MaterialToughPitchCopper, 
                                                       "VirtualChamberTarget", 0, 0, 0 );

      G4ThreeVector G4PositionOfVirtualChamberTarget 
	= G4ThreeVector( elsparam->positionVirtualChamberTarget(0)*mm,
                         elsparam->positionVirtualChamberTarget(1)*mm,
			 elsparam->positionVirtualChamberTarget(2)*mm );

      G4RotationMatrix* G4RotationOfVirtualChamberTarget = new G4RotationMatrix;
      G4RotationOfVirtualChamberTarget->rotateX( elsparam->rotationVirtualChamberTarget(0)*deg );
      G4RotationOfVirtualChamberTarget->rotateY( elsparam->rotationVirtualChamberTarget(1)*deg );
      G4RotationOfVirtualChamberTarget->rotateZ( elsparam->rotationVirtualChamberTarget(2)*deg );

      G4Transform3D Trans3DVirtualChamberTarget( *G4RotationOfVirtualChamberTarget,
 						 G4PositionOfVirtualChamberTarget );

      if( Virtual_chamber_flag == 1 ) {
	   physiVirtualChamberTarget = new G4PVPlacement( Trans3DVirtualChamberTarget,
                                                          logicVirtualChamberTarget,
                                                          "VirtualChamberTarget", logicWorld, false, 0 );
      }
      // O.K.

      // Vacuum Region 
      /*
      solidVirtualChamberCylinder_vacbase = new G4Tubs("VirtualChamberVacBase",
						       0.0*mm, 
						       elsparam->rminVirtualChamberCylinder(1)*mm,
						       elsparam->lVirtualChamberCylinder(1)*0.5*mm,
						       0.*deg, 360.*deg );
      solidVirtualChamberCylinder_vactarget = new G4Tubs( "VirtualChamberVacTarger",
							  elsparam->rminVirtualChamberTarget()*mm,
							  elsparam->rmaxVirtualChamberTarget()*mm,
							  elsparam->lVirtualChamberTarget()*0.5*mm,
							  0.*deg, 360.*deg );
       
       G4ThreeVector G4PositionOfVirtualChamberCylinderVacBase
	    = G4ThreeVector( elsparam->positionVirtualChamberCylinder(1,0)*mm,
			     elsparam->positionVirtualChamberCylinder(1,1)*mm,
			     elsparam->positionVirtualChamberCylinder(1,2)*mm );
      G4ThreeVector G4PositionOfVirtualChamberVacTarget 
   	    = G4ThreeVector( elsparam->positionVirtualChamberTarget(0)*mm,
                             elsparam->positionVirtualChamberTarget(1)*mm,
		    	     elsparam->positionVirtualChamberTarget(2)*mm );

      solidVirtualChamberCylinder_vac = new G4SubtractionSolid( "VirtualChamberVac",
						                solidVirtualChamberCylinder_vacbase,
                                                                solidVirtualChamberCylinder_vactarget,
                     					        0,
						                G4PositionOfVirtualChamberVacTarget - G4PositionOfVirtualChamberCylinderVacBase );

      logicVirtualChamberCylinder_vac = new G4LogicalVolume( solidVirtualChamberCylinder_vac, Vacuum,
                                                             "VirtualChamberVac", 0, 0, 0 );

      G4RotationMatrix* G4RotationOfVirtualChamberVac = new G4RotationMatrix;
      G4RotationOfVirtualChamberVac->rotateX( 0.0*deg );
      G4RotationOfVirtualChamberVac->rotateY( 0.0*deg );
      G4RotationOfVirtualChamberVac->rotateZ( 0.0*deg );

      G4Transform3D Trans3DVirtualChamberVac( *G4RotationOfVirtualChamberVac,
					       G4PositionOfVirtualChamberVacTarget );

      if( Virtual_chamber_flag == 1 ) {
	physiVirtualChamberCylinder_vac = new G4PVPlacement( Trans3DVirtualChamberVac,
							     logicVirtualChamberCylinder_vac,
                                                              "VirtualChamberVac", logicWorld, false, 0 );
      }
      */

     // Body parts 1-4
     char VOLUME_NAME_VirtualChamberBody[4][256];
     sprintf(VOLUME_NAME_VirtualChamberBody[0], "VirtualChamberBody");
     sprintf(VOLUME_NAME_VirtualChamberBody[1], "VirtualChamberPipe");
     sprintf(VOLUME_NAME_VirtualChamberBody[2], "VirtualChamberBody");
     sprintf(VOLUME_NAME_VirtualChamberBody[3], "VirtualChamberPipe");
     for( int i=0; i<4; i++ ){
      solidVirtualChamberBody[i] = new G4Tubs( "VirtualChamberBody",
                                               elsparam->rminVirtualChamberBody(i)*mm,
                                               elsparam->rmaxVirtualChamberBody(i)*mm,
                                               elsparam->lVirtualChamberBody(i)*0.5*mm,
	  				       0.*deg, 360.*deg );
      logicVirtualChamberBody[i]  = new G4LogicalVolume( solidVirtualChamberBody[i],
						         N2Air,
						         VOLUME_NAME_VirtualChamberBody[i], 0, 0, 0 );

      logicVirtualChamberBody[i]->SetVisAttributes( VirtualChamberVisAtt );

      G4ThreeVector G4PositionOfVirtualChamberBody 
	= G4ThreeVector( elsparam->positionVirtualChamberBody(i,0)*mm,
                         elsparam->positionVirtualChamberBody(i,1)*mm,
			 elsparam->positionVirtualChamberBody(i,2)*mm );

      G4RotationMatrix* G4RotationOfVirtualChamberBody = new G4RotationMatrix;
      G4RotationOfVirtualChamberBody->rotateX( elsparam->rotationVirtualChamberBody(i,0)*deg );
      G4RotationOfVirtualChamberBody->rotateY( elsparam->rotationVirtualChamberBody(i,1)*deg );
      G4RotationOfVirtualChamberBody->rotateZ( elsparam->rotationVirtualChamberBody(i,2)*deg );

      G4Transform3D Trans3DVirtualChamberBody( *G4RotationOfVirtualChamberBody,
						G4PositionOfVirtualChamberBody );
      
      if( Virtual_chamber_flag == 1 ) {
	/*
	  physiVirtualChamberBody[i] = new G4PVPlacement( Trans3DVirtualChamberBody,
		  			  	          logicVirtualChamberBody[i],
						          VOLUME_NAME_VirtualChamberBody[i], logicWorld, false, 0 );
	*/
      } 

     }
   //-------------------------------------------------
   // Pb Block
   //-------------------------------------------------

     G4VisAttributes* PbVisAtt = new G4VisAttributes(true,G4Colour(0.1,0.1,0.1));  
      PbVisAtt->SetForceWireframe(true);
      PbVisAtt->SetForceSolid(true);

      for( int i=0; i<5; i++ ){

        solidPbTower[i]= new G4Box( "Pb Tower", 
			            elsparam->sizePbTower(i,0)*0.5*mm,
                                    elsparam->sizePbTower(i,1)*0.5*mm,
                                    elsparam->sizePbTower(i,2)*0.5*mm );
        logicPbTower[i] = new G4LogicalVolume( solidPbTower[i], MaterialPbBlock,
					     "Pb Tower",  0, 0, 0 );
        //logicPbTower[i]->SetVisAttributes(PbVisAtt);

        G4ThreeVector G4PositionOfPbTower
	  = G4ThreeVector( elsparam->positionPbTower(i,0)*mm,
                           elsparam->positionPbTower(i,1)*mm,
			   elsparam->positionPbTower(i,2)*mm ); 
        G4RotationMatrix* G4RotationOfPbTower= new G4RotationMatrix;
  	  G4RotationOfPbTower->rotateX( 0.0*deg );
  	  G4RotationOfPbTower->rotateY( 0.0*deg );
	  G4RotationOfPbTower->rotateZ( 0.0*deg );
	G4Transform3D Trans3DPbTower( *G4RotationOfPbTower, G4PositionOfPbTower );

        physiPbTower[i] = new G4PVPlacement( Trans3DPbTower, logicPbTower[i] ,
 					     "Pb Tower", logicWorld, false, 0 );
      }


      for( int i=0; i<2; i++ ){

	// parts-1
       solidPbBlockBMparts1[i] = new G4Tubs("Pb BM",
                                            elsparam->rminPbBlockBMparts1(i)*mm, 
 			 		    elsparam->rmaxPbBlockBMparts1(i)*mm, 
                                            elsparam->lPbBlockparts1(i)*0.5*mm,
                                            elsparam->phistartPbBlockBMparts1(i)*deg, 
                                            elsparam->phideltaPbBlockBMparts1(i)*deg ); 
       logicPbBlockBMparts1[i] = new G4LogicalVolume( solidPbBlockBMparts1[i],
                                                          MaterialPbBlock,
                                                          "Pb BM", 0, 0, 0 );
       logicPbBlockBMparts1[i]->SetVisAttributes(PbVisAtt);

        G4ThreeVector G4PositionOfPbBlockBMparts1
	  = G4ThreeVector( elsparam->positionPbBlockBMparts1(i,0)*mm,
                           elsparam->positionPbBlockBMparts1(i,1)*mm,
			   elsparam->positionPbBlockBMparts1(i,2)*mm ); 
        G4RotationMatrix* G4RotationOfPbBlockBMparts1 = new G4RotationMatrix;
  	  G4RotationOfPbBlockBMparts1->rotateX( elsparam->rotationPbBlockBMparts1(i,0)*deg );
  	  G4RotationOfPbBlockBMparts1->rotateY( elsparam->rotationPbBlockBMparts1(i,1)*deg );
	  G4RotationOfPbBlockBMparts1->rotateZ( elsparam->rotationPbBlockBMparts1(i,2)*deg );
	G4Transform3D Trans3DPbBlockBMparts1( *G4RotationOfPbBlockBMparts1, G4PositionOfPbBlockBMparts1 );
   

	G4cout << " G4PositionOfPbBlockBMparts1 = " << G4PositionOfPbBlockBMparts1 << G4endl;

        physiPbBlockBMparts1[i] = new G4PVPlacement(  Trans3DPbBlockBMparts1,
                                                     logicPbBlockBMparts1[i],
                                                     "Pb BM",  logicWorld, false, 0 );

	// parts-2
	solidPbBlockBMparts2[i] = new G4Box("Pb BM",
					  elsparam->sizePbBlockBMparts2(i,0)*0.5*mm,
                                          elsparam->sizePbBlockBMparts2(i,1)*0.5*mm,
					  elsparam->sizePbBlockBMparts2(i,2)*0.5*mm );
	logicPbBlockBMparts2[i] = new G4LogicalVolume( solidPbBlockBMparts2[i],
							   MaterialPbBlock,
							   "Pb BM", 0, 0, 0 );
	logicPbBlockBMparts2[i]->SetVisAttributes(PbVisAtt);

        G4ThreeVector G4PositionOfPbBlockBMparts2
          = G4ThreeVector( elsparam->positionPbBlockBMparts2(i,0)*mm,
                           elsparam->positionPbBlockBMparts2(i,1)*mm,
                           elsparam->positionPbBlockBMparts2(i,2)*mm );
        G4RotationMatrix* G4RotationOfPbBlockBMparts2 = new G4RotationMatrix;
	  G4RotationOfPbBlockBMparts2->rotateX( 0.*deg );
  	  G4RotationOfPbBlockBMparts2->rotateY( 0.*deg );
	  G4RotationOfPbBlockBMparts2->rotateZ( 0.*deg );
        G4Transform3D Trans3DPbBlockBMparts2( *G4RotationOfPbBlockBMparts2, G4PositionOfPbBlockBMparts2 );

        physiPbBlockBMparts2[i] = new G4PVPlacement(  Trans3DPbBlockBMparts2,
                                                     logicPbBlockBMparts2[i],
                                                     "Pb BM",  logicWorld, false, 0 );

	// parts-3
        solidPbBlockBMparts3[i] = new G4Box("Pb BM",
                                          elsparam->sizePbBlockBMparts3(i,0)*0.5*mm,
                                          elsparam->sizePbBlockBMparts3(i,1)*0.5*mm,
                                          elsparam->sizePbBlockBMparts3(i,2)*0.5*mm );
        logicPbBlockBMparts3[i] = new G4LogicalVolume( solidPbBlockBMparts3[i],
                                                           MaterialPbBlock,
                                                           "Pb BM", 0, 0, 0 );
        logicPbBlockBMparts3[i]->SetVisAttributes(PbVisAtt);

        G4ThreeVector G4PositionOfPbBlockBMparts3
          = G4ThreeVector( elsparam->positionPbBlockBMparts3(i,0)*mm,
                           elsparam->positionPbBlockBMparts3(i,1)*mm,
                           elsparam->positionPbBlockBMparts3(i,2)*mm );
        G4RotationMatrix* G4RotationOfPbBlockBMparts3 = new G4RotationMatrix;
 	  G4RotationOfPbBlockBMparts3->rotateX( 0.*deg );
  	  G4RotationOfPbBlockBMparts3->rotateY( 0.*deg );
  	  G4RotationOfPbBlockBMparts3->rotateZ( 0.*deg );
        G4Transform3D Trans3DPbBlockBMparts3( *G4RotationOfPbBlockBMparts3, G4PositionOfPbBlockBMparts3 );

        physiPbBlockBMparts3[i] = new G4PVPlacement(  Trans3DPbBlockBMparts3,
                                                     logicPbBlockBMparts3[i],
                                                     "Pb BM",  logicWorld, false, 0 );
      }

   //-------------------------------------------------
   // Concrete Block
   //-------------------------------------------------
     G4VisAttributes* ConcreteVisAtt = new G4VisAttributes(true,G4Colour(1.0,1.0,1.0));  
      ConcreteVisAtt->SetForceWireframe(true);
      //ConcreteVisAtt->SetForceSolid(true);

   for( int i=0; i < elsparam->numCon(); i++ ){

   solidConcreteBlock[i] = new G4Box("ConcreteBlock",
	      			elsparam->sizeConcreteBlock(i,0)*0.5*mm,
			        elsparam->sizeConcreteBlock(i,1)*0.5*mm,
			        elsparam->sizeConcreteBlock(i,2)*0.5*mm );

   logicConcreteBlock[i] = new G4LogicalVolume( solidConcreteBlock[i],
                                                MaterialConcrete,
                                                "ConcreteBlock",
                                                0,0,0 );        
   logicConcreteBlock[i]->SetVisAttributes(ConcreteVisAtt);

   G4ThreeVector G4PositionOfConcreteBlock
     = G4ThreeVector( elsparam->positionConcreteBlock(i,0)*mm,
		      elsparam->positionConcreteBlock(i,1)*mm,
		      elsparam->positionConcreteBlock(i,2)*mm );

  G4cout << " Concrete Block " << G4endl;
       G4cout << " Xmin , Xmax " 
              <<  " " << ( elsparam->positionConcreteBlock(i,0) - elsparam->sizeConcreteBlock(i,0)/2.0 )/1000.0
              <<  " " << ( elsparam->positionConcreteBlock(i,0) + elsparam->sizeConcreteBlock(i,0)/2.0 )/1000.0 
              << G4endl;
       G4cout << " Ymin , Ymax " 
              <<  " " << ( elsparam->positionConcreteBlock(i,1) - elsparam->sizeConcreteBlock(i,1)/2.0 )/1000.0
              <<  " " << ( elsparam->positionConcreteBlock(i,1) + elsparam->sizeConcreteBlock(i,1)/2.0 )/1000.0 
              << G4endl;

   G4RotationMatrix* G4RotationOfConcreteBlock = new G4RotationMatrix;
       G4RotationOfConcreteBlock->rotateX( 0.0*deg );
       G4RotationOfConcreteBlock->rotateY( 0.0*deg );
       G4RotationOfConcreteBlock->rotateZ( 0.0*deg );

   G4Transform3D Trans3DConcreteBlock( *G4RotationOfConcreteBlock, G4PositionOfConcreteBlock );

   physiConcreteBlock[i] = new G4PVPlacement( Trans3DConcreteBlock,
                                              logicConcreteBlock[i],
                                              "ConcreteBlock",
                                              logicWorld,
                                              false, 0 );
   }
 
   //-------------------------------------------------
   // Detector Plane
   //-------------------------------------------------
     G4VisAttributes* DPVisAtt = new G4VisAttributes(true,G4Colour(1.0, 0.0, 0.0));  
   DPVisAtt->SetForceWireframe(true);
   // DPAVistt->SetForceSolid(true);
  
  for( int i=0; i < elsparam->numDP() ; i++ ){

    solidDetectorPlane[i] = new G4Box( "DPlane",
                                       elsparam->sizeDetectorPlane(i,0)*0.5*mm,
                                       elsparam->sizeDetectorPlane(i,1)*0.5*mm,
                                       elsparam->sizeDetectorPlane(i,2)*0.5*mm );
    
     G4ThreeVector G4PositionOfDetectorPlane
       = G4ThreeVector( elsparam->positionDetectorPlane(i,0)*mm,
                        elsparam->positionDetectorPlane(i,1)*mm,
			elsparam->positionDetectorPlane(i,2)*mm ); 

     G4RotationMatrix* G4RotationOfDetectorPlane = new G4RotationMatrix;
     G4RotationOfDetectorPlane->rotateX( 0.*deg );
     G4RotationOfDetectorPlane->rotateY( 0.*deg );
     G4RotationOfDetectorPlane->rotateZ( 0.*deg );

     G4Transform3D Trans3DDetectorPlane( *G4RotationOfDetectorPlane, G4PositionOfDetectorPlane );
     
     logicDetectorPlane[i] = new G4LogicalVolume( solidDetectorPlane[i],
						   Air,
					           DP_NAME[i],
						   0,0,0 );
     logicDetectorPlane[i]->SetVisAttributes(DPVisAtt);

     if( Outoputfileflag_detectorplane ){
          physiDetectorPlane[i] = new G4PVPlacement( Trans3DDetectorPlane,
	   				             logicDetectorPlane[i],
					             DP_NAME[i],
					             logicWorld,
					             false, 0 );
     }

  }

  //=================================================
  // Virtual Faraday made by Dmitri-san 
  // created in 2013.04.22 
  // Position is fixed, and geometry is also fixed 
  //-------------------------------------------------
  if( Flag_DmitriFC == 1 ){

  ///////////////////  Copper Faraday cup //////////////////////////////////
  G4double outerRadius = 30 * mm;
  G4double hz = 60 * mm;

  G4ThreeVector G4PositionOfFCDmitri 
     = G4ThreeVector( elsparam->Beam_injection_position(0)*mm,
                      elsparam->Beam_injection_position(1)*mm,
                      10000.0*mm );
  fSolidFaradayCup = new G4Tubs("FaradayCupDmitri", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicFaradayCup = new G4LogicalVolume(fSolidFaradayCup, fCopper, "FaradayCupDmitri");
  fPhysiFaradayCup = new G4PVPlacement( 0, G4PositionOfFCDmitri, 
                                       fLogicFaradayCup, "FaradayCupDmitri", logicWorld, false, 0);

  ////////////////// Upper Titanium Plate //////////////////////////////////////////
  outerRadius = 99.0 * mm;
  hz = 0.15 * mm;
  G4ThreeVector G4PositionOfFCDmitriUpperTi 
    = G4ThreeVector( G4PositionOfFCDmitri.x(), G4PositionOfFCDmitri.y(),
		     G4PositionOfFCDmitri.z() - 3.*cm - 3.*mm );
                     //G4PositionOfFCDmitri.z() - 3.*cm - 6.*mm );
  fSolidUpperTiPlate = new G4Tubs("FCDmitriUTi", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicUpperTiPlate = new G4LogicalVolume(fSolidUpperTiPlate, fTitanium, "FCDmitriUTi");
  fPhysUpperTiPlate = new G4PVPlacement( 0, G4PositionOfFCDmitriUpperTi, 
                                         fLogicUpperTiPlate, "FCDmitriUTi", logicWorld, false, 0);  

  ////////////////// Upper Copper Plate //////////////////////////////////////////
  outerRadius = 99.0 * mm;
  hz = 0.2 * mm;
  G4ThreeVector G4PositionOfFCDmitriUpperCu 
    = G4ThreeVector( G4PositionOfFCDmitri.x(), G4PositionOfFCDmitri.y(),  
                     G4PositionOfFCDmitri.z() - 3.*cm - 6.*mm );
		     //G4PositionOfFCDmitri.z() - 3.*cm - 6.*mm - 7.*mm );
  fSolidUpperCuPlate = new G4Tubs("FCDmitriUCu", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicUpperCuPlate = new G4LogicalVolume(fSolidUpperCuPlate, fCopper, "FCDmitriUCu");
  fPhysUpperCuPlate = new G4PVPlacement( 0, G4PositionOfFCDmitriUpperCu, 
                                         fLogicUpperCuPlate, "FCDmitriUCu", logicWorld, false, 0);

  
  /////////////////// Count Air Plate /////////////////////////////////////////////
  outerRadius = 10.0 * m;
  hz = 0.01 * mm;
  G4ThreeVector G4PositionOfFCdmitriCountAir
    = G4ThreeVector( G4PositionOfFCDmitri.x(), G4PositionOfFCDmitri.y(),
		     G4PositionOfFCDmitri.z() -3.*cm - 6.*mm - 7.*mm - 5.*cm );
  fSolidCountAirPlate = new G4Tubs("FCDmitriCountPlate", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicCountAirPlate = new G4LogicalVolume(fSolidCountAirPlate, Air, "FCDmitriCountPlate");
  fPhysCountAirPlate  = new G4PVPlacement( 0, G4PositionOfFCdmitriCountAir, fLogicCountAirPlate,
   					   "FCDmitriCountPlate", logicWorld, false, 0);

  ////////////////// Bottom Copper Plate //////////////////////////////////////////
  outerRadius = 99.0 * mm;
  //hz = 0.2 * mm;
  hz = 0.1 * mm;
  G4ThreeVector G4PositionOfFCDmitriBottomCu
    = G4ThreeVector( G4PositionOfFCDmitri.x(), G4PositionOfFCDmitri.y(),
		     G4PositionOfFCDmitri.z() -3.*cm - 6.*mm - 10.*cm  );
		     //G4PositionOfFCDmitri.z() -3.*cm - 6.*mm - 7.*mm - 10.*cm );
  fSolidBottomCuPlate = new G4Tubs("FCDmitriBCu", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicBottomCuPlate = new G4LogicalVolume(fSolidBottomCuPlate, fCopper, "FCDmitriBCu");
  fPhysBottomCuPlate = new G4PVPlacement( 0, G4PositionOfFCDmitriBottomCu, fLogicBottomCuPlate,
   					  "FCDmitriBCu", logicWorld, false, 0);

  ////////////////// Bottom Titanium Plate //////////////////////////////////////////
  outerRadius = 99.0 * mm;
  hz = 0.125 * mm;
  G4ThreeVector G4PositionOfFCDmitriBottomTi
    = G4ThreeVector( G4PositionOfFCDmitri.x(), G4PositionOfFCDmitri.y(),
	             G4PositionOfFCDmitri.z() -3.*cm - 6.*mm - 10.*cm - 0.1 * mm - hz / 2 );
		     //G4PositionOfFCDmitri.z() -3.*cm - 6.*mm - 7.*mm - 10.*cm - 0.1 * mm - 10.*mm - hz / 2 );
  fSolidBottomTiPlate = new G4Tubs("FCDmitriBTi", 0.*mm, outerRadius, hz / 2.0, 0.*deg , 360 * deg );
  fLogicBottomTiPlate = new G4LogicalVolume(fSolidBottomTiPlate, fTitanium, "FCDmitriBTi" );
  fPhysBottomTiPlate = new G4PVPlacement( 0, G4PositionOfFCDmitriBottomTi,
					  fLogicBottomTiPlate, "FCDmitriBTi", logicWorld, false, 0);

  }

   
  //=======================================
  // FD Dummy Mirror , added in 2013.09.10
  //=======================================  
  G4double InnerRadiusMirror = 6.067;      // m  
  G4double DepthMirror       = 5.0*0.001; // m : 11mm x 0.001 = 0.011 m
  G4double OuterRadiusMirror = InnerRadiusMirror+DepthMirror; 
  G4double DiameterMirror    = 2.00*2.0;   //1.730*2.0; // m
  G4double DistanceFromBeamToMirror = 99.931;  // m

  solidMirrorSphere = new G4Sphere("MirrorSphere",
									InnerRadiusMirror*m,  OuterRadiusMirror*m,
								    0.0*deg, 360.0*deg, 0.0*deg, 360.0*deg );

  solidMirrorTubs   = new G4Tubs("MirrorTubs", 
								 0.0*mm ,
								 DiameterMirror/2.0*m ,
								 1.0*m,
								 0.*deg, 360.*deg );
 
  G4ThreeVector DummyMirrortransVector = G4ThreeVector( 0.0*m, 0.0*m, (-1.)*InnerRadiusMirror*m );
  G4RotationMatrix* DummyMirrorRot = new G4RotationMatrix;
                    DummyMirrorRot->rotateX( 0.0*deg );
                    DummyMirrorRot->rotateY( 0.0*deg );
                    DummyMirrorRot->rotateZ( 0.0*deg );
  G4Transform3D DummyMirrortrans3D( *DummyMirrorRot, DummyMirrortransVector  );

  solidMirror
      = new G4IntersectionSolid("Mirror", solidMirrorSphere, solidMirrorTubs, DummyMirrortrans3D );

  logicMirror = new G4LogicalVolume( solidMirror, MaterialtPureAl, "Mirror", 0, 0, 0);


  G4ThreeVector G4PositionOfMirror6 = G4ThreeVector(  11.4862*m - DistanceFromBeamToMirror*m + InnerRadiusMirror*cos(25.5*deg)*m, 
                                                     -1.7182*m,
                                                      1.5267*m + InnerRadiusMirror*cos( (90.0-25.5)*deg )*m );
  G4RotationMatrix* G4RotationOfMirror6 = new G4RotationMatrix;
                    G4RotationOfMirror6->rotateZ( 0.4618644*deg );
					G4RotationOfMirror6->rotateY( (90.0-25.5)*deg );
                    G4RotationOfMirror6->rotateX( 0.0*deg );
  G4Transform3D Trans3DMirror6( *G4RotationOfMirror6, G4PositionOfMirror6 ) ;

  G4ThreeVector G4PositionOfMirror7 = G4ThreeVector(  11.4862*m - DistanceFromBeamToMirror*m + InnerRadiusMirror*cos(10.5*deg)*m,
                                                      -1.7182*m,
                                                      5.4994*m + InnerRadiusMirror*cos( (90.0-10.5)*deg )*m );
  G4RotationMatrix* G4RotationOfMirror7 = new G4RotationMatrix;
                    G4RotationOfMirror7->rotateZ( 0.4618644*deg );
                    G4RotationOfMirror7->rotateY( (90.0-10.5)*deg );
                    G4RotationOfMirror7->rotateX( 0.0*deg );
  G4Transform3D Trans3DMirror7( *G4RotationOfMirror7, G4PositionOfMirror7 );

  if( Cerenkov_flag == 1 ){
      physiMirror6 = new G4PVPlacement( Trans3DMirror6, logicMirror, "Mirror6", logicWorld, false, 0 );
      physiMirror7 = new G4PVPlacement( Trans3DMirror7, logicMirror, "Mirror7", logicWorld, false, 0 );
  }
  //==================================

  //==================================
  // ICE definition -- added 23/1/14
  //==================================
  solidIce = new G4Box("ICE",
		       elsparam->sizeIce(0)*0.5*mm,
		       elsparam->sizeIce(1)*0.5*mm,
		       elsparam->sizeIce(2)*0.5*mm );

  logicIce = new G4LogicalVolume( solidIce, ICE, "ICE", 0, 0, 0);

  G4ThreeVector G4PositionOfIce = G4ThreeVector( elsparam->positionIce(0)*mm,
						 elsparam->positionIce(1)*mm,
						 elsparam->positionIce(2)*mm);

  G4RotationMatrix* G4RotationOfIce = new G4RotationMatrix;
  G4RotationOfIce->rotateX( elsparam->rotationIce(0)*deg );
  G4RotationOfIce->rotateY( elsparam->rotationIce(1)*deg );
  G4RotationOfIce->rotateZ( elsparam->rotationIce(2)*deg );

  G4Transform3D Trans3DIce( *G4RotationOfIce, G4PositionOfIce );

  physiIce = new G4PVPlacement( Trans3DIce, logicIce, "ICE", logicWorld, false, 0 );

  

  return physiWorld;
}

void DetectorConstruction::setTargetMaterial(G4String materialName)
{
  // search the material by its name
  G4Material* pttoMaterial = G4Material::GetMaterial(materialName);
  if (pttoMaterial)
  {
    logicConcretepad->SetMaterial(pttoMaterial);   

    logicELSContainer->SetMaterial(pttoMaterial);

    for( int i=0; i<3; i++ ) logicBase[i]->SetMaterial(pttoMaterial);

    for( int i=0; i<elsparam->numBLcomp(); i++ ){
      logicBL_tubs[i]->SetMaterial(pttoMaterial);
      logicBL_vac_tubs[i]->SetMaterial(pttoMaterial);
      logicBL_FF_tubs[i]->SetMaterial(pttoMaterial);
      logicBL_FF_vac_tubs[i]->SetMaterial(pttoMaterial);
      logicBL_BF_tubs[i]->SetMaterial(pttoMaterial);
      logicBL_BF_vac_tubs[i]->SetMaterial(pttoMaterial);
    }      
    logicEGUNBlankFrange->SetMaterial(pttoMaterial);
    logicAnode->SetMaterial(pttoMaterial);

    for( int i=0; i<2; i++ ) logicSM_FF_tubs[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_FF_vac_tubs[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_BF_tubs[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_BF_vac_tubs[i]->SetMaterial(pttoMaterial);

    for( int i=0; i<2; i++ ) logicSM_duct[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_duct_vac[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_body[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicSM_body_vac[i]->SetMaterial(pttoMaterial);

    logicCuCollimeter->SetMaterial(pttoMaterial);
    logicCuCollimeter_vacuum->SetMaterial(pttoMaterial);

    logicBMF->SetMaterial(pttoMaterial);
    logicBMF_vac->SetMaterial(pttoMaterial);
    for( int i=0; i<3; i++ ) logicBMLduct[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<3; i++ ) logicBMLduct_vac[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicBMmain[i]->SetMaterial(pttoMaterial);
    logicBMmain_vac->SetMaterial(pttoMaterial);
    logicBMSub1->SetMaterial(pttoMaterial);
    logicBMSub2->SetMaterial(pttoMaterial);
    logicBMSub3->SetMaterial(pttoMaterial);
    logicBMSub4->SetMaterial(pttoMaterial);

    for( int i=0; i<5; i++ ) logicBMYork[i]->SetMaterial(pttoMaterial);
    logicBMCoil->SetMaterial(pttoMaterial);

    logicSLIT_FF_tubs->SetMaterial(pttoMaterial);
    logicSLIT_FF_vac_tubs->SetMaterial(pttoMaterial);
    logicSLIT_BF_tubs->SetMaterial(pttoMaterial);   
    logicSLIT_BF_vac_tubs->SetMaterial(pttoMaterial);   
    logicSLIT_duct->SetMaterial(pttoMaterial);
    logicSLIT_body->SetMaterial(pttoMaterial);
    logicSLIT_duct_vac->SetMaterial(pttoMaterial);
    logicSLIT_body_vac->SetMaterial(pttoMaterial);
    for( int i=0; i<2; i++ ) logicCollimeter[i]->SetMaterial(pttoMaterial);

    logicBlankFrange->SetMaterial(pttoMaterial);
    logicWindow->SetMaterial(pttoMaterial);
    logicWindow_vac->SetMaterial(pttoMaterial);
    logicTiWindow->SetMaterial(pttoMaterial);

    for( int i=0; i<3; i++ ) logicFDAlCylinder[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<4; i++ ) logicFDBody[i]->SetMaterial(pttoMaterial);

    logicGNDLayerFC4->SetMaterial(pttoMaterial);
    logicShieldLayerFC4->SetMaterial(pttoMaterial);
    logicBodyFC4->SetMaterial(pttoMaterial);
    logicIso1FC4->SetMaterial(pttoMaterial);
    logicIso2FC4->SetMaterial(pttoMaterial);
    logicTopPlateFC4->SetMaterial(pttoMaterial);

    // B.K.Shin-san7s FC5 



    logicSM3->SetMaterial(pttoMaterial);

    logicBeamAttenuator->SetMaterial(pttoMaterial);

    logicPbCollimator->SetMaterial(pttoMaterial);

    for( int i=0; i<5; i++ ) logicVirtualChamberCylinder[i]->SetMaterial(pttoMaterial);
    for( int i=0; i<4; i++ ) logicVirtualChamberBody[i]->SetMaterial(pttoMaterial);
    logicVirtualChamberTarget->SetMaterial(pttoMaterial);
    logicVirtualChamberCylinder_vac->SetMaterial(pttoMaterial);        
   
    for( int i=0; i<5; i++ )  logicPbTower[i]->SetMaterial(pttoMaterial); 
    for( int i=0; i<2; i++ ){ logicPbBlockBMparts1[i]->SetMaterial(pttoMaterial);
                              logicPbBlockBMparts2[i]->SetMaterial(pttoMaterial);
                              logicPbBlockBMparts3[i]->SetMaterial(pttoMaterial);
    }

    for( int i=0; i<20; i++ ) logicConcreteBlock[i]->SetMaterial(pttoMaterial);    

    for( int i=0; i<200; i++ ) logicDetectorPlane[i]->SetMaterial(pttoMaterial);


    fLogicFaradayCup->SetMaterial(pttoMaterial); 
    fLogicUpperTiPlate->SetMaterial(pttoMaterial);
    fLogicUpperCuPlate->SetMaterial(pttoMaterial);
    fLogicCountAirPlate->SetMaterial(pttoMaterial);    
    fLogicBottomCuPlate->SetMaterial(pttoMaterial);
    fLogicBottomTiPlate->SetMaterial(pttoMaterial);

    logicMirror->SetMaterial(pttoMaterial);

  }

}

      
  
