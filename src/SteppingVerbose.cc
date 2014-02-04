//===================================================================
//      SteppingVerbose.cc
//      Author    : T.Shibata
//      Created     : 2010.02.11
//      Last Update : 2010.11.08
//      Last Update : 2011.01.--
//      Last Update : 2011.02.02
//      Last Update : 2011.06.11
//      Last Update : 2012.10.02
//      Last Update : 2012.10.10
//===================================================================

//#include <iomanip.h>
#include <iomanip>

#include "SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include "DetectorPlaneName.hh"

//------------------------------------------------------------------------//
// Constructor
//------------------------------------------------------------------------//
SteppingVerbose::SteppingVerbose( Plist *plist, RunOption runOpt ){

  // Create ELSParam object and initialize geometry
  elsparam = new ELSParameters();
  elsparam->els_geometory_initialize(plist);

  // I think this is their detector.
  // maybe DELETE
  fdparam = new FDParameters();

  // Energy Deposit Initialize //----------
  // I want to use this but need to UPDATE
  de_mesh = new dE_Mesh();

  de_mesh->set_delta_PathLength( plist->delta_PathLength() ); // unit=m

  position_of_beaminjection.setX(elsparam->Beam_injection_position(0)/m);
  position_of_beaminjection.setY(elsparam->Beam_injection_position(1)/m);
  position_of_beaminjection.setZ(elsparam->Beam_injection_position(2)/m);
  de_mesh->set_position_of_beaminjection(position_of_beaminjection);
 
  position_of_mirror_6_BRM_coordinate.setX(plist->Position_of_mirror_6_BRM_coordinate(0));
  position_of_mirror_6_BRM_coordinate.setY(plist->Position_of_mirror_6_BRM_coordinate(1));
  position_of_mirror_6_BRM_coordinate.setZ(plist->Position_of_mirror_6_BRM_coordinate(2));
  de_mesh->set_position_of_mirror_6_BRM_coordinate( position_of_mirror_6_BRM_coordinate );

  de_mesh->set_angle_ELS_mirror_6( plist->Angle_ELS_mirror_6() );     
  de_mesh->set_distance_ELS_from_mirror6( plist->Distance_ELS_from_mirror6() );

  de_mesh->set_angle_North_Yaxis_BRM( plist->Angle_North_Yaxis_BRM() );

  mesh_dv[0]=plist->Mesh_dX();  
  mesh_dv[1]=plist->Mesh_dY();  
  mesh_dv[2]=plist->Mesh_dZ();  
  de_mesh->set_mesh_dV( mesh_dv );  
 
  mesh_vregion[0]=plist->Mesh_Xregion();
  mesh_vregion[1]=plist->Mesh_Yregion();
  mesh_vregion[2]=plist->Mesh_Zregion();
  de_mesh->set_mesh_Vregion( mesh_vregion );

  de_mesh->set_mesh_number();
  de_mesh->set_ELS_origin_BRM_coordinate();
  //-----------------------------------------

  // I don't think I need all of these files
  // Consider removing and condensing this code
  // UPDATE
  m_runOpt = runOpt;
  if( runOpt & RO_EnergyDump ) EnergyDumpFile.open("ELS_output/energyDump.dat");
  if( runOpt & RO_QuickCheck ) QuickDumpFile.open("ELS_output/quickDump.dat");

  // Set some beam constants
  CenterofBeamInjection[0] = 11486.2;  // unit=mm
  CenterofBeamInjection[1] = -1718.2;  // unit=mm
  CenterofBeamInjection[2] =  2407.0;  // unit=mm
  
  DirectionofBeamInjection[0] = 0. ;
  DirectionofBeamInjection[1] = 0. ;
  DirectionofBeamInjection[2] = 1. ;

  ChargeDeposit_in_FaradayCup=0;


  m_hackEvent = -1;
}

//------------------------------------------------------------------------//
// Destructor
//------------------------------------------------------------------------//
SteppingVerbose::~SteppingVerbose()
{
  delete elsparam;
  delete fdparam;
  delete de_mesh; 

  if(EnergyDumpFile.is_open()) EnergyDumpFile.close();
  if(QuickDumpFile.is_open())  QuickDumpFile.close();

}

//------------------------------------------------------------------------//
// Step info.  This is where we can dump the information at each 
// individual step.
//------------------------------------------------------------------------//
void SteppingVerbose::StepInfo()
{

  CopyState();

  //if( fTrack->GetKineticEnergy()/MeV < 100 )
  //if( fTrack->GetKineticEnergy()/MeV < 0.611 )
  //fTrack->SetTrackStatus(fStopAndKill);

  //if( m_runOpt & RO_EnergyDump ) energyDump();
  if( m_runOpt & RO_QuickCheck ) quickCheck();

  StepCount++;

  /*
  TrackID             : fTrack->GetTrackID()
  ParentID            : fTrack->GetParentID()
  Particle Mass       : Mass
  Current Step Number : fTrack->GetCurrentStepNumber()
  Kinetic Energy      : fTrack->GetKineticEnergy()
  Energy Diposit      : fStep->GetTotalEnergyDeposit()
  Step Length         : fStep->GetStepLength()
  Track Length        : fTrack->GetTrackLength()
  Volume Name         : fTrack->GetVolume()->GetName()
  Process Name        : fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
  */

  //-------------------------------
  // dE/dX information
  //-------------------------------
  // Output for Checking dE/dX of one ionization process ( checking the formula of ionization energy losses )
  // Modified in 2011.07.19-20
  /*
  if( OutputFile2_flag ){
    if( ParticleID == -1 && ParticleTrackID == 1 ){  // only primary electron
      if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "eIoni" ){
           OutputFile2 << " "          
	     //	       << setprecision(3) << ParticleTrackID                << " "
	     //	       << setprecision(3) << ParticleID                     << " "
	               << setprecision(6) << fTrack->GetKineticEnergy()     << " " 
	               << setprecision(6) << fStep->GetTotalEnergyDeposit() << " " 
                       << setprecision(6) << fStep->GetStepLength()/cm      << " "
           
	     //        << setprecision(6) << fTrack->GetPosition().x()/m    << " "
	     //	       << setprecision(6) << fTrack->GetPosition().y()/m    << " "
	     //	       << setprecision(6) << fTrack->GetPosition().z()/m    << " "
                       << G4endl;
      }
   }
  }
  */

}

//------------------------------------------------------------------------//
// Energy deposit
//------------------------------------------------------------------------//
void SteppingVerbose::energyDump()
{
  
  //---------------------------------
  // Energy Loss Output
  //---------------------------------
  
  int ProcessID(-100);
  if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
      == "Transportation" ) ProcessID=-1;
  /* Transportation has energy deposit ?? */
  /* Check */
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
	   == "msc" ) ProcessID=0;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() 
	   == "eIoni" ) ProcessID=1;      
  /* Ionization by primary and secondary electrons 
     The energy deposit mean total kinetic energy of secondary electron its 
     energy less than cut energy ( = 1keV )
  */
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	   == "eBrem" ) ProcessID=2;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	   == "annihil" ) ProcessID=3;
  
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	   == "phot" ) ProcessID=4;
  /* Photo-electric effect 
     The energy deposit mean total kinetic energy of photo-electrons its
     energy less than cut energy ( = 1keV )
  */
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	   == "compt" ) ProcessID=5;
  else if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	   == "conv" ) ProcessID=6;
  else ProcessID=-100;
  
  // The Energy deposit is used for yield of fluorescence photon
  if( fStep->GetTotalEnergyDeposit() > 0. &&
      fStep->GetPostStepPoint()->GetPosition().z() > 0. ){         
    if( fTrack->GetVolume()->GetName()=="World" && fTrack->GetNextVolume() != 0 ){
      
      CreatedEnergyDepositMeth created_energy_deposit_meth;
      if( de_mesh->energy_deposit_mesh_generator( ProcessID,
						  fStep->GetTotalEnergyDeposit()/MeV,
						  fStep->GetPreStepPoint()->GetPosition()/m,
						  fStep->GetPostStepPoint()->GetPosition()/m,
						  created_energy_deposit_meth ) ){
	
	for( CreatedEnergyDepositMeth::const_iterator it = created_energy_deposit_meth.begin();
	     it!=created_energy_deposit_meth.end(); it++ ){
	  
	  //-----------------------------------------------
	  // 2011.01.29 T.Shibata
	  //  Check the Field of View ( Design Geometory )
	  //  Jedge of Camera #6 and #7
	  CLHEP::Hep3Vector ShiftMirror = CLHEP::Hep3Vector( 0., 0., 0. );
	  CLHEP::Hep3Vector VertexPoint = CLHEP::Hep3Vector( (*it).vertex_brmcoord.x()*(-1),
							     (*it).vertex_brmcoord.y()*(-1),
							     (*it).vertex_brmcoord.z()      );
	  int FoV_judge[12];
	  for( int i=0; i<12; i++ ){
	    FoV_judge[i]=0;
	    if( fdparam->FieldOfViewofCamera(i, VertexPoint, ShiftMirror ) ){
	      FoV_judge[i]=1;
	    }
	  }
	  
	  //-----------------------------------------------
	  EnergyDumpFile
	    // unit = MeV --> eV
	    << setprecision(10) << (*it).energy_deposit/eV         << " "
	    // unit = cm
	    << setprecision(15) << (*it).vertex_eathcoord.x()*100. << " "  // = cm
	    << setprecision(15) << (*it).vertex_eathcoord.y()*100. << " "  // = cm
	    << setprecision(15) << (*it).vertex_eathcoord.z()*100. << " "  // = cm
	    << ProcessID        << " " 
	    << ParticleID       << " " 
	    << ParticleTrackID  << " "                       
	    << fTrack->GetTrackID() << " "
	    << FoV_judge[6]     << " " 
	    << FoV_judge[7]     << " " 
	    << setprecision(10) << PrimarElevtronEnergy    << " " 
	    << setprecision(6)  << fTrack->GetGlobalTime() << " "  // unit=nsed : added in 2012.04.12 
	    << G4endl;   
	  
	}// end loop over energy deposit meth
	
      }// end if energy deposit mesh is ok
      else{
	fprintf(stderr, "energy_deposit_mesh_generator has some problem...\n");
      }
      // clean up energy method
      created_energy_deposit_meth.clear();
	  
    }// end if volume is in world and next volume isn't outside world 

    else{   
      // The Energy deposit can not be used for yield of fluorescence photon...
      //  ... This volume is not World.
      EnergyDumpFile
	<< setprecision(15) << fStep->GetTotalEnergyDeposit()/eV << " "
	<< setprecision(15) << fStep->GetPostStepPoint()->GetPosition().x()/cm  << " "
	<< setprecision(15) << fStep->GetPostStepPoint()->GetPosition().y()/cm  << " "
	<< setprecision(15) << fStep->GetPostStepPoint()->GetPosition().z()/cm  << " "
	<< "-1113"          << " "
	<< ParticleID       << " "
	<< ParticleTrackID  << " " 
	<< -1               << " " 
	<< -1               << " " 
	<< setprecision(10) << PrimarElevtronEnergy << " " 
	<< setprecision(6)  << fTrack->GetGlobalTime() << " "  // unit=nsed : added in 2012.04.12 
	<< G4endl;
      
      
    }// end if not in world
	
  }// end if z > 0
  else if( fStep->GetTotalEnergyDeposit() > 0. &&
	   fStep->GetPostStepPoint()->GetPosition().z() < 0. ){
    
    
    // The Energy deposit can not be used for yield of fluorescence photon...                             
    // The height is under-ground.  
    EnergyDumpFile
      << setprecision(15) << fStep->GetTotalEnergyDeposit()/eV << " "
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().x()/cm  << " "
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().y()/cm  << " "
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().z()/cm  << " "
      << "-1112"          << " "
      << ParticleID       << " "
      << ParticleTrackID  << " " 
      << -1               << " " 
      << -1               << " " 
      << setprecision(10) << PrimarElevtronEnergy << " " 
      << setprecision(6)  << fTrack->GetGlobalTime() << " "  // unit=nsed : added in 2012.04.12 
      << G4endl;
    
  }
  
  if( fTrack->GetVolume()->GetName()=="World" && fTrack->GetNextVolume() == 0 ){
    EnergyDumpFile
      << setprecision(15) << fTrack->GetKineticEnergy()/eV << " " 
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().x()/cm  << " "
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().y()/cm  << " "
      << setprecision(15) << fStep->GetPostStepPoint()->GetPosition().z()/cm  << " "
      << "-1111"          << " "  
      << ParticleID       << " "
      << ParticleTrackID  << " " 
      << -1               << " " 
      << -1               << " " 
      << setprecision(10) << PrimarElevtronEnergy << " " 
      << setprecision(6)  << fTrack->GetGlobalTime() << " "  // unit=nsed : added in 2012.04.12 
      << G4endl;
    
  }
  
}
  
  
//------------------------------------------------------------------------//
// My silly quick method
//------------------------------------------------------------------------//
void SteppingVerbose::quickCheck()
{
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
  // This will be the output I want to analyze.  I am 
  // not really following what the other output is 
  // dumping, so try to build it up myself
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//


  // New output code
  // First let's try to construct the shower shape
  // Only consider e+ and e-
  //if( ParticleID == 1 || ParticleID == -1 ){
  if( fTrack->GetVolume()->GetName() == "ICE"){

    QuickDumpFile <<
      fStep->GetPostStepPoint()->GetPosition().x()/cm <<" "<<
      fStep->GetPostStepPoint()->GetPosition().y()/cm <<" "<<
      fStep->GetPostStepPoint()->GetPosition().z()/cm <<" "<<
      fStep->GetStepLength()/cm <<" "<<
      fTrack->GetKineticEnergy()/MeV<<" "<<
      fStep->GetTotalEnergyDeposit()/MeV<<" "<<
      fTrack->GetParticleDefinition()->GetPDGEncoding()<<" "<<
      fTrack->GetTrackID()<<" "<<
      fTrack->GetParentID()<<" "<<
      CreatedProcess<<" "<<
      G4endl;	       
    }// end check
  
  
}

//------------------------------------------------------------------------//
// Here we start the tracking algorithm. It seems this is called
// before the stepping takes place... maybe it is called at each step
// This should be CHECKED
//------------------------------------------------------------------------//
void SteppingVerbose::TrackingStarted()
{

  // I will leave this alone for now.

  CopyState();
  
  // Particle ID
  G4String ParticleName = fTrack->GetDefinition()->GetParticleName();
  
  if( ParticleName == "gamma" ){            ParticleID =  0; }     
  else if( ParticleName == "e-" ){          ParticleID = -1; } 
  else if( ParticleName == "e+" ){          ParticleID =  1; } 
  else if( ParticleName == "proton"){       ParticleID =  2; }
  else if( ParticleName == "anti_proton"){  ParticleID = -2; }
  else if( ParticleName == "neutron"){      ParticleID =  3; }
  else if( ParticleName == "anti_neutron"){ ParticleID = -3; }  
  else if( ParticleName == "alpha" ){       ParticleID =  4; } 
  else if( ParticleName == "deuteron" ){    ParticleID =  5; }
  
  else if( ParticleName == "opticalphoton" ){ 
    ParticleID = 50; 
  } 
  
  else if( ParticleName == "geantino" ){ ParticleID = -100; }
  
  else{ ParticleID = 100; }
  
  // Track ID : Primary particle is 1, Secondary Particle is 2;
  if( fTrack->GetTrackID() == 1 ){ 
    ParticleTrackID = 1; 
  }
  else{
    ParticleTrackID = 2;
  }
  
  // TrackPID[ fTrack->GetTrackID() ] = ParticleID; 
  
  for( int i=0; i<200; i++ ) ParticleTrackID_FLAG[i] = 0;
  
  // Initial Energy & Mometum & Position  
  InitialEnergy   = fTrack->GetKineticEnergy();
  InitialMomentum = fTrack->GetMomentum(); //fTrack->GetMomentum().mag();
  InitialPosition = fTrack->GetPosition();
  
  if( ParticleID== -1 && ParticleTrackID == 1 ) PrimarElevtronEnergy = InitialEnergy/eV;;
  
  InitialVolumeID = -1;
  StoppedVolumeID = -1;
  
  //if( ParticleTrackID == 2 &&
  //abs(ParticleID) == 1 &&
  if( InitialPosition.z()/m > 4 ){
    
    /*
      OutputFile100 
      << ParticleName << " " 
      << fTrack->GetCreatorProcess()->GetProcessName() << " "
      //             << ParticleID   << " " 
      //             << fTrack->GetTrackID()      << " "
      // 	      << fTrack->GetParentID()     << " "
      << InitialEnergy/eV          << " "
      << InitialPosition.z()/m     << " "
      //             << fTrack->GetMomentum().x()/fTrack->GetMomentum().mag() << " " 
      //	      << fTrack->GetMomentum().y()/fTrack->GetMomentum().mag() << " " 
      //	      << fTrack->GetMomentum().z()/fTrack->GetMomentum().mag() << " " 
      << G4endl;
    */
    
    /*if( fTrack->GetVolume()->GetName() == "FC4Body" ) InitialVolumeID = 1;
    if( fTrack->GetVolume()->GetName() == "FC4GNDLayer" ) InitialVolumeID = 2;
    if( fTrack->GetVolume()->GetName() == "FC4ShieldLayer" ) InitialVolumeID = 3;
    if( fTrack->GetVolume()->GetName() == "FC4TopPlate" ) InitialVolumeID = 4;
    if( fTrack->GetVolume()->GetName() == "FC4Iso1" ) InitialVolumeID = 5;
    if( fTrack->GetVolume()->GetName() == "FC4Iso2" ) InitialVolumeID = 6;*/
    
    if( fTrack->GetCreatorProcess()->GetProcessName()=="msc" )   CreatedProcess = 0;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="eIoni" ) CreatedProcess = 1;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="eBrem" ) CreatedProcess = 2;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="phot" )  CreatedProcess = 3;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="compt" ) CreatedProcess = 4;
    if( fTrack->GetCreatorProcess()->GetProcessName()=="conv" )  CreatedProcess = 5;
    
  }
  
  StepCount=0;
  
  //--------------------------------------------------------
  // Kill paritcles except gamma,e+-,opticalphoton,geantino 
  //--------------------------------------------------------
  // Why is this needed??  We don't want to simulate hadronic shower in here?
  if( ParticleID == 2 || ParticleID == -2 ||
      ParticleID == 3 || ParticleID == -3 ||
      ParticleID == 4 || ParticleID == -4 ||
      ParticleID == 5 || ParticleID == -5 ||
      ParticleID == 100 ) {
    fTrack->SetTrackStatus(fStopAndKill);
  }

  //if( fTrack->GetKineticEnergy()/MeV < 100 )
  //if( fTrack->GetKineticEnergy()/MeV < 0.611 )
  //fTrack->SetTrackStatus(fStopAndKill);
  
  if( fTrack->GetTrackID() == 1 ){
    if( m_hackEvent >= 0 )
      QuickDumpFile << "End" << G4endl;
    m_hackEvent++;
    QuickDumpFile << "Event: " << m_hackEvent << G4endl;
  }
  //G4cout.precision(prec);
}


