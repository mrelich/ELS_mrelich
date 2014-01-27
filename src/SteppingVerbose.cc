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

SteppingVerbose::SteppingVerbose( Plist *plist,
                                  char* Fname0,
                                  char* Fname1,
                                  char* Fname2,
                                  char* Fname3,
                                  char* Fname4,
				  char* mattOut){

  elsparam = new ELSParameters();
  elsparam->els_geometory_initialize(plist);

  fdparam = new FDParameters();

  // Energy Deposit Initialize //----------
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

  OutputFile0.open(Fname0);
  OutputFile1.open(Fname1);
  OutputFile2.open(Fname2);
  OutputFile3.open(Fname3);
  OutputFile4.open(Fname4);
  MattOut.open(mattOut);

  OutputFile4<<"Energy "
	     <<"X-pos "
	     <<"Y-pos "
	     <<"Z-pos "
	     <<"procID "
	     <<"partID "
	     <<"trkID "
	     <<" judge6 "
	     <<" judge7 "
	     <<"Primary E "
	     <<"Time"
	     <<endl;
  OutputFile4<<endl;

  OutputFile5.open("screenmonitor3.dat");
  OutputFile6.open("cerenkov.dat");
  OutputFile100.open("debug-output.dat");

  OutputFile0_flag=plist->Outoputfileflag_detectorplane();
  OutputFile1_flag=plist->Outoputfileflag_secondparticle();
  OutputFile2_flag=plist->Outoputfileflag_dedx();

  if( plist->Outoputfileflag_faradaycup()  == 1  || 
      plist->Outoputfileflag_faradaycup()  == 4  || 
      plist->Outoputfileflag_faradaycup()  == 5  || 
      plist->Outoputfileflag_faradaycup4() == 1  || 
      plist->Virtual_chamber_flag() == 1 
     ){ 
      OutputFile3_flag=1;
  }else{ OutputFile3_flag=0; }

  FaradayCup_version=plist->Outoputfileflag_faradaycup();
  //OutputFile3_flag=1;   

  OutputFile4_flag=plist->Outoputfileflag_energydeposit();
  OutputFile5_flag=plist->Screen_monitor3_flag();
  OutputFile6_flag=plist->Cerenkov_process_flag();

  CenterofBeamInjection[0] = 11486.2;  // unit=mm
  CenterofBeamInjection[1] = -1718.2;  // unit=mm
  CenterofBeamInjection[2] =  2407.0;  // unit=mm
  
  DirectionofBeamInjection[0] = 0. ;
  DirectionofBeamInjection[1] = 0. ;
  DirectionofBeamInjection[2] = 1. ;

  ChargeDeposit_in_FaradayCup=0;

  /*
  for( int i=0; i<1000000; i++ ) {
    TrackPID[i] = -9999;
    ParentTrackID[i] = -9999;   
  }
  TrackPID[0]      = 0;
  ParentTrackID[0] = 0;
  */

}
SteppingVerbose::~SteppingVerbose()
{
  delete elsparam;
  delete fdparam;
  delete de_mesh; 

  OutputFile0.close();
  OutputFile1.close();
  OutputFile2.close();
  OutputFile3.close();
  OutputFile4.close();

  OutputFile5.close();
  OutputFile6.close();
  OutputFile100.close();
  MattOut.close();

}

void SteppingVerbose::StepInfo_OLD()
{
}

void SteppingVerbose::StepInfo()
{
  CopyState();
  G4int prec = G4cout.precision(3);

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

  /* ====================================
     This is very temporary Kill Command 
    ====================================== */
  //if( fTrack->GetPosition().z()/m > 4.0 ){ fTrack->SetTrackStatus(fStopAndKill); }

  //-------------------------------
  // Detector Plane information
  //  Scattering Angle 
  //  Energy Distribution
  //-------------------------------
  if( OutputFile0_flag ){

   if( fStepStatus == fGeomBoundary && fTrack->GetVolume()->GetName()!="World" ){
     if( abs(ParticleID) == 1 || abs(ParticleID) == 0 || abs(ParticleID) == 50 ){  // electron or positrons or gamma
     //if( abs(ParticleID) == 1 ){  // electron or positrons
     //if( abs(ParticleID) == 1 && ParticleTrackID == 1 ){  // primary electron

       for( int i=0; i<100; i++ ){

       if( fTrack->GetVolume()->GetName() == DP_NAME[i] ){
       DetectorPlaneID = i;
       if( ParticleTrackID_FLAG[DetectorPlaneID]!=1 ){ // Avoid Double Counts
           ParticleTrackID_FLAG[DetectorPlaneID]=1;

         DetectorPlaneEnergy   = fTrack->GetKineticEnergy();
         DetectorPlaneMomentum = fTrack->GetMomentum(); ///fTrack->GetMomentum().mag();
         DetectorPlanePosition = fTrack->GetPosition();

         DetectorPlaneScatteringAngle = atan(DetectorPlaneMomentum.y()/DetectorPlaneMomentum.z());

	 //if( abs(ParticleID) == 0 && DetectorPlaneEnergy/MeV < 1 ) break;

   	    // Position on the Detector Plane plotting
            OutputFile0 << " " 
               			<< setprecision(3) << ParticleID                << " "  // e+-/gammma             
	                    << setprecision(5) << ParticleTrackID           << " "  // primary/secondary
	    	            << setprecision(4) << DetectorPlaneID           << " "  // detector ID
		   	          //<< setprecision(6) << DetectorPlaneEnergy/eV    << " "  // energy
			            << setprecision(6) << DetectorPlaneEnergy/MeV   << " "  // energy 
                        << setprecision(6) << DetectorPlanePosition.x()/m - position_of_beaminjection.x() << " " 
                        << setprecision(6) << DetectorPlanePosition.y()/m - position_of_beaminjection.y() << " "
                        << setprecision(6) << DetectorPlanePosition.z()/m << " "  
	        << G4endl;                           
          
	    /*
            OutputFile0 << " " 
	                << setprecision(3) << ParticleID                << " "  // e+-/gammma
	                << setprecision(5) << ParticleTrackID           << " "  // primary/secondary
	    	        << setprecision(4) << DetectorPlaneID           << " "  // detector plane
	                << setprecision(6) << DetectorPlaneEnergy/MeV   << " "  // energy
                                                                                // Position(x,y,z)
                        << setprecision(6) << DetectorPlanePosition.x()/m - position_of_beaminjection.x() << " " 
                        << setprecision(6) << DetectorPlanePosition.y()/m - position_of_beaminjection.y() << " "
                        << setprecision(6) << DetectorPlanePosition.z() << " "  
	          	<< setprecision(6) << DetectorPlaneScatteringAngle << " " // Scattering angle
	                << setprecision(5) << TrackPID[ fTrack->GetParentID()] << " "
	                << G4endl;                           
	    */

       }
       }
     }

    }
   }
  }

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

  // Output for Checking dE/dX of total ionization energy losses ( not one process, but all of ionization )
  // Added in 2011.07.20
  if( OutputFile2_flag )
  {
    for( int i=0; i<1; i++ ){
      DetectorPlaneID = -1;
    if( fTrack->GetVolume()->GetName() == DP_NAME[i] ) DetectorPlaneID = 0;
    if( fStepStatus != fGeomBoundary ){
      if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "eIoni" ||
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "msc"   || 
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "eBrem" ||
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "annihil" ||
	  
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "phot"  ||
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "compt" ||
	  fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() == "conv" )
	{  
	  OutputFile2 << " "
		      << setprecision(3) << ParticleID                << " "  // e+-/gammma
		      << setprecision(5) << ParticleTrackID           << " "  // primary/secondary
		      << setprecision(4) << DetectorPlaneID           << " "  // detector plane
		      << setprecision(6) << fTrack->GetKineticEnergy()/MeV     << " "   
		      << setprecision(6) << fStep->GetTotalEnergyDeposit()/MeV << " " 
		      << setprecision(6) << fStep->GetStepLength()/cm          << " " 
		      << setprecision(6) << fTrack->GetPosition().x()/cm    << " "
		      << setprecision(6) << fTrack->GetPosition().y()/cm    << " "
		      << setprecision(6) << fTrack->GetPosition().z()/cm    << " "
		      << G4endl;
	}else{
	
      }
      
    }
    }
  } 
  
  //---------------------------------
  // Faraday Cup Charge Measurement
  //---------------------------------
  if( OutputFile3_flag ){

	if( fStep->GetPostStepPoint()->GetPosition().z()/m > 5.0 ){ fTrack->SetTrackStatus(fStopAndKill); }

    int Volume_id(-1);
   ChargeDeposit_in_FaradayCup=0;  // added in 2012.10.03

   if( abs(ParticleID) == 1 ){

     if(  fTrack->GetKineticEnergy() == 0 || fTrack->GetNextVolume() == 0 ||
	   ( fStepStatus == fGeomBoundary && 
             fTrack->GetVolume()->GetName()!="World" &&
             fTrack->GetVolume()->GetName() == "VirtualChamberCylinder"
	     //fTrack->GetVolume()->GetName() == "VirtualChamberCylinderFront" 
	     )

        ){ //
 
        if( fTrack->GetVolume()->GetName() == "FCBody"    ||
	    fTrack->GetVolume()->GetName() == "FC4Body"   || 
            fTrack->GetVolume()->GetName() == "FC5Body"   ||
            fTrack->GetVolume()->GetName() == "FaradayCupDmitri" ||
            fTrack->GetVolume()->GetName() == "VirtualChamberTarget" ) StoppedVolume = 1;

        if( fTrack->GetVolume()->GetName() == "FC4ShieldLayer"  ) StoppedVolume = 2; 

        // Dmitri-san's geometry study
        if( fTrack->GetVolume()->GetName() == "FCDmitriUTi" ) StoppedVolume = 2;
        if( fTrack->GetVolume()->GetName() == "FCDmitriUCu" ) StoppedVolume = 3;
        if( fTrack->GetVolume()->GetName() == "FCDmitriBCu" ) StoppedVolume = 4;
        if( fTrack->GetVolume()->GetName() == "FCDmitriBTi" ) StoppedVolume = 5;

        double PathLength = (fStep->GetPreStepPoint()->GetPosition() - fStep->GetPostStepPoint()->GetPosition() ).mag();

	if( ( InitialVolume!=StoppedVolume || ParticleTrackID == 1 ) && PathLength > 0. ){

	    // charged particle inject ... add charge ( = particle ID ) 
            if( InitialVolume != 1 && StoppedVolume ==  1 ) ChargeDeposit_in_FaradayCup =      ParticleID;
	    // charged particle reject ... substact charge ( = -1 * particle ID ) 
            if( InitialVolume == 1 && StoppedVolume !=  1 ) ChargeDeposit_in_FaradayCup = (-1)*ParticleID;

            if( ChargeDeposit_in_FaradayCup!=0 ){

   	        CLHEP::Hep3Vector VertexPosition = fStep->GetPostStepPoint()->GetPosition();

                OutputFile3 
			    << setprecision(3) << InitialVolume << " " 
                << setprecision(3) << StoppedVolume << " "  
                << setprecision(3) << ParticleID << " "
                << setprecision(3) << ParticleTrackID << " " 
                << setprecision(3) << ChargeDeposit_in_FaradayCup << " "  // 5 
   			    << setprecision(6) << InitialPosition.x() << " "  // Changed PostPosition -> Initial Position
			    << setprecision(6) << InitialPosition.y() << " " 
			    << setprecision(6) << InitialPosition.z() << " " 
				  /*
 			    << InitialEnergy << " " 
   			    << setprecision(6) << VertexPosition.x() << " "  // Changed PostPosition -> Initial Position
			    << setprecision(6) << VertexPosition.y() << " " 
			    << setprecision(6) << VertexPosition.z() << " " 
				  */         
                            << G4endl;
	    }

	}

     }
     
   }
  }

  //For Calculation of Backscattering electron in 2013.04.17
  if( ( fTrack->GetKineticEnergy() == 0 || fStepStatus == fGeomBoundary ) 
        && fTrack->GetVolume()->GetName() == "VirtualChamberCylinder"
        && ( ParticleID == -1 || ParticleID == 1 ) 
      ){
          if( ParticleTrackID == 1 || 
              ( ParticleTrackID == 2 && InitialVolume == 1 ) ){
   	        CLHEP::Hep3Vector VertexPosition = fStep->GetPostStepPoint()->GetPosition();
 	        OutputFile100 << ParticleID           << " " 
                              << ParticleTrackID      << " " 
 			                  << fTrack->GetKineticEnergy()/MeV  << " " 
                              << VertexPosition.x()   << " " 
                              << VertexPosition.y()   << " " 
                              << VertexPosition.z() 
                              << G4endl; 
                fTrack->SetTrackStatus(fStopAndKill); 
          }
          fTrack->SetTrackStatus(fStopAndKill);
  }

  //---------------------------------
  // Dmitri-san's FaradayCup
  // added in 2013.04.23
  //---------------------------------
  G4StepPoint* point1 = fStep->GetPreStepPoint();
  G4StepPoint* point2 = fStep->GetPostStepPoint();
  
  CLHEP::Hep3Vector VertexPosition = fStep->GetPostStepPoint()->GetPosition();

  G4VPhysicalVolume* volume1 = point1->GetPhysicalVolume();
  G4VPhysicalVolume* volume2 = point2->GetPhysicalVolume();

  G4String name1 = volume1->GetName();
  G4String name2 = "nothingness";
  if (volume2) name2 = volume2->GetName();

  G4bool entering_faraday_cup;
  G4bool leaving_faraday_cup;

  //G4bool entering_faraday_cup = (name1 == "World" && name2 == "FaradayCupDmitri");
  //G4bool leaving_faraday_cup  = (name1 == "FaradayCupDmitri" && name2 == "World");

  if( FaradayCup_version == 1 ){
      entering_faraday_cup = (name1 != "FCBody" && name2 == "FCBody");
      leaving_faraday_cup  = (name1 == "FCBody" && name2 != "FCBody");

   } else if ( FaradayCup_version == 4 ){
      entering_faraday_cup = (name1 != "FC4Body" && name2 == "FC4Body");
      leaving_faraday_cup  = (name1 == "FC4Body" && name2 != "FC4Body");

   }else if( FaradayCup_version == 5 ){ 
      entering_faraday_cup = (name1 != "FC5Body" && name2 == "FC5Body");
      leaving_faraday_cup  = (name1 == "FC5Body" && name2 != "FC5Body");

      // }else if(  == 1 ){
      // entering_faraday_cup = (name1 != "VirtualChamberTarget" && name2 == "VirtualChamberTarget");
      //leaving_faraday_cup  = (name1 == "VirtualChamberTarget" && name2 != "VirtualChamberTarget");

   }else{
      entering_faraday_cup = (name1 == "World" && name2 == "FaradayCupDmitri");
      leaving_faraday_cup  = (name1 == "FaradayCupDmitri" && name2 == "World"); 
  }

  G4ParticleDefinition* prt = fStep->GetTrack()->GetDefinition();

    if ( entering_faraday_cup || leaving_faraday_cup ){

       if ( entering_faraday_cup && prt->GetPDGCharge()!=0 ){
	 OutputFile100 << 1 << " " 
                       << prt->GetPDGCharge() << " " 
                       << ParticleTrackID << " " 
                       << VertexPosition.x() << " " 
                       << VertexPosition.y() << " " 
                       << VertexPosition.z() << " " 
                       << G4endl;
    }
       if ( leaving_faraday_cup && prt->GetPDGCharge()!=0 ){
	 OutputFile100 << 2 << " " 
		       << prt->GetPDGCharge() << " "
                       << ParticleTrackID << " "
                       << VertexPosition.x() << " " 
                       << VertexPosition.y() << " " 
                       << VertexPosition.z() << " " 
                       << G4endl;
       }
    }

  //---------------------------------
  // Energy Loss Output
  //---------------------------------
  if( OutputFile4_flag ){ 

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
	  //if( fTrack->GetVolume()->GetName()=="ICE" && fTrack->GetNextVolume() != 0 ){

	  /*
	    if( ProcessID  == 0 ){
	    //if( ProcessID  == 1 && fStep->GetTotalEnergyDeposit() +  fTrack->GetKineticEnergy() < 1.0 ){
            double PreX=fStep->GetPreStepPoint()->GetPosition().x();
            double PreY=fStep->GetPreStepPoint()->GetPosition().y();
            double PreZ=fStep->GetPreStepPoint()->GetPosition().z();
	    double PosX=fStep->GetPostStepPoint()->GetPosition().x();
	    double PosY=fStep->GetPostStepPoint()->GetPosition().y();
	    double PosZ=fStep->GetPostStepPoint()->GetPosition().z();

	    double PathLength=sqrt(   (PreX-PosX)*(PreX-PosX) 
				    + (PreY-PosY)*(PreY-PosY)
				    + (PreZ-PosZ)*(PreZ-PosZ) );

           OutputFile100
             //<< fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName() << " " 
             //<< ProcessID        << " "   
	     //<< ParticleID       << " "
	     //<< fTrack->GetTrackID() << " "
	     //<< ParticleTrackID  << " "

             <<  setprecision(10) <<  PathLength << "  "
             //<<  setprecision(10) <<  fStep->GetStepLength() << " "               

 	     <<  setprecision(10) <<  fStep->GetTotalEnergyDeposit()/eV << " "
             <<  setprecision(10) <<  fTrack->GetKineticEnergy()/eV     << " "

             //<< CreatedProcess << " " 
             //<< StepCount << " " 
             //<< fTrack->GetMomentum().x()/fTrack->GetMomentum().mag() << " " 
             //<< fTrack->GetMomentum().y()/fTrack->GetMomentum().mag() << " " 
             //<< fTrack->GetMomentum().z()/fTrack->GetMomentum().mag() << " " 
	     << G4endl;
	  }
          */

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
	      OutputFile4 
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
	      
	    }
	    
	    /* Tested in 2012.04.23 No.01 */
	    /*
	      OutputFile100   << ProcessID << " " 
	      << fStep->GetTotalEnergyDeposit()/MeV << " " 
	      << G4endl;    
	    */
	    
	  }
	  else{
	    fprintf(stderr, "energy_deposit_mesh_generator has some problem...\n");
	  }
          created_energy_deposit_meth.clear();

	  }
	  else{   
	    // The Energy deposit can not be used for yield of fluorescence photon...
           //  ... This volume is not World.
	    OutputFile4 
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

            /* Tested in 2012.04.23 No.01 */
	    /*
                OutputFile100   << "-1113" << " " 
                                << fStep->GetTotalEnergyDeposit()/MeV << " " 
                               << G4endl;    
	    */ 
	  }
	  
      }else if( fStep->GetTotalEnergyDeposit() > 0. &&
	        fStep->GetPostStepPoint()->GetPosition().z() < 0. ){
	

 	  // The Energy deposit can not be used for yield of fluorescence photon...                             
          // The height is under-ground.  
	OutputFile4 
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

          /* Tested in 2012.04.23 No.01 */
          /*
             OutputFile100   << "-1112" << " " 
                             << fStep->GetTotalEnergyDeposit()/MeV << " " 
                             << G4endl;    
	  */
      }
    
      if( fTrack->GetVolume()->GetName()=="World" && fTrack->GetNextVolume() == 0 ){
       OutputFile4 
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

          /* Tested in 2012.04.23 No.01 */
          /*
              OutputFile100   << "-1111" << " " 
                              << fTrack->GetKineticEnergy()/MeV << " " 
                              << G4endl;    
	  */                   
     }

  }

  // =========================================================================
  // Top Plate
  // =========================================================================
  /*
  if( fTrack->GetVolume()->GetName()=="TopPlate" ){
          OutputFile100
	    << " " << setprecision(6) << fTrack->GetPosition().x()
	    << " " << setprecision(6) << fTrack->GetPosition().y()
	    << " " << setprecision(6) << fTrack->GetPosition().z()
	    << G4endl;
  }
  */

  // =========================================================================
  // Screen Monitor 3 Hit Position
  // =========================================================================
  if( OutputFile5_flag == 1 ){
    //if( abs(ParticleID) == 1 ){
      if( fTrack->GetVolume()->GetName()=="SM3" ){
          OutputFile5  
	       << " " << setprecision(5)  << ParticleID 
           << " " << setprecision(5)  << ParticleTrackID 
	       << " " << setprecision(5)  << fTrack->GetTrackID()
  	       << " " << setprecision(6)  << fTrack->GetPosition().x()/cm
	       << " " << setprecision(6)  << fTrack->GetPosition().y()/cm
 	       << " " << setprecision(6)  << fTrack->GetPosition().z()/cm
	       //<< " " << setprecision(10) << fTrack->GetKineticEnergy()/MeV 
	       //<< " " << setprecision(6)  << fTrack->GetGlobalTime()
               //<< " " << setprecision(10) << fStep->GetTotalEnergyDeposit()/eV 
               << G4endl;

      }
      //}
  }

  // =========================================================================
  // Beam Attenuator added in 2012.06.25
  // =========================================================================
  /*
  if( ParticleTrackID == 1 ){
      if( fTrack->GetVolume()->GetName()=="BeamAttenuator" ){
          OutputFile100
	    << " " << setprecision(6)  << fTrack->GetPosition().x()/cm
	    << " " << setprecision(6)  << fTrack->GetPosition().y()/cm
	    << " " << setprecision(6)  << fTrack->GetPosition().z()/cm
	    << " " << setprecision(10) << fTrack->GetKineticEnergy()/MeV
	    << G4endl;
      }
  }
  */

  // =========================================================================
  // Virtual Chamber
  // =========================================================================
  /*
  if( fTrack->GetVolume()->GetName()=="VirtualChamberCylinder" ||
      fTrack->GetVolume()->GetName()=="VirtualChamberTarget"  ){
        OutputFile100	  
	     << " " << setprecision(6) << fTrack->GetPosition().x()/cm
	     << " " << setprecision(6) << fTrack->GetPosition().y()/cm
	     << " " << setprecision(6) << fTrack->GetPosition().z()/cm
	     << G4endl;
  }
  */
  /*
     if( fTrack->GetVolume()->GetName()=="VirtualChamberBody" || 
         fTrack->GetVolume()->GetName()=="VirtualChamberPipe" ){
       int VID=-1;
       if( fTrack->GetVolume()->GetName()=="VirtualChamberBody" ) VID=0;
       if( fTrack->GetVolume()->GetName()=="VirtualChamberPipe" ) VID=1;
         OutputFile100   
     	       << " " << VID 
    	       << " " << setprecision(10) << fStep->GetTotalEnergyDeposit()/eV 
               << " " << setprecision(6)  << fTrack->GetPosition().x()
	       << " " << setprecision(6)  << fTrack->GetPosition().y()
 	       << " " << setprecision(6)  << fTrack->GetPosition().z()
               << G4endl;
     }
  */

     //==================================================
     // tested in 2012.10.10
     /*
     if( ParticleTrackID == 2 && 
         abs(ParticleID) == 0 &&
	 fStep->GetPostStepPoint()->GetPosition().z()/m > 2.4 &&
	 ( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="phot" ||
           fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="compt" )
         ){
       OutputFile100
  	  << fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()
	  << G4endl;
     }
     */

  //G4cout.precision(prec);

  // Check The Cerenkov Scattering Process, 2013.09.09 
  if ( ParticleID == 50 ){
	if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="OpRayleigh" ){
          RayScatNum++;
          OpticalPhotonSPoint=fStep->GetPostStepPoint()->GetPosition();
	} 
	if( fStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()=="OpMieHG" ){
          MieScatNum++;
          OpticalPhotonSPoint=fStep->GetPostStepPoint()->GetPosition();
	}
   }

  if( OutputFile6_flag ){
   if( ParticleID == 50 ){ 

     if( fTrack->GetVolume()->GetName()=="Mirror6" || fTrack->GetVolume()->GetName()=="Mirror7" ){
       
       OpticalPhotonHPoint=fStep->GetPostStepPoint()->GetPosition();
       if( RayScatNum > 0 || MieScatNum > 0 ){
	 OpticalPhotonVStart=OpticalPhotonSPoint;
       }else{
	 OpticalPhotonVStart=OpticalPhotonGPoint;
       }
       OpticalPhotonVStartEarth=de_mesh->PositionTransportatinoELStoEARTCH(OpticalPhotonVStart/m);
       OpticalPhotonVEndEarth=de_mesh->PositionTransportatinoELStoEARTCH(OpticalPhotonHPoint/m);
       OpticalPhotonVector=(OpticalPhotonVEndEarth-OpticalPhotonVStartEarth).unit();
       
       if( fTrack->GetVolume()->GetName()=="Mirror6" ) HitCamera = 6; 
       if( fTrack->GetVolume()->GetName()=="Mirror7" ) HitCamera = 7; 
       
       OutputFile6
	 << " " << setprecision(10) << PhotonEnergytoLambda(fTrack->GetKineticEnergy()/eV)
	 
	 << " " << HitCamera
	 
	 << " " << setprecision(3)  << RayScatNum  
	 << " " << setprecision(3)  << MieScatNum
	 /*
	   << " " << setprecision(8)  << OpticalPhotonHPoint.x()/cm   // = mm -> cm
	   << " " << setprecision(8)  << OpticalPhotonHPoint.y()/cm
	   << " " << setprecision(8)  << OpticalPhotonHPoint.z()/cm
	 */
	 << " " << setprecision(8)  << OpticalPhotonVStartEarth.x()*100.0 // m -> cm
	 << " " << setprecision(8)  << OpticalPhotonVStartEarth.y()*100.0
	 << " " << setprecision(8)  << OpticalPhotonVStartEarth.z()*100.0
	 
	 << " " << setprecision(8)  << OpticalPhotonVector.x() 
	 << " " << setprecision(8)  << OpticalPhotonVector.y() 
	 << " " << setprecision(8)  << OpticalPhotonVector.z() 
	 
	 << endl;                       
       
       fTrack->SetTrackStatus(fStopAndKill);       
     }
     
   }
  }
  

  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
  // This will be the output I want to analyze.  I am 
  // not really following what the other output is 
  // dumping...
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

  
  // New output code
  // First let's try to construct the shower shape
  if( fTrack->GetVolume()->GetName() == "ICE" ){
  // Only consider e+ and e-
  //if( ParticleID == 1 || ParticleID == -1 ){
    MattOut <<
      fStep->GetPostStepPoint()->GetPosition().x()/m <<" "<<
      fStep->GetPostStepPoint()->GetPosition().y()/m <<" "<<
      fStep->GetPostStepPoint()->GetPosition().z()/m <<" "<<
      fTrack->GetKineticEnergy()/MeV<<" "<<
      fStep->GetTotalEnergyDeposit()/MeV<<" "<<
      ParticleID<<" "<<
      fTrack->GetTrackID()<<" "<<
      G4endl;	       
  }// end check

      
    
    

  StepCount++;
}
void SteppingVerbose::TrackingStarted_OLD()
{
  CopyState();
}

void SteppingVerbose::TrackingStarted()
{

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
            RayScatNum=0;
            MieScatNum=0;
            } 

   else if( ParticleName == "geantino" ){ ParticleID = -100; }

   else{ ParticleID = 100; }

   // Track ID : Primary particle is 1, Secondary Particle is 2;
   if( fTrack->GetTrackID() == 1 )
  { 
       ParticleTrackID = 1; 
   }else{
       ParticleTrackID = 2;
   }

   // TrackPID[ fTrack->GetTrackID() ] = ParticleID; 

   for( int i=0; i<200; i++ ) ParticleTrackID_FLAG[i] = 0;
   
   // Initial Energy & Mometum & Position  
   InitialEnergy   = fTrack->GetKineticEnergy();
   InitialMomentum = fTrack->GetMomentum(); //fTrack->GetMomentum().mag();
   InitialPosition = fTrack->GetPosition();
   if( ParticleTrackID == 50 ) OpticalPhotonGPoint = InitialPosition;

   if( ParticleID== -1 && ParticleTrackID == 1 ) PrimarElevtronEnergy = InitialEnergy/eV;;

   InitialVolumeID = -1;
   StoppedVolumeID = -1;

   if( ParticleTrackID == 2 &&
       abs(ParticleID) == 1 &&
       InitialPosition.z()/m > 2.4 ){

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

      if( fTrack->GetVolume()->GetName() == "FC4Body" ) InitialVolumeID = 1;
      if( fTrack->GetVolume()->GetName() == "FC4GNDLayer" ) InitialVolumeID = 2;
      if( fTrack->GetVolume()->GetName() == "FC4ShieldLayer" ) InitialVolumeID = 3;
      if( fTrack->GetVolume()->GetName() == "FC4TopPlate" ) InitialVolumeID = 4;
      if( fTrack->GetVolume()->GetName() == "FC4Iso1" ) InitialVolumeID = 5;
      if( fTrack->GetVolume()->GetName() == "FC4Iso2" ) InitialVolumeID = 6;

      if( fTrack->GetCreatorProcess()->GetProcessName()=="msc" )   CreatedProcess = 0;
      if( fTrack->GetCreatorProcess()->GetProcessName()=="eIoni" ) CreatedProcess = 1;
      if( fTrack->GetCreatorProcess()->GetProcessName()=="eBrem" ) CreatedProcess = 2;
      if( fTrack->GetCreatorProcess()->GetProcessName()=="phot" )  CreatedProcess = 3;
      if( fTrack->GetCreatorProcess()->GetProcessName()=="compt" ) CreatedProcess = 4;
      if( fTrack->GetCreatorProcess()->GetProcessName()=="conv" )  CreatedProcess = 5;
      
   }

   StepCount=0;

  //-------------------
  // Generated Data
  //-------------------
  if( ParticleTrackID == 2 ){ 
   if( OutputFile1_flag ){
     /*
     if( ( abs(ParticleID) == 1 || abs(ParticleID) == 0 )  &&
           fTrack->GetParentID() == 1 && 
	   fTrack->GetVolume()->GetName() == DP_NAME[0] ){  
     */
     // if( abs(ParticleID) == 1 || abs(ParticleID) == 0 ){
     if( abs(ParticleID) == 1 ){
       //if( InitialEnergy/keV < 1.0 ){
       if( fTrack->GetCreatorProcess()->GetProcessName() == "phot" || 
           fTrack->GetCreatorProcess()->GetProcessName() == "compt" ){

            OutputFile1 << " " 
	      //<< setprecision(5) << fTrack->GetTrackID()       << " "
	      //<< setprecision(5) << fTrack->GetParentID()      << " "
              //<< setprecision(5)  << ParticleTrackID           << " "
	          << setprecision(5)  << ParticleID                << " "
              << setprecision(15) << InitialEnergy/keV         << " "
	      // << setprecision(10) << InitialPosition.x()/m     << " "
              // << setprecision(10) << InitialPosition.y()/m     << " "      
              // << setprecision(10) << InitialPosition.z()/m     << " "
	         << fTrack->GetCreatorProcess()->GetProcessName() << " " 
	          << G4endl;
       }

     }
   }
  }

  //---------------------------------
  // For Faraday Cup Charge Measurement
  //---------------------------------
  InitialVolume = 0;
  StoppedVolume = 0;
  ChargeDeposit_in_FaradayCup=0;

  if( OutputFile3_flag ){    
    if( abs(ParticleID) == 1 ){  // Paritcle is only electron and positrons           
      // Start Volume is Faraday Cup Body or not
      InitialVolumeName = fTrack->GetVolume()->GetName();
      if( InitialVolumeName == "FCBody"  || 
	  InitialVolumeName == "FC4Body" ||
          InitialVolumeName == "FC5Body" || 
          InitialVolumeName == "FaradayCupDmitri" ||  
          InitialVolumeName == "VirtualChamberTarget" ) InitialVolume = 1;
      if( InitialVolumeName == "FC4ShieldLayer" ) InitialVolume = 2;

      // Dmitri-san's geometry study
      if( InitialVolumeName == "FCDmitriUTi" ) InitialVolume = 2;
      if( InitialVolumeName == "FCDmitriUCu" ) InitialVolume = 3;
      if( InitialVolumeName == "FCDmitriBCu" ) InitialVolume = 4;
      if( InitialVolumeName == "FCDmitriBTi" ) InitialVolume = 5;
      
    }
  }

  //--------------------------------------------------------
  // Kill paritcles except gamma,e+-,opticalphoton,geantino 
  //--------------------------------------------------------
  if( ParticleID == 2 || ParticleID == -2 ||
      ParticleID == 3 || ParticleID == -3 ||
      ParticleID == 4 || ParticleID == -4 ||
      ParticleID == 5 || ParticleID == -5 ||
      ParticleID == 100 ) {
    fTrack->SetTrackStatus(fStopAndKill);
  }


  //G4cout.precision(prec);
}


