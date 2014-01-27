//===================================================================
//      SteppingVerbose.cc 
//      Author    : T.Shibata
//      Created     : 2010.02.11 
//      Last Update : 2011.06.11 
//===================================================================  
#ifndef SteppingVerbose_h
#define SteppingVerbose_h 1
                                                                                               
#include <iostream>
#include <fstream>

#include "G4SteppingVerbose.hh"

#include "ELSParameters.hh"
#include "FDParameters.hh"
#include "PhysicsParameters.hh"

#include "ReadFile.hh"

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Vector/Rotation.h>
#include "EnergyDepositMesh.hh"

class Plist;
class ELSParameters;
class EnergyDepositMeth;
class dE_Mesh;
class FDParameters;

//-----------------------------------------------
using namespace std;
//------------------------------------------------

class SteppingVerbose;
class SteppingVerbose : public G4SteppingVerbose
{
  public:
  SteppingVerbose(  Plist *plist,
                    char* Fname0, 
                    char* Fname1,             
                    char* Fname2,
                    char* Fname3,
                    char* Fname4,
		    char* mattOut);

   ~SteppingVerbose();

   void StepInfo();
   void StepInfo_OLD();
   void TrackingStarted();
   void TrackingStarted_OLD();

   ELSParameters *elsparam;

   dE_Mesh *de_mesh;
   typedef vector<EnergyDepositMeth> CreatedEnergyDepositMeth;

   FDParameters *fdparam;

  private:

    CLHEP::Hep3Vector position_of_beaminjection;
    CLHEP::Hep3Vector position_of_mirror_6_BRM_coordinate;
    double mesh_dv[3];
    double mesh_vregion[3];
  
    ofstream OutputFile0;
    ofstream OutputFile1;
    ofstream OutputFile2;
    ofstream OutputFile3;
    ofstream OutputFile4;
    ofstream OutputFile5;
    ofstream OutputFile6;
    ofstream MattOut;

    ofstream OutputFile100;   // Debug outfile

    G4int    OutputFile0_flag; 
    G4int    OutputFile1_flag; 
    G4int    OutputFile2_flag; 
    G4int    OutputFile3_flag; 
    G4int    OutputFile4_flag;
    G4int    OutputFile5_flag; 
    G4int    OutputFile6_flag; 
  
    // Center of Beam injector
    G4double CenterofBeamInjection[3];
    G4double DirectionofBeamInjection[3];

    // Output parameters 
    G4double PrimarElevtronEnergy;

    G4int    ParticleID;
    G4int    ParticleParentID;
    G4int    ParticleTrackID;
    G4int    ParticleTrackID_FLAG[200];
    G4int    TrackPID[10000];
    G4int    ParentTrackID[10000];

    G4int    CreatedProcess;    
    G4int    StepCount;

    int      FaradayCup_version;
   
    G4int    InitialVolume;
    G4int    InitialVolumeID;

    G4double InitialEnergy;
    G4ThreeVector InitialPosition;
    G4ThreeVector InitialMomentum;
    G4int          StoppedVolume;
    G4int          StoppedVolumeID;
    G4ThreeVector  StoppedPosition;
   
    G4String       InitialVolumeName;
    G4String       StoppedVolumeName;

    G4int          ChargeDeposit_in_FaradayCup;
    G4int          ChargeDeposit_in_FaradayCup_Shield;
     
    G4int          DetectorPlaneID;
    G4double       DetectorPlaneEnergy;
    G4ThreeVector  DetectorPlanePosition;
    G4ThreeVector  DetectorPlaneMomentum;
    G4double       DetectorPlaneScatteringDirection[3]; 
    G4double       DetectorPlaneScatteringMag;
    G4double       DetectorPlaneScatteringAngle;
 
    G4int          TiWindowSLitID;
    G4double       TiWindowEnergy;
    G4ThreeVector  TiWindowMomentum;
    G4ThreeVector  TiWindowPosition;
    
    G4String       GDRGeneratedParticle;
    G4double       GDRGammaEnergy;
    G4ThreeVector  GDRGPosition;
    G4double       GDRGeneratedParticleEnergy;
    G4ThreeVector  GDRGeneratedParticleMomemtum;   

    // Cerenkov Photon Study, 2013.09.12
    G4int RayScatNum;
    G4int MieScatNum;
    G4int HitCamera;

    G4ThreeVector OpticalPhotonGPoint;
    G4ThreeVector OpticalPhotonSPoint;
    G4ThreeVector OpticalPhotonHPoint;

    G4ThreeVector OpticalPhotonVStart;
    G4ThreeVector OpticalPhotonVEnd;

    G4ThreeVector OpticalPhotonVStartEarth;
    G4ThreeVector OpticalPhotonVEndEarth;
    G4ThreeVector OpticalPhotonVector;

  
};
#endif
