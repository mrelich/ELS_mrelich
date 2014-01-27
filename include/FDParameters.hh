//===================================================================
//    FD Parameter Class : FDParameters.hh
//      Author    : T.Shibata
//
//      Creation  : 2006.03.07
//
//    Last Update : 2006.04.15
//    Modified    : 2006.03.08 added PhysicsParameters.hh
//    Modified    : 2006.03.31 added FD Camera Box
//                                   FD Camera Frame
//                                   PMT & BG4 & Paraglass
//                                   FD Station
//    Modified    : 2006.04.04 added FD Station Wall & Roof Geometory                 
//    Modified    : 2006.04.15 added Segment Mirror 
//
//===================================================================
#include <iostream>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <string>

#include <CLHEP/Vector/ThreeVector.h>
#include "G4ThreeVector.hh"

#ifndef FDParameters_h
#define FDParameters_h
//-----------------------------------------------
using namespace std;
//-----------------------------------------------
class FDParameters{
public:
      //Contstruction
      FDParameters();
     //Deconstruction
     ~FDParameters();

      CLHEP::Hep3Vector GetFDStation(void);
      CLHEP::Hep3Vector GetFDStationWall( int i );
      CLHEP::Hep3Vector GetFDStationRoof( int i );

      CLHEP::Hep3Vector GetMirrorCenter(int CID);
      CLHEP::Hep3Vector GetMirror(int CID);
 
      CLHEP::Hep3Vector GetSegMirrorPos( int SegID );
      CLHEP::Hep3Vector GetSegMirrorRot( int SegID );

      CLHEP::Hep3Vector GetFrameCenter( int FID );
      CLHEP::Hep3Vector GetFrameMain( int SID, int RL );

      CLHEP::Hep3Vector GetCameraBoxCenter( int CID );
      CLHEP::Hep3Vector GetCameraBox( int CID, int PID );

      CLHEP::Hep3Vector GetParaglassCenter( int CID );

      CLHEP::Hep3Vector GetCenterPMTSurcafe( int CID );
      CLHEP::Hep3Vector GetCenterPMTs( int CID );
      CLHEP::Hep3Vector GetPMT( int CID, int PMTID );

      CLHEP::Hep3Vector GetCenterPMTSs( int CID );
      CLHEP::Hep3Vector GetPMTSen( int CID, int PMTID );

      CLHEP::Hep3Vector GetCenterBG3s( int CID );
      CLHEP::Hep3Vector GetBG3( int CID, int PMTID );

     //FD Station(Building) Geometory 
      double fdstationbackwallwidth(void){ return FDStationBackWallWidth; }
      double fdstationfrontwallwidth(void){ return FDStationFrontWallWidth; }
      double fdstationwallheight(void){ return FDStationWallHeight; }

      double fdshutterbiasheight(void){ return FDShutterBiasHeight; }
      double fdshutterwidth(void){ return FDShutterWidth; }
      double fdshutterheight(void){ return FDShutterHeight; }

      double fdstationbaseangle(void){ return FDStationBaseAngle; }

      double fdcamerar(void){ return FDCameraR; }
      double fdsamerahalfdistance(void){ return FDCameraHalfDistance; }
      double fdcamerar2(void){ return FDCameraR2; }

      double rotationangle( int i ){ return RotationAngle[i]; } 

      double fdstationwallfrontgeometory( int i ){ return FDStationWallfrontGeometory[i]; }
      double fdstationwallbackgeometory( int i ){ return FDStationWallbackGeometory[i]; }
      double fdstationshuttergeometory( int i ){ return FDStationShutterGeometory[i]; }
      double fdstationroofgeometory( int i ){ return FDStationRoofGeometory[i]; }
 
      double wallxrotationangle( int i ){ return WallXRotationAngle[i]; }
      double wallzrotationangle( int i ){ return WallZRotationAngle[i]; }
      double roofxrotationangle( int i ){ return RoofXRotationAngle[i]; }
      double roofzrotationangle( int i ){ return RoofZRotationAngle[i]; }
  
      //Mirror Geometory
      double innerradiusmirror(void){ return InnerRadiusMirror; }
      double outerradiusmirror(void){ return OuterRadiusMirror; }
      double diametermirror(void){ return DiameterMirror; }
      double depthmirror(void){ return DepthMirror; }
 
      double distancefromcameraframe( int i ){ return DistanceFromCameraFrame[i]; }
      double slopeanglemirror( int i ){ return SlopeAngleMirror[i]; }
      double heigthmirror( int i ){ return HeigthMirror[i]; }

      double mirrorphi(void){ return MirrorPhi; }
      double deltamirrorphi(void){ return DeltaMirrorPhi; }

      double posimirrorcx( int i ){ return PosiMirrorCX[i]; }
      double posimirrorcy( int i ){ return PosiMirrorCY[i]; }
      double posimirrorcz( int i ){ return PosiMirrorCZ[i]; }

      double sizesegmirror(void){ return SizeSegMirror; }

      double segmirrorposx( int i ){ return SegMirrorPosX[i]; }
      double segmirrorposy( int i ){ return SegMirrorPosY[i]; }
      double segmirrorposz( int i ){ return SegMirrorPosZ[i]; }
   
      double segmirrorrotx( int i ){ return SegMirrorRotX[i]; }
      double segmirrorroty( int i ){ return SegMirrorRotY[i]; }
      double segmirrorrotz( int i ){ return SegMirrorRotZ[i]; }

      // Camera Box
      double cameraboxlx( int i ){ return CameraBoxLX[i]; }
      double cameraboxly( int i ){ return CameraBoxLY[i]; }
      double cameraboxlz( int i ){ return CameraBoxLZ[i]; }

      double framemainwidth(void){ return FrameMainWidth; }
      double framemainheight(void){ return FrameMainHeight; }
      double framemaindepth(void){ return FrameMainDepth; }

      // Paraglass 
      double paraglasswidth(void){ return ParaglassWidth; }
      double paraglassheight(void){ return ParaglassHeight; }
      double paraglassthick(void){ return ParaglassThick; }

      // PMT
      int nopmt(void){ return NoPMT; }

      int nopmtside(void){ return NoPMTSide; }
      int nopmtzplane(void){ return NoPMTZPlanes; }
      double irpmt(int i){ return irPMT[i]; }
      double orpmt(int i){ return orPMT[i]; }
      double zpmt(int i){ return zPMT[i]; }

      int nopmtsside(void){ return NoPMTSide; }
      int nopmtszplane(void){ return NoPMTZsPlanes; }
      double irpmts(int i){ return irPMTs[i]; }
      double orpmts(int i){ return orPMTs[i]; }
      double zpmts(int i){ return zPMTs[i]; }

      int nobg3side(void){ return NoBG3Side; }
      int nobg3zplane(void){ return NoBG3ZPlanes; }
      double irbg3(int i){ return irBG3[i]; }
      double orbg3(int i){ return orBG3[i]; }
      double zbg3(int i){ return zBG3[i]; }
  
      // DUMP
      void DumpPosition(void);
      bool Cameraframe( CLHEP::Hep3Vector InputPosition, int MirrorID );
      int  PMTfame( CLHEP::Hep3Vector InputPosition, int MirrorID );

      bool MirrorHit( int MID, CLHEP::Hep3Vector HitMirrorPoint );

      double RelativeHitDistance( int CID, CLHEP::Hep3Vector HitPoint );
      CLHEP::Hep3Vector RelativeHitPoint( int CID, int PID, CLHEP::Hep3Vector HitPoint );   

      bool FieldOfViewofCamera(int CID, CLHEP::Hep3Vector CreationPoint,
	  		       CLHEP::Hep3Vector ShiftMirror );

private:

     //FD Station(Building) Geometory 
     double FDStationBackWallWidth;
     double FDStationFrontWallWidth;
     double FDStationWallHeight;
     double FDStationWallDHeight;
 
     double FDStationWallDepth;
     double FDStationRoofDepth;

     double FDShutterBiasHeight;
     double FDShutterWidth;
     double FDShutterHeight;

     double FDStationBaseAngle;
  
     double FDCameraR;
     double FDCameraHalfDistance;
     double FDCameraR2;

     double FDStationWallfrontGeometory[3];
     double FDStationShutterGeometory[3];
     double FDStationWallbackGeometory[3];
     double FDStationRoofGeometory[5];
        
     CLHEP::Hep3Vector positionFDStation; 
     CLHEP::Hep3Vector positionFDStationWall[5];
     CLHEP::Hep3Vector positionFDStationRoof[3];

     double WallXRotationAngle[5];
     double WallZRotationAngle[5];
 
     double RoofXRotationAngle[3];
     double RoofZRotationAngle[3];

     //Rotation Angle List
     double RotationAngle[12]; 

     //Mirror Geometory
     double InnerRadiusMirror;
     double OuterRadiusMirror;
     double DiameterMirror;
     double DepthMirror;

     double MirrorPhi;
     double DeltaMirrorPhi;

     double DistanceFromCameraFrame[12];
     double SlopeAngleMirror[12];
     double HeigthMirror[12]; 

     CLHEP::Hep3Vector positionMirrorCenter[12];
     CLHEP::Hep3Vector positionMirror[12];

     double PosiMirrorCX[12];
     double PosiMirrorCY[12];
     double PosiMirrorCZ[12];

     double SizeSegMirror;

     double SegMirrorPosX[18];
     double SegMirrorPosY[18];
     double SegMirrorPosZ[18];

     double SegMirrorRotX[18];
     double SegMirrorRotY[18];
     double SegMirrorRotZ[18];

     //Paraglass Geometory  
     double ParaglassWidth;
     double ParaglassHeight;
     double ParaglassThick;
     double DistanceFromCBSurface;

     CLHEP::Hep3Vector positionParaglassCenter[12];

     //Camera Box Geometory
     double CameraBoxWidth;
     double CameraBoxHeight;
     double CameraBoxThick;
     double CameraBoxThickW;   
     double CameraBoxDepth;
     double CameraDistance;
     
     double CameraBoxLX[5];
     double CameraBoxLY[5];
     double CameraBoxLZ[5];
     
     CLHEP::Hep3Vector positionCameraBoxCenter[12];
     CLHEP::Hep3Vector positionCameraBox[12][5];
 
     //Camera Frame Geometory
     double FrameMainWidth;
     double FrameMainHeight;
     double FrameMainDepth;

     double FrameTubeDiameter;
     double FrameTubeLength;

     double FrameCrossBarWidth;
     double FrameCrossBarDepth;
     double FrameCrossBarThick;
  
     CLHEP::Hep3Vector positionFrameCenter[6];

     CLHEP::Hep3Vector positionFrameMain[6][2];
     CLHEP::Hep3Vector positionFrameTube[6][10];
     CLHEP::Hep3Vector positionFrameCrossBar[6][8];

     //PMT Geometory
     int NoPMT;
     int NoPMTSide;
     int NoPMTZPlanes;

     double irPMT[6];
     double orPMT[6];
     double zPMT[6];

     int NoPMTZsPlanes;

     double irPMTs[3];
     double orPMTs[3];
     double zPMTs[3];

     double PMTShiftH;
     double PMTShiftV;
   
     double PMTSpacingH;
     double PMTSpacingV;

     double DistanceFromCBSurfacePMT;
     double DepthSensitivePMT;

     CLHEP::Hep3Vector positionCenterPMTSurface[12];

     CLHEP::Hep3Vector positionCenterPMTs[12];
     CLHEP::Hep3Vector positionPMT[12][256];

     CLHEP::Hep3Vector positionCenterPMTSs[12];
     CLHEP::Hep3Vector positionPMTSen[12][256];

     //BG3 Geometory
     int NoBG3;
     int NoBG3Side;
     int NoBG3ZPlanes;
     double BG3Depth;
     double Space_BG3toPMT;

     double irBG3[3];
     double orBG3[3];
     double zBG3[3];

     double BG3ShiftH;
     double BG3ShiftV;
   
     double BG3SpacingH;
     double BG3SpacingV;

     CLHEP::Hep3Vector positionCenterBG3s[12];
     CLHEP::Hep3Vector positionBG3[12][256];

     double CameraAngle;

};
#endif
