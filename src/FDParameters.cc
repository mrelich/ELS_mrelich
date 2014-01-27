//===================================================================
//    FD Parameter Class : FDParameters.cc
//      Author    : T.Shibata
//
//      Creation  : 2006.03.07
//
//    Last Update : 2006.01.08
//    Modified    : 2006.03.08 added PhysicsParameters.hh
//                             but Nout Used in DetectorConstructor.cc yet
//    Modified    : 2006.03.31 added FD Camera Box
//                                   FD Camera Frame
//                                   PMT & BG4 & Paraglass
//                                   FD Station
//                : 2006.04.01 bugfixed
//    Modified    : 2006.04.04 added FD Station Wall & Roof Geometory
//    Modified    : 2006.04.15 added Segment Mirror
//                : 2006.04.18 bugfixed
//    Modified    : 2006.05.24 bugfixed
//    Modified    : 2006.05.28 bugfixed 
//    Modified    : 2006.01.01 modified PMT position  
//    Modified    : 2007.01.08 changed Camera ID 
//
//===================================================================
#include "FDParameters.hh"
#include "PhysicsParameters.hh"
//-----------------------------------------------
using namespace std;
//------------------------------------------------
FDParameters::FDParameters()
{
  // -----------------------------------------------------
  //  Unit 
  //     Lenght = m;
  //     Angle  = radian;
  // -----------------------------------------------------
 
  // -----------------------------------------------------
  // FD Station(Building) Geometory
  // -----------------------------------------------------

  FDStationBackWallWidth   = 18.69363;
  FDStationFrontWallWidth  = 11.54178;
  FDStationWallHeight      = 11.12423;
                            
  FDShutterBiasHeight      = 1.2192;
  FDShutterWidth           = 8.001;
  FDShutterHeight          = 8.001; 

  FDStationWallDepth       = 0.2;
  FDStationRoofDepth       = 0.2;

  FDStationBaseAngle       = 36.0*RADI;
  FDCameraR                = 15.711;
  FDCameraHalfDistance     = 1.754;
  FDCameraR2               = sqrt( FDCameraR*FDCameraR  
				   - FDCameraHalfDistance*FDCameraHalfDistance );  

  positionFDStation.setX( 0.0 );
  //positionFDStation.setY( 100.0 );
  positionFDStation.setY( 0.0 );
  positionFDStation.setZ( 0.0 );

  FDStationWallDHeight = FDStationWallHeight - ( FDShutterBiasHeight + FDShutterHeight );
  FDStationWallDHeight = FDStationWallDHeight - FDShutterBiasHeight; // = 0.68483

  //--------------------------------
  FDStationWallfrontGeometory[0] = FDStationFrontWallWidth;
  FDStationWallfrontGeometory[1] = FDStationWallHeight + FDStationWallDHeight;
  FDStationWallfrontGeometory[2] = FDStationWallDepth;

  FDStationShutterGeometory[0] = FDShutterWidth;
  FDStationShutterGeometory[1] = FDShutterHeight;
  FDStationShutterGeometory[2] = FDStationWallDepth;

  FDStationWallbackGeometory[0]  = FDStationBackWallWidth;
  FDStationWallbackGeometory[1]  = FDStationWallHeight + FDStationWallDHeight; 
  FDStationWallbackGeometory[2]  = FDStationWallDepth;
  
  FDStationRoofGeometory[0] = FDStationBackWallWidth*sin(18.0*RADI);
  FDStationRoofGeometory[1] = 0.0;
  FDStationRoofGeometory[2] = FDStationRoofDepth;
  FDStationRoofGeometory[3] = FDStationRoofDepth;
  FDStationRoofGeometory[4] = 0.5*FDStationBackWallWidth*cos(18.0*RADI);
  //--------------------------------

  /* --------------------------------
      Numbering of FD Station Wall
      ----------------------------
            /\
           /  \
        4 /    \ 0
         /      \
        /        \
        \        /
       3 \______/ 1
            2    
  ------------------------------------
  */ 

  double WallHeight = 0.5*( FDStationWallHeight - FDStationWallDHeight );
 
  positionFDStationWall[0].setX( 0.5*FDStationBackWallWidth*sin(54.0*RADI) );
  positionFDStationWall[0].setY( positionFDStation.y() -
  				 0.5*FDStationBackWallWidth*cos(54.0*RADI) );
  positionFDStationWall[0].setZ( WallHeight );

  positionFDStationWall[1].setX( FDStationBackWallWidth*cos(18.0*RADI)*sin(36.0*RADI) );
  positionFDStationWall[1].setY( positionFDStation.y() -
                                 FDStationBackWallWidth*cos(18.0*RADI)*cos(36.0*RADI) );
  positionFDStationWall[1].setZ( WallHeight );

  positionFDStationWall[2].setX( 0.0 );
  positionFDStationWall[2].setY( positionFDStation.y() -
                                  FDStationBackWallWidth*cos(18.0*RADI) );
  positionFDStationWall[2].setZ( WallHeight );

  positionFDStationWall[3].setX( -FDStationBackWallWidth*cos(18.0*RADI)*
				  sin(36.0*RADI) );
  positionFDStationWall[3].setY( positionFDStation.y() -
                                 FDStationBackWallWidth*cos(18.0*RADI)*cos(36.0*RADI) );
  positionFDStationWall[3].setZ( WallHeight );
  
  positionFDStationWall[4].setX( -0.5*FDStationBackWallWidth*sin(54.0*RADI) );
  positionFDStationWall[4].setY( positionFDStation.y() -
                                 0.5*FDStationBackWallWidth*cos(54.0*RADI) );
  positionFDStationWall[4].setZ( WallHeight );


  double RoofHeight = 2.0*WallHeight + FDStationRoofDepth/2.0;
  positionFDStationRoof[0].setX( 0.5*FDStationBackWallWidth*cos(18.0*RADI)*sin(36.0*RADI) );   
  positionFDStationRoof[0].setY( positionFDStation.y() - 
				 0.5*FDStationBackWallWidth*cos(18.0*RADI)*cos(36.0*RADI) );
  positionFDStationRoof[0].setZ( RoofHeight );

  positionFDStationRoof[1].setX( 0.0 );
  positionFDStationRoof[1].setY( positionFDStation.y() -
                                 0.5*FDStationBackWallWidth*cos(18.0*RADI) );
  positionFDStationRoof[1].setZ( RoofHeight );
 
  positionFDStationRoof[2].setX( -0.5*FDStationBackWallWidth*cos(18.0*RADI)*sin(36.0*RADI) );
  positionFDStationRoof[2].setY( positionFDStation.y() -
                                 0.5*FDStationBackWallWidth*cos(18.0*RADI)*cos(36.0*RADI) );
  positionFDStationRoof[2].setZ( RoofHeight );

  //-----------------------------------------------------

  for( int i=0; i<5; i++ ) WallXRotationAngle[i]= 90.0*RADI;
  WallZRotationAngle[0] = WallZRotationAngle[3] = -36.0*RADI;
  WallZRotationAngle[1] = WallZRotationAngle[4] =  36.0*RADI;
  WallZRotationAngle[2] = 0.0*RADI;

  for( int i=0; i<5; i++ ) RoofXRotationAngle[i]= -90.0*RADI;
  RoofZRotationAngle[0] =  36.0*RADI;
  RoofZRotationAngle[1] = 0.0*RADI;
  RoofZRotationAngle[2] = -36.0*RADI;
 
  //------------------------------------------------------

  positionFrameCenter[1].setX( FDCameraR*sin(FDStationBaseAngle)
                               + FDCameraHalfDistance*cos(FDStationBaseAngle) );
  positionFrameCenter[1].setY( positionFDStation.y()
			       - FDCameraR*cos(FDStationBaseAngle) 
			       + FDCameraHalfDistance*sin(FDStationBaseAngle) );
  positionFrameCenter[1].setZ( 0.0 );


  positionFrameCenter[0].setX( FDCameraR*sin(FDStationBaseAngle)
                               - FDCameraHalfDistance*cos(FDStationBaseAngle) );
  positionFrameCenter[0].setY( positionFDStation.y()
			       - FDCameraR*cos(FDStationBaseAngle) 
			       - FDCameraHalfDistance*sin(FDStationBaseAngle) );
  positionFrameCenter[0].setZ( 0.0 );


  positionFrameCenter[3].setX( FDCameraHalfDistance );
  positionFrameCenter[3].setY( positionFDStation.y() - FDCameraR2 );
  positionFrameCenter[3].setZ( 0.0 );

  positionFrameCenter[2].setX( -FDCameraHalfDistance );
  positionFrameCenter[2].setY( positionFDStation.y() - FDCameraR2 );
  positionFrameCenter[2].setZ( 0.0 );


  positionFrameCenter[5].setX( - FDCameraR*sin(FDStationBaseAngle)
                               + FDCameraHalfDistance*cos(FDStationBaseAngle) );
  positionFrameCenter[5].setY( positionFDStation.y()
                               - FDCameraR*cos(FDStationBaseAngle)
                               - FDCameraHalfDistance*sin(FDStationBaseAngle) );
  positionFrameCenter[5].setZ( 0.0 );


  positionFrameCenter[4].setX( - FDCameraR*sin(FDStationBaseAngle)
                               - FDCameraHalfDistance*cos(FDStationBaseAngle) );
  positionFrameCenter[4].setY( positionFDStation.y()
			       - FDCameraR*cos(FDStationBaseAngle) 
			       + FDCameraHalfDistance*sin(FDStationBaseAngle) );
  positionFrameCenter[4].setZ( 0.0 );

  // -----------------------------------------------------
  // Rotation Angle List
  // -----------------------------------------------------
  RotationAngle[2]  = RotationAngle[3]  =  27.0*RADI;
  RotationAngle[0]  = RotationAngle[1]  =  45.0*RADI;  
  RotationAngle[6]  = RotationAngle[7]  =  -9.0*RADI;  
  RotationAngle[4]  = RotationAngle[5]  =   9.0*RADI;
  RotationAngle[10] = RotationAngle[11] = -45.0*RADI;
  RotationAngle[8]  = RotationAngle[9]  = -27.0*RADI;

  // -----------------------------------------------------
  // Mirror Geometory
  // -----------------------------------------------------
  InnerRadiusMirror        =  6.067;  
  DepthMirror              =  11.0*0.001;
  DiameterMirror           =  1.730;  
  OuterRadiusMirror        =  InnerRadiusMirror+DepthMirror;

  MirrorPhi                = asin(DiameterMirror/InnerRadiusMirror);
  DeltaMirrorPhi           = 2.0*MirrorPhi;
  
  for( int i=0; i<6; i++ ){
  DistanceFromCameraFrame[2*i]   = 2.9153;
  DistanceFromCameraFrame[2*i+1] = 2.9153 + 0.3685;

  SlopeAngleMirror[2*i]          =  25.5*RADI;
  SlopeAngleMirror[2*i+1]        =  10.5*RADI;

  HeigthMirror[2*i]              =  1.5267;
  HeigthMirror[2*i+1]            =  5.4994;
  }

  for( int i=0; i<12; i++ ){
    double PHIX(0.);
    double PHIY(0.);
    int j(0);
    if( i==2  || i==3  ) { PHIX = -sin(27.0*RADI); PHIY = cos(27.0*RADI); j=1; }
    if( i==0  || i==1  ) { PHIX = -cos(45.0*RADI); PHIY = sin(45.0*RADI); j=0; }
    if( i==6  || i==7  ) { PHIX = sin(9.0*RADI);   PHIY = cos(9.0*RADI);  j=3; }
    if( i==4  || i==5  ) { PHIX = -sin(9.0*RADI);  PHIY = cos(9.0*RADI);  j=2; }
    if( i==10 || i==11 ) { PHIX = cos(45.0*RADI);  PHIY = sin(45.0*RADI); j=5; }
    if( i==8  || i==9  ) { PHIX = sin(27.0*RADI);  PHIY = cos(27.0*RADI); j=4; }
    positionMirror[i].setX( positionFrameCenter[j].x() + DistanceFromCameraFrame[i]*PHIX );
    positionMirror[i].setY( positionFrameCenter[j].y() + DistanceFromCameraFrame[i]*PHIY );
    positionMirror[i].setZ( HeigthMirror[i] );
  }

  for( int i=0; i<12 ; i++ ){
    double PHIX(0.);
    double PHIY(0.);
    double RTMPXY(InnerRadiusMirror*cos(SlopeAngleMirror[i]));
    double RTMPZ(InnerRadiusMirror*sin(SlopeAngleMirror[i])); 
    if( i==2  || i==3  ) { PHIX = sin(27.0*RADI);  PHIY = -cos(27.0*RADI); }
    if( i==0  || i==1  ) { PHIX = cos(45.0*RADI);  PHIY = -sin(45.0*RADI); }
    if( i==6  || i==7  ) { PHIX = -sin(9.0*RADI);  PHIY = -cos(9.0*RADI);  }
    if( i==4  || i==5  ) { PHIX = sin(9.0*RADI);   PHIY = -cos(9.0*RADI);  }
    if( i==10 || i==11 ) { PHIX = -cos(45.0*RADI); PHIY = -sin(45.0*RADI); }
    if( i==8  || i==9  ) { PHIX = -sin(27.0*RADI); PHIY = -cos(27.0*RADI); }
    positionMirrorCenter[i].setX( positionMirror[i].x() + RTMPXY*PHIX );
    positionMirrorCenter[i].setY( positionMirror[i].y() + RTMPXY*PHIY );
    positionMirrorCenter[i].setZ( positionMirror[i].z() + RTMPZ ); 
  }

  CameraAngle              = 7.79*RADI;

  for( int i=0; i<12; i++ ){
        PosiMirrorCX[i]=positionMirrorCenter[i].x();
	PosiMirrorCY[i]=positionMirrorCenter[i].y();
	PosiMirrorCZ[i]=positionMirrorCenter[i].z();
  }

  // -----------------------------------------------------
  //  Segment Mirror Geometory 
  // -----------------------------------------------------
  SizeSegMirror = 0.330;
  
  double R0 = 5.900;
  double RM = InnerRadiusMirror;

  double MR1 = 0.700;
  double MR2 = 0.700/tan(30.0*RADI);
  double MR3 = 1.380;

  double X0  = R0*sin(MR3/RM);
  double X1  = R0*sin(MR1/RM);
  double X4  = (R0/RM)*MR1*cos(60.0*RADI);
  double X8  = (R0/RM)*MR3*cos(60.0*RADI);
  double X14 = (R0/RM)*MR2*cos(30.0*RADI);
  double TMPPOSIX[18]    = {  X0,  X1, -X1, -X0,
                              X4, -X4,  X4, -X4,
                              X8, -X8,  X8, -X8,
                              0.0, 0.0,
                              X14, -X14, X14, -X14 };
  
  double Y0  = R0*cos(MR3/RM);
  double Y1  = R0*cos(MR1/RM);
  double Y4  = R0*sqrt(1.0-(MR1*sin(60.0*RADI)/RM)*(MR1*sin(60.0*RADI)/RM));
  double Y8  = R0*sqrt(1.0-(MR3*sin(60.0*RADI)/RM)*(MR3*sin(60.0*RADI)/RM));
  double Y12 = R0*sqrt(1.0-(MR2/RM)*(MR2/RM));
  double TMPPOSIY[18]    = { Y0, Y1, Y1,Y0,
                             Y4, Y4, Y4, Y4,
                             Y8, Y8, Y8, Y8,
                             Y12, Y12, Y12, 
                             Y12, Y12, Y12 };

  double Z0  = 0.0;
  double Z4  = (R0/RM)*MR1*sin(60.0*RADI);
  double Z8  = (R0/RM)*MR3*sin(60.0*RADI);
  double Z12 = (R0/RM)*MR2;
  double Z14 = (R0/RM)*MR2*sin(30.0*RADI);
  double TMPPOSIZ[18]    = {  Z0, Z0, Z0, Z0,
                                Z4, Z4, -Z4, -Z4,
                                Z8, Z8, -Z8, -Z8,
                                Z12, -Z12,
                                Z14, Z14, -Z14, -Z14 };

  double RX0  = 0.0;
  double RX4  = asin( MR1*sin(60.0*RADI)/RM);
  double RX8  = asin( MR3*sin(60.0*RADI)/RM);
  double RX12 = asin( MR2/RM);
  double RX14 = asin( MR2*sin(30.0*RADI)/RM);
  double TMPROTATIONX[18] = { RX0, RX0, RX0, RX0,
                                RX4, RX4, -RX4, -RX4,
                                RX8, RX8, -RX8, -RX8,
                                RX12, -RX12,
                                RX14, RX14, -RX14, -RX14 };

  double TMPROTATIONY[18] = { 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0 };
 
  double RZ0  = asin(MR3/RM);
  double RZ1  = asin(MR1/RM);
  double RZ4  = asin(MR1*cos(60.0*RADI)/RM);
  double RZ8  = asin(MR3*cos(60.0*RADI)/RM);
  double RZ12 = asin( MR2*cos(30.0*RADI)/RM);
  double TMPROTATIONZ[18] = { -RZ0, -RZ1,  RZ1,  RZ0,
                              -RZ4,  RZ4, -RZ4,  RZ4,
                              -RZ8,  RZ8, -RZ8,  RZ8,
 			       0.0, 0.0,
                              -RZ12, RZ12, -RZ12, RZ12 };

  for( int i=0; i<18; i++ ){
    SegMirrorPosX[i]=TMPPOSIX[i];
    SegMirrorPosY[i]=TMPPOSIY[i];
    SegMirrorPosZ[i]=TMPPOSIZ[i];           

    SegMirrorRotX[i]=TMPROTATIONX[i];
    SegMirrorRotY[i]=TMPROTATIONY[i];
    SegMirrorRotZ[i]=TMPROTATIONZ[i];
  }

  // -----------------------------------------------------
  // Camera Frame Geometory
  // -----------------------------------------------------
  FrameMainWidth  = 0.05; 
  FrameMainHeight = 6.745;
  FrameMainDepth  = 0.4;

  FrameTubeDiameter = 5.0;
  FrameTubeLength   = 1.270;
 
  for( int i=0; i<6; i++ ){
    double PHIX(0.),  PHIY(0.);
    if( i==1 ){ PHIX=cos(27.0*RADI); PHIY=sin(27.0*RADI); }
    if( i==0 ){ PHIX=cos(45.0*RADI); PHIY=sin(45.0*RADI); }
    if( i==3 ){ PHIX=cos(9.0*RADI);  PHIY=-sin(9.0*RADI); }
    if( i==2 ){ PHIX=cos(9.0*RADI);  PHIY=sin(9.0*RADI);  } 
    if( i==5 ){ PHIX=cos(45.0*RADI); PHIY=-sin(45.0*RADI); }
    if( i==4 ){ PHIX=cos(27.0*RADI); PHIY=-sin(27.0*RADI); }

    positionFrameMain[i][0].setX( positionFrameCenter[i].x() + FrameTubeLength*0.5*PHIX );
    positionFrameMain[i][0].setY( positionFrameCenter[i].y() + FrameTubeLength*0.5*PHIY );
    positionFrameMain[i][0].setZ( FrameMainHeight/2.0 ); 

    positionFrameMain[i][1].setX( positionFrameCenter[i].x() - FrameTubeLength*0.5*PHIX );
    positionFrameMain[i][1].setY( positionFrameCenter[i].y() - FrameTubeLength*0.5*PHIY );
    positionFrameMain[i][1].setZ( FrameMainHeight/2.0 );
  }
  
  // -----------------------------------------------------
  // Paraglass Geometory
  // -----------------------------------------------------
  ParaglassWidth  = 1.080;
  ParaglassHeight = 0.930;
  ParaglassThick  = 0.005;
  DistanceFromCBSurface = 0.005;

  // -----------------------------------------------------
  // Camera Box Geometory
  // -----------------------------------------------------
  CameraBoxWidth  = 1.270;    // Remark !! True Box Size = 1160 
  CameraBoxHeight = 1.010;
  CameraBoxThick  = 0.040;
  CameraBoxThickW = (CameraBoxWidth-ParaglassWidth)/2.0;
  CameraBoxDepth  = 0.250;
  //CameraDistance  = 3.0;
  CameraDistance  = 3.0 - 0.001*30.5;   

  CameraBoxLX[0]  = CameraBoxWidth;
  CameraBoxLY[0]  = CameraBoxThick;
  CameraBoxLZ[0]  = CameraBoxHeight-CameraBoxThick*2.0;

  CameraBoxLX[1]  = CameraBoxThickW;
  CameraBoxLY[1]  = CameraBoxDepth-CameraBoxThick*2.0;
  CameraBoxLZ[1]  = CameraBoxHeight-CameraBoxThick*2.0;

  CameraBoxLX[2]  = CameraBoxLX[1];
  CameraBoxLY[2]  = CameraBoxLY[1];
  CameraBoxLZ[2]  = CameraBoxLZ[1];

  CameraBoxLX[3]  = CameraBoxWidth;
  CameraBoxLY[3]  = CameraBoxDepth;
  CameraBoxLZ[3]  = CameraBoxThick;

  CameraBoxLX[4]  = CameraBoxLX[3];
  CameraBoxLY[4]  = CameraBoxLY[3];
  CameraBoxLZ[4]  = CameraBoxLZ[3];

  for( int i=0; i<12 ; i++ ){
    double PHIX(0.);
    double PHIY(0.);

    double RTMPXY( ( CameraBoxDepth*0.5 + CameraDistance  )*cos(SlopeAngleMirror[i])  );
    double RTMPZ(  ( CameraBoxDepth*0.5 + CameraDistance  )*sin(SlopeAngleMirror[i])  ); 
    double RTMPXYP(( DistanceFromCBSurface + ParaglassThick*0.5 + CameraDistance )
                       *cos(SlopeAngleMirror[i]) );
    double RTMPZP( ( DistanceFromCBSurface + ParaglassThick*0.5 + CameraDistance )
		       *sin(SlopeAngleMirror[i]) );

    if( i==2  || i==3  ) { PHIX = sin(27.0*RADI);  PHIY = -cos(27.0*RADI); }
    if( i==0  || i==1  ) { PHIX = cos(45.0*RADI);  PHIY = -sin(45.0*RADI); }
    if( i==6  || i==7  ) { PHIX = -sin(9.0*RADI);  PHIY = -cos(9.0*RADI);  }
    if( i==4  || i==5  ) { PHIX = sin(9.0*RADI);   PHIY = -cos(9.0*RADI);  }
    if( i==10 || i==11 ) { PHIX = -cos(45.0*RADI); PHIY = -sin(45.0*RADI); }
    if( i==8  || i==9  ) { PHIX = -sin(27.0*RADI); PHIY = -cos(27.0*RADI); }

   positionCameraBoxCenter[i].setX( positionMirror[i].x() + RTMPXY*PHIX );
   positionCameraBoxCenter[i].setY( positionMirror[i].y() + RTMPXY*PHIY );
   positionCameraBoxCenter[i].setZ( positionMirror[i].z() + RTMPZ ); 

   positionParaglassCenter[i].setX( positionMirror[i].x() + RTMPXYP*PHIX );
   positionParaglassCenter[i].setY( positionMirror[i].y() + RTMPXYP*PHIY );
   positionParaglassCenter[i].setZ( positionMirror[i].z() + RTMPZP );
  }

  for( int i=0; i<12 ; i++ ){
    double PHIX(0.);
    double PHIY(0.);

    double RTMPXY( (CameraBoxDepth+CameraBoxThick)*0.5*cos(SlopeAngleMirror[i])  );
    double RTMPZ(  (CameraBoxDepth+CameraBoxThick)*0.5*sin(SlopeAngleMirror[i])  ); 

    if( i==2  || i==3  ) { PHIX = sin(27.0*RADI);  PHIY = -cos(27.0*RADI); }
    if( i==0  || i==1  ) { PHIX = cos(45.0*RADI);  PHIY = -sin(45.0*RADI); }
    if( i==6  || i==7  ) { PHIX = -sin(9.0*RADI);  PHIY = -cos(9.0*RADI); }
    if( i==4  || i==5  ) { PHIX = sin(9.0*RADI);   PHIY = -cos(9.0*RADI); }
    if( i==10 || i==11 ) { PHIX = -cos(45.0*RADI); PHIY = -sin(45.0*RADI); }
    if( i==8  || i==9  ) { PHIX = -sin(27.0*RADI); PHIY = -cos(27.0*RADI); }

    positionCameraBox[i][0].setX( positionCameraBoxCenter[i].x() + RTMPXY*PHIX );
    positionCameraBox[i][0].setY( positionCameraBoxCenter[i].y() + RTMPXY*PHIY );
    positionCameraBox[i][0].setZ( positionCameraBoxCenter[i].z() + RTMPZ );
  }
   
  for( int i=0; i<12; i++ ){
    double PHIX(0.),  PHIY(0.);
    double XWIDTH = ( CameraBoxWidth - CameraBoxThickW )*0.5;

    double XTMP   = positionCameraBoxCenter[i].x();
    double YTMP   = positionCameraBoxCenter[i].y();
    double ZTMP   = positionCameraBoxCenter[i].z();

    if( i==2  || i==3  ){ PHIX=cos(27.0*RADI); PHIY=sin(27.0*RADI); }
    if( i==0  || i==1  ){ PHIX=cos(45.0*RADI); PHIY=sin(45.0*RADI); }
    if( i==6  || i==7  ){ PHIX=cos(9.0*RADI);  PHIY=-sin(9.0*RADI); }
    if( i==4  || i==5  ){ PHIX=cos(9.0*RADI);  PHIY=sin(9.0*RADI);  } 
    if( i==10 || i==11 ){ PHIX=cos(45.0*RADI); PHIY=-sin(45.0*RADI); }
    if( i==8  || i==9  ){ PHIX=cos(27.0*RADI); PHIY=-sin(27.0*RADI); }

    positionCameraBox[i][1].setX( XTMP + XWIDTH*PHIX );
    positionCameraBox[i][1].setY( YTMP + XWIDTH*PHIY );
    positionCameraBox[i][1].setZ( ZTMP );

    positionCameraBox[i][2].setX( XTMP - XWIDTH*PHIX );
    positionCameraBox[i][2].setY( YTMP - XWIDTH*PHIY );
    positionCameraBox[i][2].setZ( ZTMP );
  }

  for( int i=0; i<12 ; i++ ){
    double PHIX(0.);
    double PHIY(0.);
    double RTMPXY( (CameraBoxHeight-CameraBoxThick)*0.5*sin(SlopeAngleMirror[i])  );
    double RTMPZ(  (CameraBoxHeight-CameraBoxThick)*0.5*cos(SlopeAngleMirror[i])  ); 

    if( i==2  || i==3  ) { PHIX = -sin(27.0*RADI); PHIY = cos(27.0*RADI); }
    if( i==0  || i==1  ) { PHIX = -cos(45.0*RADI); PHIY = sin(45.0*RADI); }
    if( i==6  || i==7  ) { PHIX = sin(9.0*RADI);   PHIY = cos(9.0*RADI);  }
    if( i==4  || i==5  ) { PHIX = -sin(9.0*RADI);  PHIY = cos(9.0*RADI);  }
    if( i==10 || i==11 ) { PHIX = cos(45.0*RADI);  PHIY = sin(45.0*RADI); }
    if( i==8  || i==9  ) { PHIX = sin(27.0*RADI);  PHIY = cos(27.0*RADI); }

    positionCameraBox[i][3].setX( positionCameraBoxCenter[i].x() + RTMPXY*PHIX );
    positionCameraBox[i][3].setY( positionCameraBoxCenter[i].y() + RTMPXY*PHIY );
    positionCameraBox[i][3].setZ( positionCameraBoxCenter[i].z() + RTMPZ );

    positionCameraBox[i][4].setX( positionCameraBoxCenter[i].x() - RTMPXY*PHIX );
    positionCameraBox[i][4].setY( positionCameraBoxCenter[i].y() - RTMPXY*PHIY );
    positionCameraBox[i][4].setZ( positionCameraBoxCenter[i].z() - RTMPZ );
  }

  // -----------------------------------------------------
  // PMT Geometory
  // -----------------------------------------------------
  NoPMT = 256;

  DistanceFromCBSurfacePMT = 0.0305;
  DepthSensitivePMT        = 0.001;

  NoPMTSide = 6;
  NoPMTZPlanes = 6;
   
  irPMT[0] = irPMT[1] = irPMT[2] = irPMT[3] = irPMT[4] = irPMT[5] = 0.0;

  //orPMT[0] = 30.25*0.001/cos(30.0*RADI);
  orPMT[0] = 30.25*0.001;
  orPMT[1] = orPMT[0];
  orPMT[2] = 0.05;
  orPMT[3] = 0.05;
  orPMT[4] = 0.044;
  orPMT[5] = 0.044;

  zPMT[0]  = 0.0;
  zPMT[1]  = 0.02;
  zPMT[2]  = 0.04;
  zPMT[3]  = 0.100;
  zPMT[4]  = 0.100;
  zPMT[5]  = 0.124;

  NoPMTZsPlanes=3;
  irPMTs[0] = irPMTs[1] = irPMTs[2] = 0.0;

  //orPMTs[0] = 29.25*0.001/cos(30.0*RADI);  
  orPMTs[0] = 30.25*0.001;
  orPMTs[1] = orPMTs[0];
  orPMTs[2] = orPMTs[0];

  zPMTs[0]  = 0.0;
  zPMTs[1]  = DepthSensitivePMT*0.5;
  zPMTs[2]  = DepthSensitivePMT;

  //PMTShiftH = 0.4805;
  PMTShiftH = -0.4805;
  PMTShiftV =  0.4026;

  PMTSpacingH = 0.062;
  PMTSpacingV = 0.05369;

  // -----------------------------------------------------
  // BG3 Geometory
  // -----------------------------------------------------
  NoBG3        = NoPMT;
  NoBG3Side    = NoPMTSide;
  NoBG3ZPlanes = 3;
  BG3Depth     = 0.004;     
  Space_BG3toPMT = 0.0001;
                                                                                      
  irBG3[0]     = 0.0;
  irBG3[1]     = 0.0;
  irBG3[2]     = 0.0;

  orBG3[0]     = orPMT[0];
  orBG3[1]     = orPMT[0];
  orBG3[2]     = orPMT[0];
 
  zBG3[0]      = 0.0;
  zBG3[1]      = BG3Depth*0.5;
  zBG3[2]      = BG3Depth;

  BG3ShiftH    = PMTShiftH;
  BG3ShiftV    = PMTShiftV; 
                                                                                       
  BG3SpacingH  = PMTSpacingH;
  BG3SpacingV  = PMTSpacingV;

  for( int i=0; i<12 ; i++ ){
    double PHIX(0.);
    double PHIY(0.);

    double DistanceSurface = CameraDistance;
    double DistanceBias = DistanceFromCBSurfacePMT+CameraDistance;

    double RTMPXYSurface( DistanceSurface*cos(SlopeAngleMirror[i])  );
    double RTMPZSurface(  DistanceSurface*sin(SlopeAngleMirror[i])  );

    double RTMPXY( ( DistanceBias + DepthSensitivePMT )*cos(SlopeAngleMirror[i])  );
    double RTMPZ(  ( DistanceBias + DepthSensitivePMT )*sin(SlopeAngleMirror[i])  ); 

    double RTMPXYS( ( DistanceBias )*cos(SlopeAngleMirror[i])  );
    double RTMPZS(  ( DistanceBias )*sin(SlopeAngleMirror[i])  );

    // Added Space_BG3toPMT ( 0.0001 m )  2006.10.25 by T.Shibata
    double RTMPXYBG3( ( DistanceBias - BG3Depth - Space_BG3toPMT )*cos(SlopeAngleMirror[i])  );
    double RTMPZBG3(  ( DistanceBias - BG3Depth - Space_BG3toPMT )*sin(SlopeAngleMirror[i])  );

    if( i==2  || i==3  ) { PHIX = sin(27.0*RADI);  PHIY = -cos(27.0*RADI); }
    if( i==0  || i==1  ) { PHIX = cos(45.0*RADI);  PHIY = -sin(45.0*RADI); }
    if( i==6  || i==7  ) { PHIX = -sin(9.0*RADI);  PHIY = -cos(9.0*RADI);  }
    if( i==4  || i==5  ) { PHIX = sin(9.0*RADI);   PHIY = -cos(9.0*RADI);  }
    if( i==10 || i==11 ) { PHIX = -cos(45.0*RADI); PHIY = -sin(45.0*RADI); }
    if( i==8  || i==9  ) { PHIX = -sin(27.0*RADI); PHIY = -cos(27.0*RADI); }

   positionCenterPMTSurface[i].setX( positionMirror[i].x() + RTMPXYSurface*PHIX ); 
   positionCenterPMTSurface[i].setY( positionMirror[i].y() + RTMPXYSurface*PHIY );
   positionCenterPMTSurface[i].setZ( positionMirror[i].z() + RTMPZSurface*PHIY );

   positionCenterPMTs[i].setX( positionMirror[i].x() + RTMPXY*PHIX );
   positionCenterPMTs[i].setY( positionMirror[i].y() + RTMPXY*PHIY );
   positionCenterPMTs[i].setZ( positionMirror[i].z() + RTMPZ ); 

   positionCenterPMTSs[i].setX( positionMirror[i].x() + RTMPXYS*PHIX );
   positionCenterPMTSs[i].setY( positionMirror[i].y() + RTMPXYS*PHIY );
   positionCenterPMTSs[i].setZ( positionMirror[i].z() + RTMPZS );

   positionCenterBG3s[i].setX( positionMirror[i].x() + RTMPXYBG3*PHIX );
   positionCenterBG3s[i].setY( positionMirror[i].y() + RTMPXYBG3*PHIY );
   positionCenterBG3s[i].setZ( positionMirror[i].z() + RTMPZBG3 );
  }			       

  for( int CID=0; CID<12 ; CID++ ){

     double PHIX(0.), PHIX2(0.);
     double PHIY(0.), PHIY2(0.);
     if( CID==2  || CID==3  ) { PHIX = -sin(27.0*RADI); PHIY = cos(27.0*RADI); 
                                PHIX2 = cos(27.0*RADI); PHIY2 = sin(27.0*RADI); }
     if( CID==0  || CID==1  ) { PHIX = -cos(45.0*RADI); PHIY = sin(45.0*RADI); 
                                PHIX2 = cos(45.0*RADI); PHIY2 = sin(45.0*RADI); }
     if( CID==6  || CID==7  ) { PHIX = sin(9.0*RADI);   PHIY = cos(9.0*RADI);  
                                PHIX2 = cos(9.0*RADI);  PHIY2 = -sin(9.0*RADI); }
     if( CID==4  || CID==5  ) { PHIX = -sin(9.0*RADI);  PHIY = cos(9.0*RADI);  
                                PHIX2 = cos(9.0*RADI);  PHIY2 = sin(9.0*RADI);  }
     if( CID==10 || CID==11 ) { PHIX = cos(45.0*RADI);  PHIY = sin(45.0*RADI); 
                                PHIX2 = cos(45.0*RADI); PHIY2 = -sin(45.0*RADI);}
     if( CID==8  || CID==9  ) { PHIX = sin(27.0*RADI);  PHIY = cos(27.0*RADI); 
                                PHIX2 = cos(27.0*RADI); PHIY2 = -sin(27.0*RADI);}

     for( int copyNo=0; copyNo<NoPMT; copyNo++ ){
       int     i=copyNo/16;      // 0 - 15 : along X-axis
       int     j=copyNo-i*16;    // 0 - 15 : along Z-axis
       double  dH;
       if( j/2*2 == j ){ dH = 0.;
       }else{            dH = PMTSpacingH/2.0; }

       //double ShiftH = PMTShiftH - i*PMTSpacingH - dH;
       double ShiftH = PMTShiftH + i*PMTSpacingH + dH;  // <-- changed by T.Shibata 2007.01.01
       double ShiftV = PMTShiftV - j*PMTSpacingV;

       double NewCenterX=positionCenterPMTs[CID].x() + ShiftV*sin(SlopeAngleMirror[CID])*PHIX;
       double NewCenterY=positionCenterPMTs[CID].y() + ShiftV*sin(SlopeAngleMirror[CID])*PHIY;  

       positionPMT[CID][copyNo].setX( NewCenterX + ShiftH*PHIX2 );
       positionPMT[CID][copyNo].setY( NewCenterY + ShiftH*PHIY2 );
       positionPMT[CID][copyNo].setZ( positionCenterPMTs[CID].z() + ShiftV*cos(SlopeAngleMirror[CID]) );
     
       double NewCenterXS=positionCenterPMTSs[CID].x() + ShiftV*sin(SlopeAngleMirror[CID])*PHIX;
       double NewCenterYS=positionCenterPMTSs[CID].y() + ShiftV*sin(SlopeAngleMirror[CID])*PHIY;

       positionPMTSen[CID][copyNo].setX( NewCenterXS + ShiftH*PHIX2 );
       positionPMTSen[CID][copyNo].setY( NewCenterYS + ShiftH*PHIY2 );
       positionPMTSen[CID][copyNo].setZ( positionCenterPMTSs[CID].z() + ShiftV*cos(SlopeAngleMirror[CID]) );

       double NewCenterXBG3=positionCenterBG3s[CID].x() + ShiftV*sin(SlopeAngleMirror[CID])*PHIX;
       double NewCenterYBG3=positionCenterBG3s[CID].y() + ShiftV*sin(SlopeAngleMirror[CID])*PHIY;

       positionBG3[CID][copyNo].setX( NewCenterXBG3 + ShiftH*PHIX2 );
       positionBG3[CID][copyNo].setY( NewCenterYBG3 + ShiftH*PHIY2 );
       positionBG3[CID][copyNo].setZ( positionCenterBG3s[CID].z() + ShiftV*cos(SlopeAngleMirror[CID]) );
     }

  } 

}
//-------------------------------------------------------------------------------------------
FDParameters::~FDParameters(){}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetFDStation(void){ return positionFDStation; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetFDStationWall( int i ){ return positionFDStationWall[i]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetFDStationRoof( int i ){ return positionFDStationRoof[i]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetMirrorCenter(int CID){ return positionMirrorCenter[CID];}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetMirror(int CID){ return positionMirror[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetSegMirrorPos( int SegID ){
  CLHEP::Hep3Vector tmpVector;
       tmpVector.setX(SegMirrorPosX[SegID]);
       tmpVector.setY(SegMirrorPosY[SegID]);
       tmpVector.setZ(SegMirrorPosZ[SegID]);
  return tmpVector;
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetSegMirrorRot( int SegID ){
  CLHEP::Hep3Vector tmpVector;
       tmpVector.setX(SegMirrorRotX[SegID]);
       tmpVector.setY(SegMirrorRotY[SegID]);
       tmpVector.setZ(SegMirrorRotZ[SegID]);
  return tmpVector;
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetFrameCenter( int FID ){ return positionFrameCenter[FID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetFrameMain( int SID, int RL ){ return positionFrameMain[SID][RL]; } 
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCameraBoxCenter( int CID ){ return positionCameraBoxCenter[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCameraBox( int CID, int PID ){ return positionCameraBox[CID][PID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetParaglassCenter( int CID ){ return positionParaglassCenter[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCenterPMTSurcafe( int CID ){ return positionCenterPMTSurface[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCenterPMTs( int CID ){ return positionCenterPMTs[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetPMT( int CID, int PMTID ){ return positionPMT[CID][PMTID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCenterPMTSs( int CID ){ return positionCenterPMTSs[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetPMTSen( int CID, int PMTID ){ return positionPMTSen[CID][PMTID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetCenterBG3s( int CID ){ return positionCenterBG3s[CID]; }
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::GetBG3( int CID, int PMTID ){ return positionBG3[CID][PMTID]; }
//-------------------------------------------------------------------------------------------
void FDParameters::DumpPosition(void)
{ 
  ofstream dump;
  dump.open("./dump.dat"); 

  /*
  dump << GetFDStation().x() << " " << GetFDStation().y() << " " << GetFDStation().z() << endl;
  */

  /*
  for( int i=0; i<12; i++ )
    dump << i<< " " << GetMirror(i).x() << " " << GetMirror(i).y() << " " << GetMirror(i).z() << endl;
  */

  /*
  for( int i=0; i<12; i++ )
    dump << GetMirrorCenter(i).x() << " " << GetMirrorCenter(i).y() << " " << GetMirrorCenter(i).z() << endl;  

  for( int i=0; i<6; i++ )
     dump << GetFrameCenter(i).x() << " " << GetFrameCenter(i).y() << " " << GetFrameCenter(i).z() << endl;

  for( int i=0; i<6; i++ ){
    for( int j=0; j<2; j++ ) 
      dump << GetFrameMain(i,j).x() << " " << GetFrameMain(i,j).y() << " " << GetFrameMain(i,j).z() << endl;
  }

  for( int i=0; i<12; i++ )
     dump << GetCameraBoxCenter(i).x() << " " << GetCameraBoxCenter(i).y() << " " << GetCameraBoxCenter(i).z() << endl;
 
  for( int i=0; i<12; i++ ){
    for( int j=0; j<5; j++ )
      dump << GetCameraBox(i,j).x() << " " << GetCameraBox(i,j).y() << " " << GetCameraBox(i,j).z() << endl; 
  }

  for( int i=0; i<12; i++ )
    dump << GetParaglassCenter(i).x() << " " << GetParaglassCenter(i).y() << " " << GetParaglassCenter(i).z() << endl;

  */
  for( int i=0; i<12; i++ )
    dump << GetCenterPMTs(i).x() << " " << GetCenterPMTs(i).y() << " " << GetCenterPMTs(i).z() << endl;
  /*

  for( int i=0; i<12; i++ ){
    for( int j=0; j<256; j++ )
      dump << GetPMT(i,j).x() << " " << GetPMT(i,j).y() << " " << GetPMT(i,j).z() << endl;
  }

  for( int i=0; i<12; i++ )
    dump << GetCenterPMTSs(i).x() << " " << GetCenterPMTSs(i).y() << " " << GetCenterPMTSs(i).z() << endl;
  */

  /*
  for( int i=0; i<12; i++ ){
    for( int j=0; j<256; j++ )
      dump << i << " " << j << " " << GetPMTSen(i,j).x() << " " << GetPMTSen(i,j).y() << " " << GetPMTSen(i,j).z() << endl;
  } 
  */

  /*
  for( int i=0; i<12; i++ )
    dump << GetCenterBG3s(i).x() << " " << GetCenterBG3s(i).y() << " " << GetCenterBG3s(i).z() << endl;

  */

  /*
  for( int i=0; i<12; i++ ){
    for( int j=0; j<256; j++ )
      dump << GetBG3(4,j).x() << " " << GetBG3(4,j).y() << " " << GetBG3(4,j).z() << endl;
  }
  */

  dump.close();

 return;
}
//-------------------------------------------------------------------------------------------
bool FDParameters::Cameraframe( CLHEP::Hep3Vector InputPosition, int MirrorID )
{

  if( InputPosition.z() <  HeigthMirror[MirrorID] ) return false;

  double tangentz 
    = (InputPosition.z()-HeigthMirror[MirrorID])/(InputPosition - positionMirror[MirrorID]).perp();

  CLHEP::Hep3Vector rmirror = positionMirrorCenter[MirrorID] - positionMirror[MirrorID];
  CLHEP::Hep3Vector rinput  = InputPosition                  - positionMirror[MirrorID];

  double cosrphi
    = (rmirror.x()*rinput.x() + rmirror.y()*rinput.y())/(rmirror.perp()*rinput.perp());

  if( tan(SlopeAngleMirror[MirrorID]-CameraAngle) < tangentz      &&
      tangentz < tan(SlopeAngleMirror[MirrorID]+CameraAngle)      &&
      cosrphi > cos(CameraAngle) ) return true; 

  return false;
}
//-------------------------------------------------------------------------------------------
int FDParameters:: PMTfame( CLHEP::Hep3Vector InputPosition, int MirrorID )
{
  int rPMT(-1);
  if( Cameraframe(InputPosition, MirrorID ) ){     
     
      rPMT = 1;
  }else{
      rPMT =-1;
  }

  return rPMT;
}
//-------------------------------------------------------------------------------------------
bool FDParameters::MirrorHit( int MID, CLHEP::Hep3Vector HitMirrorPoint )
{

  double delta = MirrorPhi; 

  CLHEP::Hep3Vector rc = positionMirror[MID];
  CLHEP::Hep3Vector r0 = positionMirrorCenter[MID];
  CLHEP::Hep3Vector ht = HitMirrorPoint*0.001;  // mm --> m

  CLHEP::Hep3Vector r1 = rc - r0;
  CLHEP::Hep3Vector r2 = ht - r0;

  r1 = r1.unit();
  r2 = r2.unit();
 
  double crossAngle = r1.dot(r2);
         crossAngle = acos(crossAngle);

  if( crossAngle < delta ) return true;
 
  return false;
}
//-------------------------------------------------------------------------------------------
double FDParameters::RelativeHitDistance( int CID, CLHEP::Hep3Vector HitPoint )
{
  CLHEP::Hep3Vector r0 = GetCenterPMTs(CID);
  CLHEP::Hep3Vector ht = HitPoint*0.001 - r0;  // mm --> m  
  return ht.mag();
}
//-------------------------------------------------------------------------------------------
CLHEP::Hep3Vector FDParameters::RelativeHitPoint( int CID, int PID, 
                                                  CLHEP::Hep3Vector HitPoint )
{
  CLHEP::Hep3Vector r0 = GetCenterPMTSs(CID);
  //CLHEP::Hep3Vector r0 = GetPMTSen(CID, PID);
  CLHEP::Hep3Vector ht  = HitPoint*0.001 - r0;  // mm --> m  
  CLHEP::Hep3Vector RelativePosition = ht;

  RelativePosition.rotateZ(-RotationAngle[CID]);
  RelativePosition.rotateX(SlopeAngleMirror[CID]);

  return RelativePosition;    
}
//-------------------------------------------------------------------------------------------
bool FDParameters::FieldOfViewofCamera( int CID, CLHEP::Hep3Vector CreationPoint,
					CLHEP::Hep3Vector ShiftMirror )
{

  double HorizontalAngleTheta(-100.);
  double VerticalAnglePsi(-100.);

  CLHEP::Hep3Vector PhotonDirectionfromMirror
    = ( CreationPoint - ( GetMirror(CID) + ShiftMirror ) ).unit();
  CLHEP::Hep3Vector MirrorOC
    = ( GetMirrorCenter(CID) - GetMirror(CID) ).unit();

  // Horizontal Angle(degree)
  double costh = MirrorOC.x()*PhotonDirectionfromMirror.x()
    + MirrorOC.y()*PhotonDirectionfromMirror.y();

  HorizontalAngleTheta=acos(costh/(PhotonDirectionfromMirror.rho()*MirrorOC.rho()))/RADI;

  double judge
    =MirrorOC.x()*PhotonDirectionfromMirror.y()-MirrorOC.y()*PhotonDirectionfromMirror.x();
  if( judge > 0 ) HorizontalAngleTheta*=(-1);

  // Vertical Angle(degree)
  VerticalAnglePsi=asin(PhotonDirectionfromMirror.z()/PhotonDirectionfromMirror.mag())/RADI;
  VerticalAnglePsi-=(SlopeAngleMirror[CID]/RADI);

  // if( abs(HorizontalAngleTheta)<=9 && abs(VerticalAnglePsi)<=8 ) return true;
  if( abs(HorizontalAngleTheta)<=10 && abs(VerticalAnglePsi)<=8 ) return true;

  return false;
}
//-------------------------------------------------------------------------------------------



