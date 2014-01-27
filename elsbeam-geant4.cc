////////////////////////////////////////////////////////
//    LINAC CONFIGURATION --- BEAM BACKGROUND RADIATION
//                            Version 2
//       Creation : 2009/11/30 by T.ShiBbata
//       Last Updata : 2010/11/08 by T.Shibata
//       Last Updata : 2011/12/13 by T.Shibata
////////////////////////////////////////////////////////
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"
 
#include "ReadFile.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIExecutive.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#ifdef G4UI_USE_XM
#include "G4UIXm.hh"
#endif

#ifdef G4UI_USE_WIN32
#include "G4UIWin32.hh"
#endif

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

//----------------------------------------------------------------------
int main(int argc,char** argv) {

  // Set the random seed to the timer                                                                                                       
  int randseed = (int) time(NULL);
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine(randseed));

  //---------------------------
  // GEANT EVENT SET 
  int visualization = 0;
  int number_of_Events = 1;
  char* GeneratorFilename="gene.dat";
  char* OutputFilename0="detector-plane.dat";
  char* OutputFilename1="generated.dat";
  char* OutputFilename2="dedx.dat";
  char* OutputFilename3="faradaycup.dat";
  char* OutputFilename4="energydeposit.dat";
  char* mattOutput="matt_test.dat";
  char* SetupFilename="setup-file.dat";
  //---------------------------

  int c(-1);
  while(1){
    c = getopt(argc, argv, "d:i:o:g:k:f:e:s:v:m:r");
    if(c==-1) break;
    switch(c){
    case 'd':
      number_of_Events  = atoi(optarg);   
      break;
    case 'i':
      GeneratorFilename = optarg;
      break;
    case 'o':
      OutputFilename0   = optarg;   // detector-plane.dat
      break;
    case 'g':
      OutputFilename1   = optarg;   // generated.dat
      break;
    case 'k':
      OutputFilename2   = optarg;   // dedx.dat
      break;
    case 'f':
      OutputFilename3   = optarg;   // faradaycup.dat
      break;
    case 'e':
      OutputFilename4   = optarg;   // energydeposit.dat
      break;
    case 'm':
      mattOutput = optarg;
      break;
    case 's':
      SetupFilename = optarg;
      break;
    case 'v':
      visualization = atoi(optarg);
      break;
    case 'r':
      break;
    default:
      exit(1);
    }
  }

  //--------------------------------
  // Read Setup File
  Plist *plist;
  plist = new Plist();
  plist->readfile(SetupFilename);
  //--------------------------------

  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose( plist,
                                                       OutputFilename0,
                                                       OutputFilename1,
                                                       OutputFilename2,
                                                       OutputFilename3,
                                                       OutputFilename4,
						       mattOutput));

  // Run manager
  G4RunManager * runManager = new G4RunManager;

  // UserInitialization classes (mandatory)
  DetectorConstruction* FCdetector = new DetectorConstruction;
  FCdetector->setup_detector_parameter(plist);

  runManager->SetUserInitialization(FCdetector);
  //runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new PhysicsList(plist));

  //========================
  G4UIsession* session=0;

#if defined(G4UI_USE_XM)
      session = new G4UIXm(argc,argv);
#elif defined(G4UI_USE_WIN32)
      session = new G4UIWin32();
#elif defined(G4UI_USE_TCSH)
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif

#ifdef G4VIS_USE
      // Visualization, if you choose to have it!
      G4VisManager* visManager = new G4VisExecutive;       
      visManager->Initialize();
#endif

  //========================

  // UserAction classes
  runManager->SetUserAction(new PrimaryGeneratorAction(FCdetector,
                                                       plist,
                                                       1, //number_of_Events,
                                                       GeneratorFilename));

  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);
  runManager->SetUserAction(new SteppingAction);

  //Initialize G4 kernel
  runManager->Initialize();

  // get the pointer to the User Interface manager
  G4UImanager * UI = G4UImanager::GetUIpointer();

  //=================================================

  //=======================================

  UI->ApplyCommand("/run/verbose 1");
  UI->ApplyCommand("/tracking/verbose 1");
  UI->ApplyCommand("/event/verbose 1");

  //======================================
  //if( visualization == 1 ){
  //#ifdef G4VIS_USE
    //delete visManager;
  //#endif
    //  delete runManager;
    //}
  //======================================
  
  // BEAM ON -- EXECUTE RUN --
  if( visualization != 1 ){
   runManager->BeamOn(number_of_Events);
   delete runManager;
  }
  else{
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    ui->SessionStart();
    delete UI;
  }

  return 0;
}
//----------------------------------------------------------------------
