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

//----------------------------------------------------------------------//
// Help Menu
//----------------------------------------------------------------------//
void help()
{
  
  cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  cout<<endl;
  cout<<"Options: "<<endl;
  cout<<"\t-ne <int>"<<endl;
  cout<<"\t\tSets the number of events"<<endl;
  cout<<"\t-np <int>"<<endl;
  cout<<"\t\tSets the number of particles"<<endl;
  cout<<"\t-i <filename>"<<endl;
  cout<<"\t\tSpecify the generator filename"<<endl;
  cout<<"\t-s <filename>"<<endl;
  cout<<"\t\tSpecify the setup filename"<<endl;
  cout<<"\t-v <int>"<<endl;
  cout<<"\t\tSpecify visulization on or off (0 off, 1 on)"<<endl;
  cout<<"\t-e"<<endl;
  cout<<"\t\tTurn on energy dump option"<<endl;
  cout<<"\t-i"<<endl;
  cout<<"\t\tTurn on quick check method"<<endl;
  cout<<"\t-h"<<endl;
  cout<<"\t\tPrint this menu"<<endl;
  cout<<endl;
  cout<<"-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-="<<endl;
  cout<<endl;
}

//----------------------------------------------------------------------//
// Main
//----------------------------------------------------------------------//
int main(int argc,char** argv) {

  // Set the random seed to the timer                                                                                                       
  int randseed = (int) time(NULL);
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine(randseed));

  //---------------------------
  // GEANT EVENT SET 
  int visualization = 0;
  int number_of_Events = 1;
  int number_of_Particles = 1; // # generated particles
  char* GeneratorFilename="gene.dat";
  char* SetupFilename="setup-file.dat";
  RunOption runOpt = RO_NULL;
  
  //---------------------------
  for(int i=1; i<argc; ++i){
    if( strcmp(argv[i], "-ne") == 0 )
      number_of_Events = atoi( argv[++i] );
    else if( strcmp(argv[i], "-np") == 0 )
      number_of_Particles = atoi( argv[++i] );
    else if( strcmp(argv[i], "-i") == 0 )
      GeneratorFilename = argv[++i];
    else if( strcmp(argv[i], "-s") == 0 )
      SetupFilename = argv[++i];
    else if( strcmp(argv[i], "-v") == 0 )
      visualization = atoi(argv[++i]);
    else if( strcmp(argv[i], "-e") == 0 )
      runOpt = (RunOption) (runOpt | RO_EnergyDump);
    else if( strcmp(argv[i], "-q") == 0 )
      runOpt = (RunOption) (runOpt | RO_QuickCheck);
    else{
      help();
      return 0;
    }
  }//end loop over arguments

  
  //--------------------------------
  // Read Setup File
  Plist *plist;
  plist = new Plist();
  plist->readfile(SetupFilename);
  //--------------------------------
  
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose( plist, runOpt ) );

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
                                                       GeneratorFilename,
						       number_of_Particles));

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

