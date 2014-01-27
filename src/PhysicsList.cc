#include <stdlib.h>

#include "globals.hh"

#include "PhysicsList.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"

// added in 2012.09.24
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

// added in 2013.05.07
#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"
#include "G4SystemOfUnits.hh"

#include "G4PhotoNuclearProcess.hh"
#include "G4GammaNuclearReaction.hh"
            
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4LENeutronInelastic.hh"
#include "G4NeutronInelasticProcess.hh"

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPCapture.hh"
          
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScatteringData.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPCaptureData.hh"

//#include "G4CrossSectionDataStore.hh"

#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"

// Comment out in 2011.12.30 because of update g4.9.4p02
//  #include "G4MultipleScattering.hh" 

#include "G4CoulombScattering.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4UrbanMscModel93.hh"

// --- added in 2013.04.23 for Dmitri-san's EM
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"

//-----------------------
// Geant4.9.2.p02
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"

// Low Energy Physics Process
#include "G4EmLivermorePhysics.hh"
// Livemore Data Library 
  // Photo-electric effect (class G4LivermorePhotoElectricModel)
  // Compton scattering (class G4LivermoreComptonModel)
  // Gamma conversion (also called pair production, class G4LivermoreGammaConversionModel)
  // Nuclear gamma conversion (class G4LivermoreNuclearGammaConversionModel)

  // Bremsstrahlung (class G4LivermoreBremsstrahlungModel)
  // Ionisation and delta ray production (class G4LivermoreIonisationModel)

// ICRU73 Based Ion Model
// Ionisation and delta ray production (class G4IonParametrisedLossModel)
//  G4ionIonisation* ionIoni = new G4ionIonisation();
//  ionIoni -> SetEmModel(new G4IonParametrisedLossModel());

#include "G4EmPenelopePhysics.hh"
// Penelope MC code
  // Compton scattering (class G4PenelopeComptonModel)
  // Rayleigh scattering (class G4PenelopeRayleighModel)
  // Gamma conversion (also called pair production, class GPenelopeGammaConversionModel)
  // Photo-electric effect (class G4PenelopePhotoElectricModel)

  // Bremsstrahlung (class G4PenelopeBremsstrahlungModel)
  // Ionisation and delta ray production (class G4PenelopeIonisationModel)
  // Positron annihilation (class class G4PenelopeAnnihilationModel)
#include "G4EmDNAPhysics.hh"


#include "G4ProductionCutsTable.hh"

//-----------------------

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

//-- added in 2012.05.10 
#include "G4EmProcessOptions.hh"

#include "G4MscStepLimitType.hh"
#include "G4KleinNishinaModel.hh"
//---


//PhysicsList::PhysicsList():  G4VUserPhysicsList()
PhysicsList::PhysicsList(Plist *plist):  G4VUserPhysicsList()
{

  defaultCutValue = 0.1*cm;  // changed by T.Shibata 2010.11.08
  SetVerboseLevel(1);
 
  Flag_PhotoNuclear_Process = 0;
  Flag_PhotoNuclear_Process=plist->Photo_nuclear_process_flag();

  Flag_Cerenkov_Process = 0;
  Flag_Cerenkov_Process=plist->Cerenkov_process_flag();

  cutlenght_particle = plist->Cutlenght_particle();
  defaultCutValue = cutlenght_particle*mm;

  // Dmitri-san defined it as 1e-4 * mm

}
PhysicsList::~PhysicsList()
{
  
}

void PhysicsList::ConstructParticle()
{
  ConstructLeptons();
  ConstructBosons();
  //ConstructHadrons();
}

void PhysicsList::ConstructHadrons()
{
  // Proton
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  
  // Neutron  
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // Meson
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  // Ions
  G4Alpha::Alpha(); 
  G4Deuteron::DeuteronDefinition();
  G4GenericIon::GenericIonDefinition();
  G4He3::He3Definition();
  G4Triton::TritonDefinition();

}
void PhysicsList::ConstructBosons()
{
  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();

}
void PhysicsList::ConstructLeptons()
{
  // leptons : e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}
void PhysicsList::ConstructProcess()
{
  AddTransportation();
  
  // Standrd Processcd 
  ConstructEM(); 
  // ConstructMyEM();
  // ConstructEM_Dmitri();

  // Lowenergy Process
  //ConstructLEEM();

  // ConstructHD();  

  ConstructOp();

}

void PhysicsList::ConstructLEEM()
{

  //emPhysicsList =  new G4EmStandardPhysics();   // default 

  // emPhysicsList =  new G4EmStandardPhysics_option1();  // HEP fast, not precise
  //                                                         minimal MSC step limitaion a=0.8 for e+-
  // emPhysicsList =  new G4EmStandardPhysics_option2();  // Experimental ?
  //                                                         
  // emPhysicsList =  new G4EmStandardPhysics_option3();  // Medical, Spcaes

  emPhysicsList = new G4EmLivermorePhysics();
  
  //emPhysicsList = new G4EmDNAPhysics();
  
  //emPhysicsList = new G4EmPenelopePhysics(); 

  /* added in 2012.10.04 */
  //G4EmProcessOptions emOptions;
  //emOptions.SetDeexcitationActiveRegion("World",true, true, true );

  //emOptions.SetDeexcitaionActive(true);
  //emOptions.SetAugerActive(true);
  //emOptions.SetPIXEActive(true);

  //emOptions.SetFluo(true);
  //emOptions.SetAuger(true);
  //emOptions.SetPIXE(true);

  emPhysicsList->ConstructProcess(); 
}

void PhysicsList::ConstructMyEM() // added in 2012.05.10 
{ 
 
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add standard EM Processes                                                                                  
  //                                                                                                            
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {

      ////ph->RegisterProcess(new G4RayleighScattering, particle);                                              
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      G4ComptonScattering* cs   = new G4ComptonScattering;
      cs->SetModel(new G4KleinNishinaModel());
      ph->RegisterProcess(cs, particle);
      ph->RegisterProcess(new G4GammaConversion, particle);

    } else if (particleName == "e-") {

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
   
      G4eIonisation* eIoni = new G4eIonisation();

      // eIoni->SetStepFunction(0.2, 1.0*mm);     // same as standars process 
      // eIoni->SetStepFunction(0.2, 100.0*um);      // same as lowenergy process
         eIoni->SetStepFunction(0.2, 1.0*mm);   

      ph->RegisterProcess(eIoni, particle);

      ph->RegisterProcess(new G4eBremsstrahlung(), particle);

    } else if (particleName == "e+") {

      ph->RegisterProcess(new G4eMultipleScattering(), particle);
 
      G4eIonisation* eIoni = new G4eIonisation();

      // eIoni->SetStepFunction(0.2, 1.0*mm);  // same as standars process 
      // eIoni->SetStepFunction(0.2, 100.0*um);   // same as lowenergy process
      eIoni->SetStepFunction(0.2, 1.0*mm);   

      ph->RegisterProcess(eIoni, particle);

      ph->RegisterProcess(new G4eBremsstrahlung(), particle);

      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    }

  }

    G4EmProcessOptions emOptions;
    // multiple coulomb scattering 
    // emOptions.SetMscStepLimitation(fMinimal);  // same as standars process 
    emOptions.SetMscStepLimitation(fUseSafety);   // same as standars process 
    // emOptions.SetMscStepLimitation(fUseDistanceToBoundary); // same as lowenergy process 

    // try 01 & 02 & 12
    emOptions.SetMscRangeFactor(0.04);
    emOptions.SetMscGeomFactor(2.5);
    emOptions.SetSkin(3.0);    

}

void PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);

      // Photo-Nuclear Process added by 2009.04.07-08
      // This process is Inelastics photon process
      // This model generates the final state for gamma-nuclear inelastic scattering using the CHIPS model. 
      //  For details, see the chapter on the Chiral Invariant Phase Space Model (CHIPS) 
      //  in the Geant4 Physics Reference Manual. 
      if( Flag_PhotoNuclear_Process ){
        G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess();

        G4GammaNuclearReaction* GDNmodel = new G4GammaNuclearReaction();

        thePhotoNuclearProcess->RegisterMe(GDNmodel);
        pmanager->AddDiscreteProcess(thePhotoNuclearProcess);
      }
      //

    } else if (particleName == "e-") {
      //electron

      // EM Standard 
      // Comment out in 2011.12.30 because of update g4.9.4p02
      // pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      // G4eMultipleSacattering uses G4UrbanMscModel2 ( 93 ? )
      // G4WentzelVIModel, G4CoulombScattering                 
      // ... But 
      //     G4EmLivermorePhysics G4eMultipleScattering uses G4GoudsmithSaundersonMscModel for e+/e-    

      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);

      /* use as default */      
      //pmanager->AddProcess(new G4eBremsstrahlung,   -1, -1,3);

    } else if (particleName == "e+") {
      // positron

      // EM Standard
      // Comment out in 2011.12.30 because of update g4.9.4p02
      //pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4eIonisation,         -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);

      /* use as default */      
      //pmanager->AddProcess(new G4eBremsstrahlung,   -1, -1,3);    
      pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 4);

    } else {

       // EM Standard
       // Comment out in 2011.12.30 because of update g4.9.4p02
       // pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
         pmanager->AddProcess(new G4eMultipleScattering,-1, 1,1);

    }

  }

  // Comment out = try 07
  // use fUseDistanceToBoundary = try 06
  // use fUseDistanceToBoundary = try 10
  // G4EmProcessOptions emOptions;
  // emOptions.SetMscStepLimitation(fUseSafety);  //default  
  // emOptions.SetMscStepLimitation(fUseDistanceToBoundary); 

}

void PhysicsList::ConstructEM_Dmitri()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  theParticleIterator->reset();
  while ((*theParticleIterator)())
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4String particleName = particle->GetParticleName();

      if (particleName == "gamma")
        {
          // gamma                                                                                                                      
          ph->RegisterProcess(new G4PhotoElectricEffect, particle);
          ph->RegisterProcess(new G4ComptonScattering, particle);
          ph->RegisterProcess(new G4GammaConversion, particle);

        }
      else if (particleName == "e-")
        {
          //electron                                                                                                                    
          ph->RegisterProcess(new G4eMultipleScattering, particle);
          ph->RegisterProcess(new G4eIonisation, particle);
          ph->RegisterProcess(new G4eBremsstrahlung, particle);

	}
      else if (particleName == "e+")
	{
          //positron                                                                                                                    
          ph->RegisterProcess(new G4eMultipleScattering, particle);
          ph->RegisterProcess(new G4eIonisation, particle);
          ph->RegisterProcess(new G4eBremsstrahlung, particle);
          ph->RegisterProcess(new G4eplusAnnihilation, particle);

	}
      else if (particleName == "mu+" || particleName == "mu-")
	{
          //muon                                                                                                                        
          ph->RegisterProcess(new G4MuMultipleScattering, particle);
          ph->RegisterProcess(new G4MuIonisation, particle);
          ph->RegisterProcess(new G4MuBremsstrahlung, particle);
          ph->RegisterProcess(new G4MuPairProduction, particle);

        }
      else if (particleName == "proton" || particleName == "pi-" || particleName == "pi+")
        {
          //proton                                                                                                                      
          ph->RegisterProcess(new G4hMultipleScattering, particle);
          ph->RegisterProcess(new G4hIonisation, particle);
          ph->RegisterProcess(new G4hBremsstrahlung, particle);
          ph->RegisterProcess(new G4hPairProduction, particle);

	}
      else if (particleName == "alpha" || particleName == "He3")
        {
          //alpha                                                                                                                       
          ph->RegisterProcess(new G4hMultipleScattering, particle);
          ph->RegisterProcess(new G4ionIonisation, particle);

	}
      else if (particleName == "GenericIon")
        {
          //Ions                                                                                                                        
          ph->RegisterProcess(new G4hMultipleScattering, particle);
          ph->RegisterProcess(new G4ionIonisation, particle);

        }
      else if ((!particle->IsShortLived()) && (particle->GetPDGCharge() != 0.0)
	       && (particle->GetParticleName() != "chargedgeantino"))
        {
          //all others charged particles except geantino                                                                                
          ph->RegisterProcess(new G4hMultipleScattering, particle);
          ph->RegisterProcess(new G4hIonisation, particle);
        }
    }
}



void PhysicsList::ConstructHD()
{
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if ( particleName == "proton" ){

       G4HadronElasticProcess* protonElastic = new G4HadronElasticProcess();    
       G4LElastic *theProtonElasticModel = new G4LElastic();
       protonElastic->RegisterMe(theProtonElasticModel);  

       pmanager->AddDiscreteProcess( protonElastic ); // Add Process
 
    } else if ( particleName == "neutron" ) {     

      // Neutron Elastic Scattering Process and Model < 20MeV
      G4HadronElasticProcess* neutronElastic           = new G4HadronElasticProcess();    

      /* Comment out by T.Shibata 2010.02.08
      G4LElastic *theElasticModel = new G4LElastic(); // Model for Elastic Scattering  
      neutronElastic->RegisterMe(theElasticModel);  //Registration of Model
      */

      // Neutron HP w/ Thermal Neutron
      // Cross Section 
      G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData(); // Cross Section Data
      neutronElastic->AddDataSet(theHPElasticData);   // Registration of Cross Section Data 

      G4NeutronHPThermalScatteringData* theHPTHermalScatteringData 
                                    	= new G4NeutronHPThermalScatteringData();
      neutronElastic->AddDataSet(theHPTHermalScatteringData);

      // Model w/ Thermal Neutron
      G4NeutronHPElastic* theHPElasticModel  = new G4NeutronHPElastic(); // Neutron Elastic Model   
      theHPElasticModel->SetMinEnergy(4.0*eV);      // Set Minimum Energy
      neutronElastic->RegisterMe(theHPElasticModel);  // Registration of Neutron Elastic Model

      G4NeutronHPThermalScattering* theNeutronTHermalElasticModel 
                                             = new G4NeutronHPThermalScattering();  
      theNeutronTHermalElasticModel->SetMaxEnergy(4.0*eV); 
      neutronElastic->RegisterMe(theNeutronTHermalElasticModel);
  
      pmanager->AddDiscreteProcess( neutronElastic ); // Add Process      
       
      // pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);

      // Neutron Inelastic Proces
      G4NeutronInelasticProcess* neutronInelastic = new G4NeutronInelasticProcess();

      // Added by T.Shibata 2010.02.12 
      //-------
      G4LENeutronInelastic* theLEInelasticModel = new G4LENeutronInelastic();
      theLEInelasticModel->SetMinEnergy(19.0*MeV);
      neutronInelastic->RegisterMe(theLEInelasticModel);
      //-------
      
      G4NeutronHPInelastic* theHPInelasticModel = new G4NeutronHPInelastic();
      neutronInelastic->RegisterMe(theHPInelasticModel);

      G4NeutronHPInelasticData* theHPInelasticData = new G4NeutronHPInelasticData();
      neutronInelastic->AddDataSet(theHPInelasticData);

      pmanager->AddDiscreteProcess( neutronInelastic );

      // Neutron Fission Process
      G4HadronFissionProcess* neutronFission = new G4HadronFissionProcess();

      G4LFission* theLFissionModel = new G4LFission();
      theLFissionModel->SetMaxEnergy(20.*TeV);
      neutronFission->RegisterMe(theLFissionModel);
      /*
      G4NeutronHPFission* theHPFissionModel = new G4NeutronHPFission();
      G4NeutronHPFissionData* theHPFissionData = new G4NeutronHPFissionData();
      neutronFission->RegisterMe(theHPFissionModel); 
      neutronFission->AddDataSet(theHPFissionData); 
      */
      pmanager->AddDiscreteProcess( neutronFission );

      // Neutron Capture Process
      G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
      G4LCapture* theLCaptureModel = new G4LCapture();
      theLCaptureModel->SetMaxEnergy(20.*TeV);
      neutronCapture->RegisterMe(theLCaptureModel);
      /*
      G4NeutronHPCapture* theHPCaptureModel = new G4NeutronHPCapture();
      G4NeutronHPCaptureData* theHPCaptureData = new G4NeutronHPCaptureData();
      neutronCapture->RegisterMe(theHPCaptureModel);
      neutronCapture->AddDataSet(theHPCaptureData);
      */

      pmanager->AddDiscreteProcess( neutronCapture );
 
    }else if( particleName == "alpha" ||
              particleName == "GenericIon" ||
              particleName == "deuteron" ||
              particleName == "triton" ){
      pmanager->AddProcess( new G4hMultipleScattering, -1, 1, 1 );
      pmanager->AddProcess( new G4ionIonisation, -1, 2, 2 );

    }

  }

}

void PhysicsList::ConstructOp()
{

  if( !Flag_Cerenkov_Process ) return;

  G4Cerenkov* theCerenkovProcess =  new G4Cerenkov("Cerenkov");
  //G4Scintillation* theScintillationProcess = new G4Scintillation("Scintillation");

  theCerenkovProcess->SetMaxNumPhotonsPerStep(20);
  theCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  theCerenkovProcess->SetTrackSecondariesFirst(true);

  //theScintillationProcess->SetScintillationYieldFactor(1.);
  //theScintillationProcess->SetTrackSecondariesFirst(true);

  G4OpRayleigh*    theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpMieHG*         theMieScatteringProcess = new G4OpMieHG();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (theCerenkovProcess->IsApplicable(*particle)) {
        pmanager->AddProcess(theCerenkovProcess);
        pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
	if (particleName == "opticalphoton") {
   	    pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	    pmanager->AddDiscreteProcess(theMieScatteringProcess);
	}
	/*
    if (theScintillationProcess->IsApplicable(*particle)) {
        pmanager->AddProcess(theScintillationProcess);
        pmanager->SetProcessOrderingToLast(theScintillationProcess, idxAtRest);
        pmanager->SetProcessOrderingToLast(theScintillationProcess, idxPostStep);
    }
	*/
  }

}

// added in 2012.09.24
void PhysicsList::ConstructOp_OLD()
{
  /*
  G4OpAbsorption*  theAbsorptionProcess         = new G4OpAbsorption();
  G4OpRayleigh*    theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpBoundaryProcess* theBoundaryProcess       = new G4OpBoundaryProcess();
  */

  G4int MaxNumPhotons(300);
  G4Cerenkov* theCerenkovProcess = new G4Cerenkov("Cerenkov");

  theCerenkovProcess->SetTrackSecondariesFirst(true);
  theCerenkovProcess->SetMaxNumPhotonsPerStep(MaxNumPhotons);

  /*
  G4OpticalSurfaceModel themodel = unified;
  theBoundaryProcess->SetModel(themodel);
  */

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){

    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if( theCerenkovProcess->IsApplicable(*particle) ){
      //pmanager->AddContinuousProcess(theCerenkovProcess);
      // changed in 2013.05.07
      pmanager->AddProcess(theCerenkovProcess);
      pmanager->SetProcessOrdering( theCerenkovProcess, idxPostStep );
    }

    /*
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(theAbsorptionProcess);
      pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(theBoundaryProcess);
    }
    */

  }
  
}

void PhysicsList::SetCuts()
{

  SetCutsWithDefault();
  if (verboseLevel>0) DumpCutValuesTable();

  //-----------------------
  // Geant4.9.2.p02
  // -->Geant4.9.5
  //-----------------------
  /* 250eV and 990eV have same results */
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(250*eV, 1*GeV);  // not be used 

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(990.0*eV, 1*GeV); // <-- we can use this

  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(500*eV, 100*GeV); // <-- we can use this
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 100*GeV); // <-- we can use this
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(50*eV, 100*GeV); // <-- we can use this
  
  //G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(10*eV, 100*GeV); // <-- we can use this

  G4cout << "CUT VALUES: (below)" << G4endl;
  DumpCutValuesTable();
  G4cout << "CUT VALUES: (above)" << G4endl;

  /*
  G4cout << 
    " GetLowEnergy  = " << G4ProductionCutsTable::GetProductionCutsTable()->GetLowEdgeEnergy()/eV  <<
    " GetHighEnergy = " << G4ProductionCutsTable::GetProductionCutsTable()->GetHighEdgeEnergy()/eV << G4endl;
  */
  //-----------------------

}
