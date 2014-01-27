#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

#include "ELSParameters.hh"

#include "ReadFile.hh"

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsList: public G4VUserPhysicsList
{
public:
  PhysicsList(Plist *plist);
  ~PhysicsList();

protected:
  //void ConstructParticle();
  void ConstructParticle();
  void ConstructProcess();

  void SetCuts();

protected:
  void ConstructHadrons();
  void ConstructBosons();
  void ConstructLeptons();

protected:
  void ConstructEM();
  void ConstructMyEM(); // added in 2012.05.10 
  void ConstructEM_Dmitri(); // added in 2013.04.23

  // Added in 2011.12.31 , G4.9.4p02
  void ConstructLEEM();

  void ConstructHD();

protected:  // added in 2012.09.24
  void ConstructOp();
  void ConstructOp_OLD();  

private: 

  // Added in 2011.12.31 , G4.9.4p02
  G4VPhysicsConstructor*  emPhysicsList;

  int Flag_PhotoNuclear_Process;  
  int Flag_Cerenkov_Process;

  double cutlenght_particle;
};
#endif

  
