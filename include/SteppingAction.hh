#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include <fstream>

class SteppingAction : public G4UserSteppingAction
{
public:
  
  // Constructor
  SteppingAction(std::ofstream* output);

  // Destructor
  ~SteppingAction(){
    m_output = NULL;
  };

  // Action!
  void UserSteppingAction(const G4Step*);
  
private:
  std::ofstream* m_output;

};
#endif

