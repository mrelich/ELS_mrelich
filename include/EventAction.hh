#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include <fstream>

class G4Event;
// Event Action !!
class EventAction : public G4UserEventAction
{

 public:
  EventAction(std::ofstream* trkOut, std::ofstream* stepOut);
  ~EventAction();
  
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
 private:
  
  std::ofstream* m_trkOut;
  std::ofstream* m_stepOut;
  
};
#endif
