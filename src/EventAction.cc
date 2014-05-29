#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

//------------------------------------------------//
// Constructor
//------------------------------------------------//
EventAction::EventAction(std::ofstream* trkOut,
			 std::ofstream* stepOut) :
  m_trkOut(NULL),
  m_stepOut(NULL)
{
  
  m_trkOut = trkOut;
  m_stepOut = stepOut;
  
}

//------------------------------------------------//
// Destructor
//------------------------------------------------//
EventAction::~EventAction()
{

  m_trkOut  = NULL;
  m_stepOut = NULL;
}

//------------------------------------------------//
// Begin
//------------------------------------------------//
void EventAction::BeginOfEventAction(const G4Event* evt)
{

  G4int event_id = evt->GetEventID();
  if(m_trkOut)  (*m_trkOut)  << "Event " << event_id << G4endl;
  if(m_stepOut) (*m_stepOut) << "Event " << event_id << G4endl;

}

//------------------------------------------------//
// End
//------------------------------------------------//
void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4int event_id = evt->GetEventID();

  // periodic printing
  if( event_id % 50 == 0 )
    G4cout << "Event# = " << event_id << " :  " << G4endl;


  // Save that event is over
  if(m_trkOut)  (*m_trkOut)  << "End " << event_id << G4endl;
  if(m_stepOut) (*m_stepOut) << "End " << event_id << G4endl;
}
