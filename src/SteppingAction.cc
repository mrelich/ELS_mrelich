#include "SteppingAction.hh"
#include "G4SteppingManager.hh"


//-----------------------------------//
// Constructor
//-----------------------------------//
SteppingAction::SteppingAction(std::ofstream* output) :
  m_output(NULL)
{ 
  
  m_output = output;

}

//-----------------------------------//
// Action!
//-----------------------------------//
void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  // Decide quantities want to write to output file
  // Save PDG, pre-step x,y,z, dX, E, dE loss, 
  // dE loss from ionization only, trk ID
  G4Track* aTrack = aStep->GetTrack();
  int pdg     = aTrack->GetParticleDefinition()->GetPDGEncoding();
  float x     = aStep->GetPreStepPoint()->GetPosition().x()/cm;
  float y     = aStep->GetPreStepPoint()->GetPosition().y()/cm;
  float z     = aStep->GetPreStepPoint()->GetPosition().z()/cm;
  float dX    = aStep->GetStepLength()/cm;
  float E     = aTrack->GetKineticEnergy()/MeV;
  float dE    = aStep->GetTotalEnergyDeposit()/MeV;
  float dEion = dE - aStep->GetNonIonizingEnergyDeposit()/MeV;
  int trkId   = aTrack->GetTrackID();

  // Write quantites
  (*m_output) <<
    pdg   << " " <<
    x     << " " << 
    y     << " " <<
    z     << " " <<
    dX    << " " <<
    E     << " " <<
    dE    << " " <<
    dEion << " " <<
    trkId << " " <<
    G4endl;
  

}

