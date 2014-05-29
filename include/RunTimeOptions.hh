#ifndef RunTimeOptions_hh
#define RunTimeOptions_hh

//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//
// Runtime flags and useful constants.  This will be updated  //
// periodically as more functionality is added.               //
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=//

//#include "G4String.hh"

//------------------------------//
// Runtime options
//------------------------------//

enum RunOption
{
  RO_NULL       = 0,
  RO_Tracking   = 1<<0,
  RO_Stepping   = 1<<1
};

//------------------------------//
// Quick Particle Reference
//------------------------------//
enum PDGID
{
  PDG_E     = 11,
  PDG_nuE   = 12,
  PDG_MU    = 13,
  PDG_nuMU  = 14,
  PDG_TAU   = 15,
  PDG_nuTAU = 16,
  PDG_GAMMA = 22
};

#endif
