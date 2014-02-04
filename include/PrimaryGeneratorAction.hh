#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <CLHEP/Vector/ThreeVector.h>

#include "PrimaryEbeam.hh"
#include "BeamWaveform.hh"       

#include "G4VUserPrimaryGeneratorAction.hh"
#include "ReadFile.hh"
                                                                                             
class Plist;
class Beam;
class Waveform;
class DetectorConstruction;
class G4ParticleGun;
class G4Event;

// Primary Particle ----
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*,
                         Plist *plist,
                         char* Fname,
			 G4int Nparticles);
  ~PrimaryGeneratorAction();

  Ebeam *ebeam;
  Waveform *waveform; 

public:
  void GeneratePrimaries(G4Event*);

private:

  ofstream EventNumberFile;  
  ofstream OutputFileBeamEnergy;

  G4ParticleGun* particleGun;
  DetectorConstruction* myDetector;

  // Beam Energy Parameter
  int beam_energy_mode;

  int beam_particle;

  int beam_waveform_flag;

    // constant energy and gaussian function
    double primary_energy[4];
    double constant_energy;
    double gaussian_mean;
    double gaussian_sigma; 

    //crystal ball function
    double crystalball_mean[4];
    double crystalball_sigma[4];
    double crystalball_alpha[4];
    double crystalball_n[4];
    double crystalball_xmin;
    double crystalball_xmax;
  
    // Read Generator File
    ifstream BeamEnergyFile;
  
  // Beam Emmitance mode 
  int beam_emmitance_mode;

  // Beam Position Parameter
  double injection_position_x[4];
  double injection_position_y[4];
  double injection_position_z[4];

  double injection_transverse_phi0x[4];
  double injection_transverse_phi0y[4];

  double transverse_phix;
  double transverse_phiy;
  double injection_transverse_momentum_x;
  double injection_transverse_momentum_y;
  double injection_transverse_momentum_z;

  double injection_position_shift_totalx;
  double injection_position_shift_totaly;
  double injection_position_shift_x[4];
  double injection_position_sigma_x[4];
  double injection_position_shift_y[4];
  double injection_position_sigma_y[4];

  double injection_position_sigma[4];
  double injection_position_delta_h;
  double injection_position_delta_v;

  double injection_position_r;
  double injection_position_theta;

  // Beam Direction
  double injection_direction;
  double injection_direction_x[4];
  double injection_direction_y[4];
  double injection_direction_z[4];

  double injection_beamemmitance_g[2][4];
  double injection_beamemmitance_b[2][4];
  double injection_beamemmitance_a[2][4];
  double injection_beamemmitance_e[2][4];

  double injection_beamemmitance_h[4];
  double injection_beamemmitance_v[4];

  // Time 
  double particle_time;

};
#endif
