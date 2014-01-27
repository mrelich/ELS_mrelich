#include <iomanip>

#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
               
#include "PhysicsParameters.hh"
        
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* myDC,
                                               Plist *plist,
                                               G4int Nevents,
                                               char* Fname )
  :myDetector(myDC)  
{
  EventNumberFile.open("eventnumber.dat");
  OutputFileBeamEnergy.open("initialbeamenergy.dat");

  // Electron Beam Object
  ebeam = new Ebeam();
  waveform = new Waveform();

  //-----------------------------------------------------
  // Set up Parameters

  //Beam Energy Mode 
  beam_energy_mode = plist->beamenergy_mode();    

  beam_particle    = plist->beam_particle();
 
  // Read waveform 
  beam_waveform_flag = plist->ElectronBeam_Waveform_flag();
  if( beam_waveform_flag == 1 ) waveform->ReadWaveformFile(plist);
  
  // Constrant mode and Gaussian Mode
  for( int i=0; i<4; i++ ){
    primary_energy[i]=plist->primaryenergy(i);
  }
  constant_energy = primary_energy[0];
  gaussian_mean   = primary_energy[0];
  gaussian_sigma  = primary_energy[1];
   
  // Crystal Ball Fucntion Mode
  for( int i=0; i<4; i++ ){
    crystalball_mean[i]  = plist->primaryenergy_crystalball_mean(i);
    crystalball_sigma[i] = plist->primaryenergy_crystalball_sigma(i);
    crystalball_alpha[i] = plist->primaryenergy_crystalball_alpha(i);
    crystalball_n[i]     = plist->primaryenergy_crystalball_n(i);
  }
  crystalball_xmin       = plist->primaryenergy_crystalball_xmin(0);
  crystalball_xmax       = plist->primaryenergy_crystalball_xmax(0);

  // Generator File 
  BeamEnergyFile.open(Fname);  

  // Beam Emmitance Mode 
  beam_emmitance_mode = plist->beam_emmitance_mode();    

  // Beam Position and Direction 
  for( int i=0; i<4; i++ ){
    injection_position_x[i]=plist->injection_positionx(i);
    injection_position_y[i]=plist->injection_positiony(i);
    injection_position_z[i]=plist->injection_positionz(i);
  
    injection_position_sigma[i]=plist->injection_positionsigma(i); 

    injection_position_shift_x[i]=plist->injection_position_shiftx(i);
    injection_position_sigma_x[i]=plist->injection_position_sigmax(i);
    injection_position_shift_y[i]=plist->injection_position_shifty(i);
    injection_position_sigma_y[i]=plist->injection_position_sigmay(i);

    injection_direction_x[i]=plist->injection_directionx(i);
    injection_direction_y[i]=plist->injection_directiony(i);
    injection_direction_z[i]=plist->injection_directionz(i);

    injection_transverse_phi0x[i]=plist->injection_transverse_phi0x(i);  
    injection_transverse_phi0y[i]=plist->injection_transverse_phi0y(i);  
  }

  injection_direction=sqrt( injection_direction_x[0]*injection_direction_x[0]
                            +injection_direction_y[0]*injection_direction_y[0]
			                +injection_direction_z[0]*injection_direction_z[0] );

  injection_direction_x[0] = injection_direction_x[0]/injection_direction;
  injection_direction_y[0] = injection_direction_y[0]/injection_direction;
  injection_direction_z[0] = injection_direction_z[0]/injection_direction;

  for( int i=0; i<2; i++ ){ 
   for( int j=0; j<4; j++ ){ 
     injection_beamemmitance_g[i][j] = plist->injection_beamemmitance_g(i,j);
     injection_beamemmitance_b[i][j] = plist->injection_beamemmitance_b(i,j);
     injection_beamemmitance_a[i][j] = plist->injection_beamemmitance_a(i,j);
     injection_beamemmitance_e[i][j] = plist->injection_beamemmitance_e(i,j);
   }
  }
  injection_beamemmitance_h[0] = injection_beamemmitance_g[0][0];
  injection_beamemmitance_h[1] = injection_beamemmitance_b[0][0];
  injection_beamemmitance_h[2] = injection_beamemmitance_a[0][0];
  injection_beamemmitance_h[3] = injection_beamemmitance_e[0][0];

  injection_beamemmitance_v[0] = injection_beamemmitance_g[1][0];
  injection_beamemmitance_v[1] = injection_beamemmitance_b[1][0];
  injection_beamemmitance_v[2] = injection_beamemmitance_a[1][0];
  injection_beamemmitance_v[3] = injection_beamemmitance_e[1][0];
  //-----------------------------------------------------

  //-----------------------------------------------------
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  /* Time Offset by Waveform */ 
  // read waveform function ... 
  particle_time=0.0;
  particleGun->SetParticleTime( particle_time ); // Global Time Setting ( Unit = nsec )

}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  EventNumberFile.close();
  OutputFileBeamEnergy.close();

  BeamEnergyFile.close();

  delete particleGun;
  delete ebeam;
  delete waveform;
}
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //----------------------------------
  // Primary Particle 
  //----------------------------------
  G4String PrimaryParticle;
  switch( beam_particle ){
  case 0:  
           PrimaryParticle = "gamma";
           G4cout << "beam_particle is gamma" << G4endl;
           break;
  case -1:
           PrimaryParticle = "e-";
           G4cout << "beam_particle is electron" << G4endl;
           break;
  case  1:
           PrimaryParticle = "e+"; 
           G4cout << "beam_particle is positron" << G4endl;
           break;
  case  50:
           PrimaryParticle = "opticalphoton"; 
           G4cout << "beam_particle is " << G4endl;
           break;
  case -100: 
           PrimaryParticle = "geantino";  // Test particle
           G4cout << "beam_particle is geantino" << G4endl;
           break;
  default:
           PrimaryParticle = "e-";
           G4cout << "beam_particle is default electron" << G4endl;
  }  
  //----------------------------------
  //:::: Energy Setting ::::::
  G4double      Ebeam_energy(0.);
  //:::: Position Setting ::::
  G4ThreeVector Ebeam_position;
  //:::: Momentum Setting ::::::
  G4ThreeVector Ebeam_momentum;

  int NEVENT = anEvent->GetEventID();
  if( NEVENT%10==0 ){ EventNumberFile << NEVENT << G4endl; }

  // Beam Energy  
  switch( beam_energy_mode ){
  case 0:
    Ebeam_energy = constant_energy; 
    break;
  case 1:
    Ebeam_energy = ebeam->gaussian( gaussian_mean, gaussian_sigma );   
    break;
  case 2:
    Ebeam_energy = ebeam->crystal_ball( crystalball_mean[0],
                                        crystalball_sigma[0],
                                        crystalball_alpha[0],
                                        crystalball_n[0],
                                        crystalball_xmin, crystalball_xmax );
    break;
  case 3:
    BeamEnergyFile >> Ebeam_energy; 
    //G4cout << "BeamEnergy " << Ebeam_energy << G4endl;  
    break;  
  default:
    Ebeam_energy = primary_energy[0];
  }

  // Beam Emmitance 
  // Beam Position & Momentum
  injection_position_delta_h = 0.;
  injection_position_delta_v = 0.;

  injection_position_shift_totalx = 0.0;
  injection_position_shift_totaly = 0.0;

  transverse_phix = 0.0;
  transverse_phiy = 0.0;
  injection_transverse_momentum_x = 0;
  injection_transverse_momentum_y = 0;
  injection_transverse_momentum_z = 0;

  switch( beam_emmitance_mode ){

  case 0: 
    injection_position_delta_h = ebeam->pulse( (-1)*injection_position_sigma[0], injection_position_sigma[0]);
    injection_position_delta_v = ebeam->pulse( (-1)*injection_position_sigma[0], injection_position_sigma[0]);  
    break;

  case 1:
    injection_position_delta_h = 0.;
    injection_position_delta_v = 0.;   

    if( injection_position_sigma_x[0] != 0. ){
      injection_position_shift_totalx = ebeam->gaussian( injection_position_shift_x[0], injection_position_sigma_x[0] );
    }else{
      injection_position_shift_totalx = injection_position_shift_x[0];
    }

    if( injection_position_sigma_y[0] != 0. ){
      injection_position_shift_totaly = ebeam->gaussian( injection_position_shift_y[0], injection_position_sigma_y[0] );
    }else{
      injection_position_shift_totaly = injection_position_shift_y[0];
    }

    transverse_phix = injection_transverse_phi0x[0]
                      *(injection_position_shift_totalx-injection_position_shift_x[0])/injection_position_sigma_x[0];
    transverse_phiy = injection_transverse_phi0y[0]
                      *(injection_position_shift_totaly-injection_position_shift_y[0])/injection_position_sigma_y[0];

    injection_transverse_momentum_x = sin( transverse_phix*RADI );
    injection_transverse_momentum_y = sin( transverse_phiy*RADI );
    //injection_transverse_momentum_z = cos( transverse_phix*RADI ) + cos( transverse_phiy*RADI );
    injection_transverse_momentum_z = cos( transverse_phiy*RADI );
    break;


  case 2:
    injection_position_r       = ebeam->gaussian( 0.0, injection_position_sigma[0] );       
    injection_position_theta   = ebeam->rn_radian();           
    injection_position_delta_h = injection_position_r*cos(injection_position_theta); 
    injection_position_delta_v = injection_position_r*sin(injection_position_theta); 
    break;


  case 3:
    injection_position_delta_h = ebeam->emittance_position(injection_beamemmitance_h, 0.0);
    injection_position_delta_v = ebeam->emittance_position(injection_beamemmitance_v, 0.0);
    injection_direction_y[0]   = ebeam->emittance_direction(injection_beamemmitance_h, injection_position_delta_h);
    injection_direction_z[0]   = ebeam->emittance_direction(injection_beamemmitance_v, injection_position_delta_v);    
    injection_position_delta_h*=1E+3;
    injection_position_delta_v*=1E+3;
    break;


  case 4: 
    injection_position_delta_h = 0.;
    injection_position_delta_v = 0.;
    injection_position_shift_totalx = ebeam->pulse( (-1)*injection_position_sigma_x[0], injection_position_sigma_x[0] );  
    injection_position_shift_totaly = ebeam->pulse( (-1)*injection_position_sigma_y[0], injection_position_sigma_y[0] );  
    injection_position_shift_totalx+=injection_position_shift_x[0];
    injection_position_shift_totaly+=injection_position_shift_y[0];
    break;

  case 101102:
    // 2010.11.02 Special Version 
    injection_position_delta_h = ebeam->pulse( (-1)*injection_position_sigma[0], injection_position_sigma[0]);
    injection_position_delta_v = ebeam->pulse( (-1)*injection_position_sigma[0], injection_position_sigma[0]);
    double tmp_injection_position[3];
    tmp_injection_position[0]=injection_position_x[0];
    tmp_injection_position[1]=injection_position_y[0]+injection_position_delta_h;
    tmp_injection_position[2]=injection_position_z[0]+injection_position_delta_v;
    ebeam->special_beam_direction_20101102( tmp_injection_position, 
                                            injection_direction_x[0],
                                            injection_direction_y[0],
                                            injection_direction_z[0] );
    break;

  case 110530:
    // 2011.05.30 Special Version 

    injection_position_delta_h = 0.;
    injection_position_delta_v = 0.;

    injection_position_r       = 17.0;    
    injection_position_theta   = ebeam->rn_radian();           

    injection_position_shift_totalx = injection_position_r*cos(injection_position_theta); 
    injection_position_shift_totaly = injection_position_r*sin(injection_position_theta); 
  
    injection_direction_x[0] = sin(0.0319384)*cos(injection_position_theta);
    injection_direction_y[0] = sin(0.0319384)*sin(injection_position_theta);
    injection_direction_z[0] = cos(0.0319384);

    break;

  case 130912: 
    injection_position_delta_h = 0.;
    injection_position_shift_totaly = ebeam->pulse( (-1)*injection_position_sigma_y[0], injection_position_sigma_y[0] );  
    injection_position_delta_v      = ebeam->pulse( -4250.0, 4250.0 );  
    break;

  default:    
    injection_position_delta_h = 0.;
    injection_position_delta_v = 0.;
  }

  Ebeam_position.setX( injection_position_x[0] + injection_position_shift_totalx );
  Ebeam_position.setY( injection_position_y[0] + injection_position_shift_totaly 
                                               + injection_position_delta_h );
  Ebeam_position.setZ( injection_position_z[0] + injection_position_delta_v );

  Ebeam_momentum.setX( injection_direction_x[0] + injection_transverse_momentum_x );
  Ebeam_momentum.setY( injection_direction_y[0] + injection_transverse_momentum_y );
  Ebeam_momentum.setZ( injection_direction_z[0] + injection_transverse_momentum_z );
  /*
  Ebeam_momentum.setX( injection_transverse_momentum_x );         
  Ebeam_momentum.setY( injection_transverse_momentum_y ); 
  Ebeam_momentum.setZ( injection_transverse_momentum_z );
  */
      
  // Waveform Time
  if( beam_waveform_flag == 1 ){
      particle_time=waveform->Output_time();
      particleGun->SetParticleTime(particle_time);
  }else{
      particle_time=0.0;
  }

  OutputFileBeamEnergy 
    << setprecision(10) << " " << Ebeam_energy 
    << setprecision(10) << " " << Ebeam_position.x() 
    << setprecision(10) << " " << Ebeam_position.y() 
    << setprecision(10) << " " << Ebeam_position.z() 
    << setprecision(10) << " " << Ebeam_momentum.x() 
    << setprecision(10) << " " << Ebeam_momentum.y() 
    << setprecision(10) << " " << Ebeam_momentum.z()
    << setprecision(10) << " " << particle_time
    << G4endl;

  particleGun->SetParticleEnergy(Ebeam_energy*MeV);
  particleGun->SetParticlePosition( Ebeam_position );
  particleGun->SetParticleMomentumDirection(Ebeam_momentum.unit());

  if( PrimaryParticle == "opticalphoton" ){
  //::::: Photon Polarization Setting :::::                                                                                                   
    G4ThreeVector polar = G4ThreeVector( 0, 0, 1.0 );
    particleGun->SetParticlePolarization(polar);
  }

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle(PrimaryParticle);
  particleGun->SetParticleDefinition(particle);

  particleGun->GeneratePrimaryVertex(anEvent);


}

