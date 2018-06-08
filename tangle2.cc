// *******************************************************************
// * License and Disclaimer                                          *
// *                                                                 *
// * This software is copyright of Geant4 Associates International   *
// * Ltd (hereafter 'G4AI'). It is provided under the terms and      *
// * conditions described in the file 'LICENSE' included in the      *
// * software system.                                                *
// * Neither the authors of this software system nor G4AI make any   *
// * representation or warranty, express or implied, regarding this  *
// * software system or assume any liability for its use.            *
// * Please see the file 'LICENSE' for full disclaimer and the       *
// * limitation of liability.                                        *
// *******************************************************************
// $Id$
// John Allison  22nd May 2017
// Gary Smith    27th Nov 2017

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif
#include "G4PhysListFactory.hh"
#include "Tangle2DetectorConstruction.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "Tangle2ActionInitialization.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

#include "Tangle2Data.hh"

int main(int argc,char** argv)
{
  
  //------------------------
  // Graphics
  G4bool  useGraphics = false;
  
  //--------------------------
  // Beam

  // A - axis - photons emitted along x-axis if true or in random directions if false; must be false if positrons is true
  Tangle2::fixedAxis = false;
  
  // B - beam particle - entagled photons produced by positron annihilation if true, unentangled back-to-back photons if false
  Tangle2::positrons = true;
  
  // polarisation direction; perpPol must be false if positrons is true
  Tangle2::perpPol   = false;  //determines whether polarisation is orthogonal or random for back-to-back photons
  Tangle2::polYZ     = false;  //for use with fixed beam in x direction; if true, one photon is y-polarised and the other z-polarised

  // C - polarised/unpolarised Compton scattering
  G4bool  usePolarisedCompton = true;

  //--------------------------
  // Detector 
  // True: arrays 90 cm apart (45 cm from source)  
  // False: arrays 6 cm apart (3 cm from source)  
  Tangle2::fullPET   = false;

  // Do this first to capture all output
  G4UIExecutive* ui = nullptr;
  
  if(useGraphics)
    ui = new G4UIExecutive(argc, argv);
  
  // Choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  
  // Set seed using system time
  G4long seeds[2];
  time_t systime = time(NULL);
  seeds[0] = (long) systime;
  seeds[1] = (long) (systime*G4UniformRand());
  G4Random::setTheSeeds(seeds);
  
#ifdef G4MULTITHREADED
  G4MTRunManager* runManager = new G4MTRunManager;
#else
  G4RunManager* runManager = new G4RunManager;
#endif
 
  runManager->SetUserInitialization(new Tangle2DetectorConstruction);
    
  G4int verbose;
  G4PhysListFactory factory;
  
  G4VModularPhysicsList* physList = factory.GetReferencePhysList("FTFP_BERT");
  physList->SetVerboseLevel(verbose = 1);
  
  if(usePolarisedCompton)
    physList->ReplacePhysics(new G4EmLivermorePolarizedPhysics);
  else
    physList->ReplacePhysics(new G4EmLivermorePhysics); 

  runManager->SetUserInitialization(physList);

  runManager->SetUserInitialization(new Tangle2ActionInitialization);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  
  if(useGraphics)
    UImanager->ApplyCommand("/control/execute visGraph.mac");
  else
    UImanager->ApplyCommand("/control/execute visNoGraph.mac");
  
  G4cout << " ------------------------------------------ " << G4endl;
  G4cout << " QETlab simulation " << G4endl;
  G4cout <<  G4endl;

  if(Tangle2::fullPET){
    G4cout << " Human PET diameter " << G4endl;
  }
  else{
    G4cout << " Lab experiment diameter " << G4endl;
  }
    
  G4cout <<  G4endl;
  G4cout << " Generated : " << G4endl;
  
  // Print beam choices to screen
  if(Tangle2::positrons)
    G4cout << " Positrons " << G4endl;
  else{
    G4cout << " Back to back gammas. " << G4endl;
    if(Tangle2::fixedAxis)
      G4cout << " On a fixed axis, " << G4endl;
    else
      G4cout << " Covering detector area, " << G4endl;
    if(Tangle2::perpPol)
      G4cout << " with perpendicular polarisation. " << G4endl;
    else
      G4cout << " with random relative polarisation. " << G4endl;
  }
  if(usePolarisedCompton)
    G4cout << " Polarized Compton scattering " << G4endl;
  else
    G4cout << " UnPolarized Compton scattering " << G4endl;
     
  G4cout << " ------------------------------------------ " << G4endl;
  
     if(useGraphics)
       ui->SessionStart();
  
  delete ui;
  delete visManager;
  delete runManager;
}
