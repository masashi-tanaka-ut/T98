//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4::PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"


#include "TROOT.h"     
#include "TH1.h"       
#include "TFile.h"     
#include "TMath.h"     
#include "TTree.h"     
#include "TVector3.h"  
#include "TRandom3.h"  
#include "TString.h"   
#include <vector>

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(int nevent)
{
  Nevent = nevent;
  
  G4int nofParticles = 1;
  fParticleGun = new G4ParticleGun(nofParticles);
  // default particle kinematic
  //
  auto particleDefinition
    = G4ParticleTable::GetParticleTable()->FindParticle("anti_proton");
  fParticleGun->SetParticleDefinition(particleDefinition);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));
  fParticleGun->SetParticleMomentum(700.*MeV);  

  inputRootFile = TFile::Open("../input/T98_00518.root"); 
  if (!inputRootFile || inputRootFile->IsZombie()) {  
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction()",  
		"FileError", FatalException, 
		"Cannot open ROOT file"); 
  }
  tree = (TTree*)inputRootFile->Get("tree"); 
  if (!tree) { 
    G4Exception("PrimaryGeneratorAction::PrimaryGeneratorAction()",  
		"TreeError", FatalException, 
		"Cannot find tree in ROOT file");  
  }

  tree->SetBranchAddress("x_bpc2", &x_bpc2); 
  tree->SetBranchAddress("y_bpc2", &y_bpc2); 
  tree->SetBranchAddress("dxdz_bpc2", &dxdz_bpc2); 
  tree->SetBranchAddress("dydz_bpc2", &dydz_bpc2); 
  tree->SetBranchAddress("trig_kaon3", &trig_kaon3);
  tree->SetBranchAddress("beam_tof", &beam_tof);
  tree->SetBranchAddress("beam_dp", &beam_dp);
  fNEntries = tree->GetEntries();

  // G4cout << "Nevent : " << Nevent << G4endl; 
  // G4cout << "Entry : " << fEntry << G4endl;
  fEntry = 0;
  int ievent=0;
  while(Nevent>ievent){
    tree->GetEntry(fEntry);
    //if(abs(beam_tof)<5.&&abs(x_bpc2)<100)ievent++; //pion
    if(abs(beam_tof+12)<7.&&abs(x_bpc2)<100&&trig_kaon3)ievent++; //pbar/proton
    //G4cout << "Skipped Pion Event : " << ievent << " Entry : " <<  fEntry << G4endl;
    fEntry++;  	
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  if (inputRootFile) { 
    inputRootFile->Close();  
    delete inputRootFile;  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // This function is called at the begining of event

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get world volume
  // from G4LogicalVolumeStore
  //

  TRandom3*tr3;
  
  auto analysisManager = G4AnalysisManager::Instance();
  G4double worldZHalfLength = 0.;
  auto worldLV = G4LogicalVolumeStore::GetInstance()->GetVolume("World");

  // Check that the world volume has box shape
  G4Box* worldBox = nullptr;
  if (  worldLV ) {
    worldBox = dynamic_cast<G4Box*>(worldLV->GetSolid());
  }

  if ( worldBox ) {
    worldZHalfLength = worldBox->GetZHalfLength();
  }
  else  {
    G4ExceptionDescription msg;
    msg << "World volume of box shape not found." << G4endl;
    msg << "Perhaps you have changed geometry." << G4endl;
    msg << "The gun will be place in the center.";
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()",
      "MyCode0002", JustWarning, msg);
  }

  //T98 Beam Profile 
  tr3=new TRandom3();
  tr3->SetSeed(0);

  if (fEntry >= fNEntries) {  
    G4Exception("PrimaryGeneratorAction::GeneratePrimaries()", 
		"EOF", JustWarning, 
		"Reached end of ROOT tree, restarting from first entry");  
    fEntry = 0;
  }
  
  G4cout << "Start Fentry : " << fEntry << G4endl; 
  
  //fEntry = Nevent;
  //Max Pbar Event:714771
  while(1){
    fEntry = (int)tr3->Uniform(714771);
    tree->GetEntry(fEntry);
    //if(abs(beam_tof)<5.&&abs(x_bpc2)<100)ievent++; //pion
    if(abs(beam_tof+12)<7.&&abs(x_bpc2)<100&&trig_kaon3)break; //pbar/proton
    //fEntry++;
  } 
  G4cout << "Fentry : " << fEntry << " X_PBC2 : " << x_bpc2 <<  G4endl;
  
  TVector3 vecPos(1.0, 0.0, 5.0);
  TVector3 vecMom(1.0, 0.0, 0.0);
  
  //vecPos[0] = 400.; 
  vecPos[0] = 2600.; 
  vecPos[1] = y_bpc2;  
  vecPos[2] = -1*x_bpc2;
  //vecPos[2] = -1*x_bpc2*10;
  vecMom[0] = -1.;
  vecMom[1] = dydz_bpc2;
  vecMom[2] = -1*dxdz_bpc2;
  // vecMom[0] = -1.;
  // vecMom[1] = 0;
  // vecMom[2] = 0;
  vecMom = vecMom.Unit();
  
  //double Momentum = 700;
  double Momentum = 700.*(1.+beam_dp/100.);
  ////Random Beam
  // TVector3 vecPosSigma(0, 3.0, 3.0);
  // TVector3 vecMomSigma(0, 0.002, 0.002);
  // vecPos[0] = 400.; 
  // vecPos[1] = 0.;
  // vecPos[2] = 0.;
  // G4double r_z,r_y;
  // r_y = tr3->Gaus(0, vecMomSigma[1]);
  // r_z = tr3->Gaus(0, vecMomSigma[2]); 
  // vecMom[0] = -1.;
  // vecMom[1] = r_y;
  // vecMom[2] = r_z;
  // vecMom = vecMom.Unit();
  
  
  G4cout << "X Pos : " << vecPos[0] << " Y Pos : " << vecPos[1] << " Z Pos : " << vecPos[2] << G4endl;
  G4cout << "X Mom : " << vecMom[0] << " Y Mom : " << vecMom[1] << " Z Mom : " << vecMom[2] << G4endl;
  G4cout << "X In : " << vecPos[0]+vecPos[0]*vecMom[0] << " Y In : " << vecPos[1]+vecPos[0]*vecMom[1] << " Z In : " << vecPos[2]+vecPos[0]*vecMom[2] << G4endl;
  
  // Set gun position
  fParticleGun->SetParticlePosition(G4ThreeVector(vecPos[0]*mm,vecPos[1]*mm,vecPos[2]*mm)); 
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(vecMom[0],vecMom[1],vecMom[2]));
  fParticleGun->SetParticleMomentum(Momentum*MeV);  
  
  analysisManager->FillNtupleIColumn(2,0, anEvent->GetEventID());
  analysisManager->FillNtupleDColumn(2,1, vecPos[0]);
  analysisManager->FillNtupleDColumn(2,2, vecPos[1]);
  analysisManager->FillNtupleDColumn(2,3, vecPos[2]);
  analysisManager->FillNtupleDColumn(2,4, dxdz_bpc2);
  analysisManager->FillNtupleDColumn(2,5, dydz_bpc2);
  analysisManager->FillNtupleDColumn(2,6, Momentum);
  analysisManager->AddNtupleRow(2);
  
  fParticleGun->GeneratePrimaryVertex(anEvent);
  fEntry++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
