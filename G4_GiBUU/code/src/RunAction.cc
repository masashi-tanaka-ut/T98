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
/// \file B4/B4a/src/RunAction.cc
/// \brief Implementation of the B4::RunAction class

#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

namespace B4
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);
\
  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();

  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  analysisManager->CreateH1("cap","Process", 11, -0.5, 10.5);
  
  // Creating ntuple
  //
  analysisManager->CreateNtuple("trk", "track");
  analysisManager->CreateNtupleIColumn("event");
  analysisManager->CreateNtupleIColumn("pid");
  analysisManager->CreateNtupleDColumn("mass");
  analysisManager->CreateNtupleIColumn("trkid");
  analysisManager->CreateNtupleIColumn("mat");
  analysisManager->CreateNtupleDColumn("t");
  analysisManager->CreateNtupleDColumn("x");
  analysisManager->CreateNtupleDColumn("y");
  analysisManager->CreateNtupleDColumn("z");
  analysisManager->CreateNtupleDColumn("e");
  analysisManager->CreateNtupleDColumn("l");
  // analysisManager->CreateNtupleDColumn("xdir");
  // analysisManager->CreateNtupleDColumn("ydir");
  // analysisManager->CreateNtupleDColumn("zdir");
  // analysisManager->CreateNtupleDColumn("k");
  // analysisManager->CreateNtupleIColumn("parid");
  // analysisManager->CreateNtupleIColumn("gparid");
  // analysisManager->CreateNtupleIColumn("g2parid");
  // analysisManager->CreateNtupleIColumn("parpid");
  // analysisManager->CreateNtupleIColumn("gparpid");
  // analysisManager->CreateNtupleIColumn("g2parpid");
  // analysisManager->CreateNtupleIColumn("killflag");
  
  analysisManager->CreateNtuple("Header", "Header");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("CaptureFlag");
  analysisManager->CreateNtupleDColumn("InjecPosx");
  analysisManager->CreateNtupleDColumn("InjecPosy");
  analysisManager->CreateNtupleDColumn("InjecPosz");
  analysisManager->CreateNtupleDColumn("StopPosx");
  analysisManager->CreateNtupleDColumn("StopPosy");
  analysisManager->CreateNtupleDColumn("StopPosz");
  analysisManager->CreateNtupleDColumn("InjecDirx");
  analysisManager->CreateNtupleDColumn("InjecDiry");
  analysisManager->CreateNtupleDColumn("InjecDirz");
  analysisManager->CreateNtupleDColumn("E_in");
  // analysisManager->CreateNtupleIColumn("NumPi0");
  // analysisManager->CreateNtupleIColumn("NumSecPi0");

  analysisManager->CreateNtuple("Beam", "Beam");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleDColumn("x_bpc2");
  analysisManager->CreateNtupleDColumn("y_bpc2");
  analysisManager->CreateNtupleDColumn("z_bpc2");
  analysisManager->CreateNtupleDColumn("dxdz_bpc2");
  analysisManager->CreateNtupleDColumn("dydz_bpc2");
  analysisManager->CreateNtupleDColumn("beam_Mom");

  analysisManager->CreateNtuple("Par", "Par");
  analysisManager->CreateNtupleIColumn("EventID");
  analysisManager->CreateNtupleIColumn("pProcess");
  analysisManager->CreateNtupleIColumn("pmat");
  analysisManager->CreateNtupleIColumn("ptrkid");
  analysisManager->CreateNtupleIColumn("ppid");
  analysisManager->CreateNtupleIColumn("pparid");
  analysisManager->CreateNtupleDColumn("pMom_x");
  analysisManager->CreateNtupleDColumn("pMom_y");
  analysisManager->CreateNtupleDColumn("pMom_z");
  analysisManager->CreateNtupleDColumn("pMom");
  analysisManager->CreateNtupleDColumn("pk");
  
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  G4String fileName = "root/Geant4/PBAR/B4.root";
  //G4String fileName = "root/Geant4/PION/B4.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
  // G4String fileName = "B4.xml";
  analysisManager->OpenFile();
  //analysisManager->OpenFile(fileName);
  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();
  // if ( analysisManager->GetH1(1) ) {
  //   G4cout << G4endl << " ----> print histograms statistic ";
  //   if(isMaster) {
  //     G4cout << "for the entire run " << G4endl << G4endl;
  //   }
  //   else {
  //     G4cout << "for the local thread " << G4endl << G4endl;
  //   }
   
  //   G4cout << " EAbs : mean = "
  //      << G4BestUnit(analysisManager->GetH1(0)->mean(), "Energy")
  //      << " rms = "
  //      << G4BestUnit(analysisManager->GetH1(0)->rms(),  "Energy") << G4endl;

  //   G4cout << " EGap : mean = "
  //      << G4BestUnit(analysisManager->GetH1(1)->mean(), "Energy")
  //      << " rms = "
  //      << G4BestUnit(analysisManager->GetH1(1)->rms(),  "Energy") << G4endl;

  //   G4cout << " LAbs : mean = "
  //     << G4BestUnit(analysisManager->GetH1(2)->mean(), "Length")
  //     << " rms = "
  //     << G4BestUnit(analysisManager->GetH1(2)->rms(),  "Length") << G4endl;

  //   G4cout << " LGap : mean = "
  //     << G4BestUnit(analysisManager->GetH1(3)->mean(), "Length")
  //     << " rms = "
  //     << G4BestUnit(analysisManager->GetH1(3)->rms(),  "Length") << G4endl;            
  // }

  // save histograms & ntuple
  //
  
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
