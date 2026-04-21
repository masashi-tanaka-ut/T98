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
/// \file B4/B4a/src/SteppingAction.cc
/// \brief Implementation of the B4a::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4ios.hh"
#include "G4DynamicParticle.hh"

#include <vector>
#include <iostream>
#include <map>
#include <fstream>

using namespace B4;
G4int PreTID, lastEventID, num_pi0=0, num_secpi0=0;G4String Premat;
G4int CapFlag = 0;
bool injec = false;
bool stop = false;
G4double StopPos[3]={300,300,300} ,InjecPos[3]={300,300,300}, InjecDir[3]={10,10,10}, E_in=0;
std::map<int,std::pair<int,int>> ParentID;
std::map<int,std::pair<int,int>> GrandParentID;

namespace B4a
{

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(const DetectorConstruction* detConstruction,
                               EventAction* eventAction)
  : fDetConstruction(detConstruction),
    fEventAction(eventAction)
{
  
  //gfile = TFile::Open("EventOutput.Real.00000001.root");
  for (int i=0;i<10;i++){
  
    gfile[i] = TFile::Open(Form("../input/GiBUU/%iMeV/EventOutput.Real.00000001.root",(i+1)*25));
    gtree[i] = (TTree*)gfile[i]->Get("RootTuple");
    printf("GiBUU root file open ../input/GiBUU/%iMeV/EventOutput.Real.00000001.root\n",(i+1)*25);
    gtree[i]->SetBranchAddress("weight", &weight);
    gtree[i]->SetBranchAddress("barcode", &barcode);
    gtree[i]->SetBranchAddress("Px", &Px);
    gtree[i]->SetBranchAddress("Py", &Py);
    gtree[i]->SetBranchAddress("Pz", &Pz);
    gtree[i]->SetBranchAddress("E", &E);
    gtree[i]->SetBranchAddress("b", &b);
    ncount[i]=0;
  }
  tr = new TRandom();
  tr->SetSeed(1);
  TFile *hfile = TFile::Open("../input/GiBUU/prob.root");
  TH1D *heprob = (TH1D*) hfile->Get("heprob");
  for (int i=0;i<250;i++){
    eprob[i]=heprob->GetBinContent(i);
  }
}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  
void SteppingAction::UserSteppingAction(const G4Step* step)
{
  // Collect energy and track length step by step
  
  auto analysisManager = G4AnalysisManager::Instance();
  const G4Event* event =
    G4RunManager::GetRunManager()->GetCurrentEvent();
  
  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();
  
  // step length
  G4double stepLength = 0.;
  int pid = step->GetTrack()->GetDefinition()->GetPDGEncoding();
  int eventID = event->GetEventID();
  auto track = step->GetTrack();
  int trackid = track->GetTrackID();
  int parentid = track->GetParentID();
  G4String MatName  = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  G4String PName  =   track->GetDefinition()->GetParticleName();
  G4String Prname =   step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();

  int mat = 0;
  double kene = step->GetTrack()->GetKineticEnergy();
  if (pid == -2212) {
    int ikene = (int)kene;
    double ppp = eprob[ikene]*edep;
    
    //printf("*****Pbar****%i %f MeV iE %i dE %f Prob %f P %f\n",eventID,kene,ikene,edep,eprob[ikene],ppp);
    if (tr->Rndm()<ppp || kene<25.) {

      
      printf("*****Pbar1***%i %f MeV iE %i dE %f Prob %f P %f\n",eventID,kene,ikene,edep,eprob[ikene],ppp);	

      int itree=kene/25;

      // Exclude GiBUU events with final state pbar
      int npbar=1;      
      while(npbar!=0){
	npbar=0;
	gtree[itree]->GetEntry(ncount[itree]);
	int npar = barcode->size();
	for (int ip=0;ip<npar;ip++){
	  if (barcode->at(ip)==-2212) {
	    npbar++;
	  }
	}
	ncount[itree]++;
      }

      printf("*****Pbar2***%i %i %i %i\n",eventID,itree,ncount[itree],npbar);	

      // Rotation Matrix GiBUU (0,0,1) -> pbar track direction
      G4ThreeVector z(0,0,1);
      G4ThreeVector dir = track->GetMomentumDirection();
      G4ThreeVector axis = z.cross(dir);
      double angle = z.angle(dir);
      G4RotationMatrix rot;
      rot.rotate(angle, axis); 
      
      int npar = barcode->size();
      for (int ip=0;ip<npar;ip++){
	
    	G4ThreeVector gmom(Px->at(ip), Py->at(ip), Pz->at(ip)); 
	G4ThreeVector mom = rot * gmom; // Rotate
	
    	double p = mom.mag();
    	double ene = E->at(ip);
    	double mass = sqrt((ene+p)*(ene-p));
    	double ke  = ene - mass;

    	//printf("GenGen: %i %i %i mass:%f p:%f E:%f KE:%f\n",ip,npar2,barcode->at(ip),mass,p,ene,ke);
    	auto particleDef =G4ParticleTable::GetParticleTable()->FindParticle(barcode->at(ip));
        G4DynamicParticle* dynamicParticle = new G4DynamicParticle(particleDef, mom, ke);
        G4StepPoint* postStepPoint = step->GetPostStepPoint();
        G4Track* newTrack = new G4Track(dynamicParticle, postStepPoint->GetGlobalTime(), postStepPoint->GetPosition());
        
        newTrack->SetParentID(track->GetTrackID());
    	G4Step* nonConstStep = const_cast<G4Step*>(step);
    	nonConstStep->GetfSecondary()->push_back(newTrack);
      }
      step->GetTrack()->SetTrackStatus(fStopAndKill); // Delete pbar
    }
  }  

  
  if ( (step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. && edep>1.0e-6) || (pid==111 || pid==221 || pid==331 || pid==22)) {
    stepLength = step->GetStepLength();
    if (MatName == "Detector")  mat = 1;
    if (MatName == "LAr")       mat = 2;
    if (MatName == "Vessel")    mat = 3;
    analysisManager->FillNtupleIColumn(0,0, fEventAction->Geteventid());
    analysisManager->FillNtupleIColumn(0,1, step->GetTrack()->GetDefinition()->GetPDGEncoding());
    analysisManager->FillNtupleDColumn(0,2, (double) step->GetTrack()->GetDefinition()->GetPDGMass());
    analysisManager->FillNtupleIColumn(0,3, trackid);//Track ID
    analysisManager->FillNtupleIColumn(0,4, mat);
    analysisManager->FillNtupleDColumn(0,5, step->GetPreStepPoint()->GetGlobalTime());
    
    analysisManager->FillNtupleDColumn(0,6, step->GetPreStepPoint()->GetPosition().x());
    analysisManager->FillNtupleDColumn(0,7, step->GetPreStepPoint()->GetPosition().y());
    analysisManager->FillNtupleDColumn(0,8, step->GetPreStepPoint()->GetPosition().z());
    analysisManager->FillNtupleDColumn(0,9, edep);
    analysisManager->FillNtupleDColumn(0,10, stepLength);
    // analysisManager->FillNtupleDColumn(0,10, step->GetTrack()->GetMomentumDirection().x());//Direction
    // analysisManager->FillNtupleDColumn(0,11, step->GetTrack()->GetMomentumDirection().y());
    // analysisManager->FillNtupleDColumn(0,12, step->GetTrack()->GetMomentumDirection().z());
    // analysisManager->FillNtupleDColumn(0,13, step->GetTrack()->GetKineticEnergy());//Kinetic Energy
    analysisManager->AddNtupleRow(0);
  }
  //Particle Information
  int cap=0;
  if (Prname == "hFritiofCaptureAtRest") cap=1;
  if (Prname == "anti_protonInelastic") cap=2;
  if (cap > 0) {
    G4cout << "Stopped** " << PName << " " << Prname << G4endl;    
    const std::vector<const G4Track*>* secondaries =
        step->GetSecondaryInCurrentStep();
    
    for (const auto* secTrack : *secondaries) {
      G4double ke  = secTrack->GetKineticEnergy();      // KE      
      if (ke > 1.0 * CLHEP::MeV) {
	analysisManager->FillNtupleIColumn(3,0, eventID);
	analysisManager->FillNtupleIColumn(3,1, cap);
	analysisManager->FillNtupleIColumn(3,2, mat);
	analysisManager->FillNtupleIColumn(3,3, secTrack->GetTrackID());
	analysisManager->FillNtupleIColumn(3,4, secTrack->GetDefinition()->GetPDGEncoding());
	analysisManager->FillNtupleIColumn(3,5, secTrack->GetParentID());
	analysisManager->FillNtupleDColumn(3,6, secTrack->GetMomentumDirection().x());
	analysisManager->FillNtupleDColumn(3,7, secTrack->GetMomentumDirection().y());
	analysisManager->FillNtupleDColumn(3,8, secTrack->GetMomentumDirection().z());
	analysisManager->FillNtupleDColumn(3,9, secTrack->GetMomentum().mag());
	analysisManager->FillNtupleDColumn(3,10, ke);
	analysisManager->AddNtupleRow(3);
      }
    }
  }
  
  //Header Information
  if(Premat!="Detector"&&MatName == "Detector"&&trackid==1){
    InjecPos[0] =  step->GetPreStepPoint()->GetPosition().x();//Injec
    InjecPos[1] =  step->GetPreStepPoint()->GetPosition().y();
    InjecPos[2] =  step->GetPreStepPoint()->GetPosition().z();
    
    InjecDir[0] =  step->GetTrack()->GetMomentumDirection().x();//InjecDire
    InjecDir[1] =  step->GetTrack()->GetMomentumDirection().y();
    InjecDir[2] =  step->GetTrack()->GetMomentumDirection().z();
    E_in = step->GetTrack()->GetKineticEnergy();
  }else if(Premat!="Detector"&&step->GetTrack()->GetTrackStatus() == fStopAndKill&&trackid==1){
    InjecPos[0] =  -300;//Injec
    InjecPos[1] =  -300;
    InjecPos[2] =  -300;

    InjecPos[0] =  10;//Injec
    InjecPos[1] =  10;
    InjecPos[2] =  10;
    E_in = 0;
  }
  if(trackid==1&&step->GetTrack()->GetTrackStatus() == fStopAndKill){
    StopPos[0] =  step->GetPreStepPoint()->GetPosition().x();//Injec
    StopPos[1] =  step->GetPreStepPoint()->GetPosition().y();
    StopPos[2] =  step->GetPreStepPoint()->GetPosition().z();
  }
  if(trackid==1&&Prname=="hFritiofCaptureAtRest"&&MatName == "Detector"){
    CapFlag=1;
    G4cout << "Captured " << G4endl;
  }else if(trackid==1&&Prname=="anti_protonInelastic"&&MatName == "Detector"){
    CapFlag=2;
    G4cout << "Inelastic " << G4endl;
  }else if(trackid==1&&Prname=="hFritiofCaptureAtRest"&&MatName != "Detector"){
    CapFlag=4;
    G4cout << "Capture Before Detector " << G4endl;
  }else if(trackid==1&&step->GetTrack()->GetTrackStatus() == fStopAndKill&&MatName != "Detector"){
    CapFlag=3;
    G4cout << "Out of Detector " << G4endl;
  }
  if(trackid==1&&step->GetTrack()->GetTrackStatus() == fStopAndKill){
    analysisManager->FillNtupleIColumn(1,0,eventID);
    analysisManager->FillNtupleIColumn(1,1,CapFlag);
    analysisManager->FillNtupleDColumn(1,2,InjecPos[0]);
    analysisManager->FillNtupleDColumn(1,3,InjecPos[1]);
    analysisManager->FillNtupleDColumn(1,4,InjecPos[2]);
    analysisManager->FillNtupleDColumn(1,5,StopPos[0]);
    analysisManager->FillNtupleDColumn(1,6,StopPos[1]);
    analysisManager->FillNtupleDColumn(1,7,StopPos[2]);
    analysisManager->FillNtupleDColumn(1,8,InjecDir[0]);
    analysisManager->FillNtupleDColumn(1,9,InjecDir[1]);
    analysisManager->FillNtupleDColumn(1,10,InjecDir[2]);
    analysisManager->FillNtupleDColumn(1,11,E_in);

    analysisManager->AddNtupleRow(1);
    CapFlag=0;
    StopPos[0]=300;StopPos[1]=300;StopPos[2]=300;
    InjecPos[0]=300;InjecPos[1]=300;InjecPos[2]=300;
    InjecDir[0]=10;InjecDir[1]=10;InjecDir[2]=10;
    E_in = 0;
  }
 
  // if(trackid==1&&step->GetTrack()->GetTrackStatus() == fStopAndKill){
  //   ParentID.clear();  GrandParentID.clear();
  // }
  // if(PreTID!=trackid){
  //   ParentID[trackid] = std::make_pair(parentid,pid);
  //   GrandParentID[trackid] = std::make_pair(ParentID[parentid].first, ParentID[parentid].second);
  // }

  //bool Print = true;
  bool Print = false;
  if(Print){
    //comment out check
    if(step->GetTrack()->GetKineticEnergy()>10&&parentid==1&&PreTID!=trackid&&CapFlag&&MatName =="Detector"){
      G4cout << "\033[31m"
	     << "PartileName: "<< step->GetTrack()->GetDefinition()->GetParticleName()
	     << "  TrackID: "  << step->GetTrack()->GetTrackID()
	     << "  ParentID: " << step->GetTrack()->GetParentID()
	     << "  Energy: "   << step->GetTrack()->GetKineticEnergy() << " MeV"
	     << " Length: " << step->GetStepLength()
	     << "\033[0m "
	     <<  G4endl;
    }
    //grand particle in Detector
    if(step->GetTrack()->GetKineticEnergy()>10&&ParentID[parentid].first==1&&step->GetTrack()->GetTrackStatus() == fStopAndKill){
      //if(step->GetTrack()->GetKineticEnergy()>10&&PreTID!=trackid&&ParentID[parentid].first==1&&MatName == "Detector"){
      G4cout << "\033[32m"
	     << "PartileName: "<< step->GetTrack()->GetDefinition()->GetParticleName()
	     << "  TrackID: "  << step->GetTrack()->GetTrackID()
	     << "  ParentID: " << step->GetTrack()->GetParentID()
	     << "  GrandParentID: " << ParentID[parentid].first
	     << "  Energy: "   << step->GetTrack()->GetKineticEnergy() << " MeV"
	     << " Pos: " << step->GetPreStepPoint()->GetPosition()
	     << " Length: " << step->GetStepLength()
	     << "\033[0m "
	     <<  G4endl;
    }

    //Out of Detector
    if(step->GetTrack()->GetKineticEnergy()>10&&PreTID!=trackid&&ParentID[parentid].first==1&&MatName != "Detector"){
      G4cout << "\033[34m"
	     << "PartileName: "<< step->GetTrack()->GetDefinition()->GetParticleName()
	     << "  TrackID: "  << step->GetTrack()->GetTrackID()
	     << "  ParentID: " << step->GetTrack()->GetParentID()
	     << "  GrandParentID: " << ParentID[parentid].first
	     << "  Energy: "   << step->GetTrack()->GetKineticEnergy() << " MeV"
	     << "  Edepo: "      << edep
	     << " Length: " << step->GetStepLength()
	     << "\033[0m "
	     <<  G4endl;
    }
    
    //electron from gamma 
    //if(step->GetTrack()->GetKineticEnergy()>10&&PreTID!=trackid&&GrandParentID[parentid].first==1&&MatName != "Detector"){
    if(step->GetTrack()->GetKineticEnergy()>0.001&&PreTID!=trackid&&GrandParentID[parentid].first==1&&MatName != "Detector"){
      G4cout << "\033[36m"
	     << "PartileName: "<< step->GetTrack()->GetDefinition()->GetParticleName()
	     << "  TrackID: "  << step->GetTrack()->GetTrackID()
	     << "  ParentID: " << step->GetTrack()->GetParentID()
	     << "  GrandParentID: " << ParentID[parentid].first
	     << "  GrandParentID2: "<< GrandParentID[parentid].first
	     << "  Energy: "   << step->GetTrack()->GetKineticEnergy() << " MeV"
	     << "  Edepo: "      << edep
	     << " Length: " << step->GetStepLength()
	     << "\033[0m "
	     <<  G4endl;
    }
    
    //Pi+- and gamma check
    if(parentid!=1&&step->GetTrack()->GetKineticEnergy()>10&&PreTID!=trackid&&MatName == "Detector"&&(PName=="gamma"||PName=="pi+"||PName=="pi-"||PName=="pi0"||PName=="e+"||PName=="e-")){
      G4cout << "\033[37m"
	     << "PartileName: "<< step->GetTrack()->GetDefinition()->GetParticleName()
	     << "  TrackID: "  << step->GetTrack()->GetTrackID()
	     << "  ParentID: " << step->GetTrack()->GetParentID()
	     << "  GrandParentID: " << ParentID[parentid].first
	     << "  Energy: "   << step->GetTrack()->GetKineticEnergy() << " MeV"
	     << "  Edepo: "      << edep
	     << " Length: " << step->GetStepLength()
	     << "\033[0m "
	     <<  G4endl;
    }
  }
  PreTID = step->GetTrack()->GetTrackID();
  Premat = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  lastEventID =  event->GetEventID();
}
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
