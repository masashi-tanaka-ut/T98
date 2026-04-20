#define readr_cxx
#include "readr.h"
#include "myheader.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void readr::Loop(char *fname)
{
//   In a ROOT session, you can do:
//      root> .L readr.C
//      root> readr t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   TH1D *hid = new TH1D("hid",Form("%s; Partile ID; Number of Particles per Events",fname),5000,-2500,2500);
   TH1D *hid2 = new TH1D("hid2","hid2; Partile ID; Number of Particles per Events",5000,-2500,2500);

   TH1D *hepbar = new TH1D("hepbar","hepbar; Kinetic Energy (GeV); Particles/Event/10 MeV",1000,0,1);
   TH1D *hepip = new TH1D("hepip","hepip; Kinetic Energy (GeV); Particles/Event/10 MeV",100,0,1);
   TH1D *hepim = new TH1D("hepim","hepim",100,0,1);
   TH1D *hepi0 = new TH1D("hepi0","hepi0",100,0,1);
   TH1D *hep = new TH1D("hep","hep",100,0,1);
   TH1D *hen = new TH1D("hen","hen",100,0,1);
  
   TH1D *hnpi  = new TH1D("hnpi","hnpi; Number of Particles (GeV); Fraction of Event",11,-0.5,10.5);
   TH1D *hnpip = new TH1D("hnpip","hnpip; Number of Particles (GeV); Fraction of Event",11,-0.5,10.5);
   TH1D *hnpim = new TH1D("hnpim","hnpim",11,-0.5,10.5);
   TH1D *hnpi0 = new TH1D("hnpi0","hnpi0",11,-0.5,10.5);
   TH1D *hnp = new TH1D("hnp","hnp",11,-0.5,10.5);
   TH1D *hnn = new TH1D("hnn","hnn",11,-0.5,10.5);
   
   
   TNtuple *n1 = new TNtuple("n1","n1","iev:ip:id:m:p:e:ke");

   int iev=0;
   int ieva=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      ieva++;
      int npar = barcode->size();
      int npip=0;
      int npim=0;
      int npi0=0;
      int np=0;
      int nn=0;
      int npbar = 0;
      for (int ip=0;ip<npar;ip++){
	if (barcode->at(ip) == -2212) {
	  npbar++;
	  hepbar->Fill(E->at(ip)-0.93827);
	}
      }
      if (npbar>0) continue;

      for (int ip=0;ip<npar;ip++){
	float par[7];
	par[0]=ientry;
	par[1]=ip;
	
	par[2]=barcode->at(ip)*1.0;
	par[4]=sqrt(Px->at(ip)*Px->at(ip)+Py->at(ip)*Py->at(ip)+Pz->at(ip)*Pz->at(ip));
	par[3]=sqrt(E->at(ip)*E->at(ip)-par[4]*par[4]);
	par[5]=E->at(ip);
	par[6]=par[5]-par[3];
	hid->Fill(par[2]);
	if (par[6]>0.05) hid2->Fill(par[2]);
	if (barcode->at(ip)==211) {
	  hepip->Fill(par[6]);
	  npip++;
	}
	if (barcode->at(ip)==-211) {
	  hepim->Fill(par[6]);
	  npim++;
	}
	if (barcode->at(ip)==111) {
	  hepi0->Fill(par[6]);
	  npi0++;
	}
	if (barcode->at(ip)==2212) {
	  hep->Fill(par[6]);
	  if (par[6]>0.05) np++;
	}
	if (barcode->at(ip)==2112) {
	  hen->Fill(par[6]);
	  if (par[6]>0.05) nn++;
	}
	//n1->Fill(par);
	
      }

      hnpi->Fill(npip+npim+npi0);
      hnpip->Fill(npip);
      hnpim->Fill(npim);
      hnpi0->Fill(npi0);
      hnp->Fill(np);
      hnn->Fill(nn);
      iev++;
   }

   printf("********** %i\n",iev);
   hid->SetTitle(Form("%s Pbar Ann. %f",fname, (iev+0.0)/(ieva+0.1)));
   hid->Scale(1./iev);
    hid2->Scale(1./iev);
    hepip->Scale(1./iev);
    hepim->Scale(1./iev);
    hepi0->Scale(1./iev);
    hep->Scale(1./iev);
    hen->Scale(1./iev);

    
    hnpi->Scale(1./iev);
    hnpip->Scale(1./iev);
    hnpim->Scale(1./iev);
    hnpi0->Scale(1./iev);
    hnp->Scale(1./iev);
    hnn->Scale(1./iev);

    hid2->SetLineColor(2);
    hepip->SetLineColor(2);
    hepim->SetLineColor(4);
    hepi0->SetLineColor(1);
    hep->SetLineColor(6);
    hen->SetLineColor(7);
    
    TCanvas *c1 = new TCanvas("c1","c1",1200,800);
    c1->Divide(2,2);
    c1->cd(1)->SetLogy();
    //hid->GetXaxis()->SetRangeUser(-500,2500);
    hid->Draw("histo");
    hid2->Draw("histosame");
    TLegend *l0 = new TLegend(0.2,0.8,0.4,0.9);
    l0->AddEntry(hid,"All");
    l0->AddEntry(hid2,"KE > 50 MeV");
    l0->Draw();
    
    //c1->cd(2)->SetLogy();
    //hepip->SetMaximum(10);
    //hepip->SetMinimum(1e-4);
    c1->cd(2)->SetLogy(0);
    hepip->SetMaximum(0.1);
    //hepip->SetMinimum(1e-4);
    hepip->Draw("hist");
    hepim->Draw("histosame");
    hepi0->Draw("histosame");
    hep->Draw("histosame");
    hen->Draw("histosame");
    TLegend *l1 = new TLegend(0.5,0.6,0.8,0.9);
    l1->AddEntry(hep,"Proton");
    l1->AddEntry(hen,"Neutron");
    l1->AddEntry(hepip,"#pi^{+}");
    l1->AddEntry(hepim,"#pi^{-}");
    l1->AddEntry(hepi0,"#pi^{0}");
    l1->Draw();

    hnpi->SetMaximum(.5);
    hnpip->SetMaximum(.5);
    
    hnpip->SetLineColor(2);
    hnpim->SetLineColor(4);
    hnpi0->SetLineColor(1);
    hnp->SetLineColor(6);
    hnn->SetLineColor(7);
    
    c1->cd(3);
    hnpip->Draw("hist");
    hnpim->Draw("histsame");
    hnpi0->Draw("histsame");
    TLegend *l2 = new TLegend(0.5,0.6,0.8,0.9);
    l2->AddEntry(hnpip,"#pi^{+}");
    l2->AddEntry(hnpim,"#pi^{-}");
    l2->AddEntry(hnpi0,"#pi^{0}");
    l2->Draw();

    c1->cd(4);
    hnpi->Draw("hist");
    hnp->Draw("histsame");
    hnn->Draw("histsame");
    TLegend *l3 = new TLegend(0.6,0.6,0.8,0.9);
    l3->AddEntry(hnp,"Proton");
    l3->AddEntry(hnn,"Neutron");
    l3->AddEntry(hnpi,"#pi");
    l3->Draw();
   
    c1->Print(Form("%s.png",fname));


}
