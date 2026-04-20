//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb  4 11:57:18 2026 by ROOT version 6.36.04
// from TTree RootTuple/RootTuple
// found on file: 10MeV/EventOutput.Real.00000001.root
//////////////////////////////////////////////////////////

#ifndef readr_h
#define readr_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class readr {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        weight;
   vector<int>     *barcode;
   vector<double>  *Px;
   vector<double>  *Py;
   vector<double>  *Pz;
   vector<double>  *E;
   Double_t        b;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_barcode;   //!
   TBranch        *b_Px;   //!
   TBranch        *b_Py;   //!
   TBranch        *b_Pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_b;   //!

  readr(char *fname, int nfile, TTree *tree=0);
   virtual ~readr();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(char *fname);
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef readr_cxx
readr::readr(char *fname, int nfile,TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

     TChain *chain = new TChain("RootTuple","RootTuple");
     for (int ifile = 1; ifile <= nfile; ifile++){
       chain->Add(Form("./%s/EventOutput.Real.0000000%i.root",fname,ifile));
       printf("File Added %s",Form("./%s/EventOutput.Real.0000000%i.root\n",fname,ifile));
     }
     //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("10MeV/EventOutput.Real.00000001.root");
     //if (!f || !f->IsOpen()) {
     //    f = new TFile("10MeV/EventOutput.Real.00000001.root");
     // }
     // f->GetObject("RootTuple",tree);
     tree=chain;
   }
   Init(tree);
}

readr::~readr()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t readr::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t readr::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void readr::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   barcode = 0;
   Px = 0;
   Py = 0;
   Pz = 0;
   E = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("barcode", &barcode, &b_barcode);
   fChain->SetBranchAddress("Px", &Px, &b_Px);
   fChain->SetBranchAddress("Py", &Py, &b_Py);
   fChain->SetBranchAddress("Pz", &Pz, &b_Pz);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("b", &b, &b_b);
   Notify();
}

bool readr::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void readr::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t readr::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef readr_cxx
