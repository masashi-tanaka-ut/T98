#include <iostream>
#include <fstream>
#include <vector>
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTree.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TLegend.h"  
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraph2D.h" 
#include "TVirtualFFT.h"     
#include "TLine.h"    
#include <filesystem>
//ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit");
void LoadRootFile(Int_t RunNumber, Int_t FileNumber);
void DefineHist();
void LoopHist();
void MakeWaveform(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber);
void MakeSimWaveform(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber);
void MakeWaveform2(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber);
void MakeBox();
void Tracking(Bool_t Rawtrack, Bool_t Rebintrack, Bool_t NoiseReductrack);
void Clustering();
void DrawHough(Int_t EventNumber);
void ChargeCalc();
void Hough();
void Ped();
void DrawPlot(Bool_t a, Bool_t b);
void Signal(Int_t RunNumber, Int_t FileNumber, Int_t EventNumber, Int_t nRebin, Bool_t DrawFlag);
Int_t Cut();
void Loop(Int_t RunNumber, Int_t StartFileNumber, Int_t EndFileNumber);


//For Merge analysis code
void MakeLTARSWaveform(Int_t RunNumber, Int_t FileNumber, Int_t EventNumber);
void MakeSISWaveform(Int_t RunNumber, Int_t FileNumber, Int_t EventNumber);
void DrawMergePlot(bool RawDisp,bool TrackDisp,bool MargePlot);

//For Peaking Search
void Peaking(Int_t EventNumber,Bool_t ShowPlot);
void Charge(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber);

//NR Waveform
void MakeSingleRootFile(Int_t RunNumber, Int_t FileNumber);
void LoadNRRootFile(Int_t RunNumber, Int_t FileNumber); 
void Tracking();

//BeamLine
void DrawBLPlot();
void DrawMaxPlot();
void DrawStopPlot();
void DrawChargePlot();
void Drawdedx();
void LoadBLRootFile(Int_t RunNumber);
void LoadChargeRootFile(Int_t RunNumber);
void LoadMaxRootFile(Int_t RunNumber,Int_t StartFileNumber, Int_t EndFileNumber);
void LoadStopRootFile(Int_t RunNumber,Int_t StartFileNumber, Int_t EndFileNumber);

//Loop
void LoopMaxpos(Int_t BLRunNumber, Int_t RunNumber, Int_t StartFileNumber, Int_t EndFileNumber, Bool_t DrawFlag);
void Loop(Int_t BLRunNumber, Int_t RunNumber, Int_t StartFileNumber, Int_t EndFileNumber, Bool_t DrawFlag);
void MakeCut();

void dEdx(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber, Bool_t DrawFlag);


//ana_Pion
void DrawPion();
void DrawCalib();
void DrawComp();
void AnaPion(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber, Bool_t DrawFlag);
void Plot(Int_t FileNumber,Int_t Sigma,Int_t Amp,Int_t EventNumber, Int_t Dx=0, Int_t Dy=0);
void AnaPion(Int_t EventNumber, Bool_t DrawFlag);
void LoopHist();
void LoopCalib();
void Landau();


//MC analysis  
struct track_info{
  double p0;
  double p1;
  double energy;  
  double tstart;  
  double tend; 
  double xstart;  
  double xend; 
  bool Chbool; 
  bool xChbool;
  bool zChbool;
  int type; // 0: primary, 1: secondary  
}; 
void HoughPrimary();
void HoughSecondary();
void SecondaryAnalysis(int plane_id, int particle_id);
void AngleSearch();
