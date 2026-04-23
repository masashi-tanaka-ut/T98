#include "myheader.h"


//const about Plot
const Int_t nch=256;
const Int_t nhch=128;
const Int_t nbin=1000;
const Int_t ptrig=125;
const Double_t timeconv = 0.4;
const Double_t ADCconv = (Double_t)2./pow(2,12);

const Double_t YMAX=500., YMIN=-500.;
Double_t TEND=350., ZMAX=300.;
const Double_t ymin1=-300., ymax1=300.;
const Double_t ymin2=-300.,ymax2=300.;

//const about analysis
const Int_t bgch=16;
const Int_t PED_STARTBIN=0, PED_STOPBIN=100;
Double_t HGRMSbin=120.,HGRMSbin_nr=50., LGRMSbin=5., LGRMSbin_nr=5.;


TString inFileName;
TFile *fin;
TTree *tree;
TBranch *b_trigCounter;
TBranch *b_trvecADC;
TBranch *b_CLKCounter;
unsigned int trigCounter[4];
unsigned int CLKCounter[4];
vector<vector<double>>*trvecADC,*Position,*Direction, *GammaStopPos;
vector<double> *Kinetic_Energy, *Length;
vector<int> *PID, *trkID,*ParID,*ParPID,*GParID,*GParPID,*G2ParID,*G2ParPID;
Int_t Header[20],MHeaderI[20], Process;Double_t HeaderD[20],MHeaderD[20],StopPos[3], InjecPos[3];


TLatex*latexhre;
Int_t fitbin=1;
//Global variables
Int_t runnumber;
Int_t filenumber;
Int_t eventnumber;

//TH1D 
TH1D *hped;
TH1D *hPedN;
TH1D *hWaveform[nch];
TH1D *hWaveform_raw[nch];
TH1D *hWaveform_nr[nch];

//TH2D
TH2D *hxz;
TH2D *hyz;
TH2D *hxz_nr;
TH2D *hyz_nr;

//Canvas
TCanvas *cTrack={NULL};
TCanvas *cWaveform[8]={NULL};

void LoadRootFile(char *fname){

  inFileName = Form("%s",fname);

  fin = TFile::Open(inFileName);
  delete tree;
  tree = (TTree*)fin->Get("NRWaveTree");
  tree->SetBranchAddress("trigCounter", &trigCounter, &b_trigCounter);
  //tree->SetBranchAddress("CLKCounter", &CLKCounter, &b_CLKCounter);
  tree->SetBranchAddress("trvecADC", &trvecADC, &b_trvecADC);
  tree->SetBranchAddress("Header", Header);
  tree->SetBranchAddress("HeaderD", HeaderD);
  tree->SetBranchAddress("MHeaderI", MHeaderI);
  tree->SetBranchAddress("MHeaderD", MHeaderD);
  // tree->SetBranchAddress("Position", &Position);
  // tree->SetBranchAddress("Direction", &Direction);
  // tree->SetBranchAddress("Kinetic_Energy", &Kinetic_Energy);
  // tree->SetBranchAddress("Length", &Length);
  // tree->SetBranchAddress("trkID", &trkID);
  // tree->SetBranchAddress("PID", &PID);
  // tree->SetBranchAddress("ParID", &ParID);
  // tree->SetBranchAddress("ParPID", &ParPID);
  // tree->SetBranchAddress("GParID", &GParID);
  // tree->SetBranchAddress("GParPID", &GParPID);
  // tree->SetBranchAddress("G2ParID", &G2ParID);
  // tree->SetBranchAddress("G2ParPID", &G2ParPID);
  // tree->SetBranchAddress("InjecPos", &InjecPos);
  // tree->SetBranchAddress("StopPos", &StopPos);
  // tree->SetBranchAddress("GammaPosition", &GammaStopPos);
  // tree->SetBranchAddress("Process", &Process);

  cout << "Load " << inFileName << endl;
}

double LiqVel(double EField)
{
    double vd = 0.0;
    double e = EField;
    double t = 89.0;
    double t0_ = 90.371;

    double p1 = -0.01481;
    double p2 = -0.0075;
    double p3 =  0.141;
    double p4 = 12.4;
    double p5 = 1.627;
    double p6 = 0.317;

    double po[6] = { -0.03229, 6.231, -10.62, 12.74, -9.112, 2.83 };

    if (e > 0.5) {
        vd = (p1 * (t - t0_) + 1.0)
           * (p3 * e * std::log(1.0 + p4 / e) + p5 * std::pow(e, p6))
           + p2 * (t - t0_);
    } else {
        for (int i = 0; i < 6; i++) {
            vd += po[i] * std::pow(e, i);
        }
        double etmp = 0.5;
        double tmp1 = 0.0;
        for (int i = 0; i < 6; i++) {
            tmp1 += po[i] * std::pow(etmp, i);
        }
        double tmp2 = (p1 * (t - t0_) + 1.0)
            * (p3 * etmp * std::log(1.0 + p4 / etmp) + p5 * std::pow(etmp, p6))
            + p2 * (t - t0_);
        vd *= (tmp2 / tmp1);
    }
    return vd;
}

void Plot(Int_t EventNumber){
  eventnumber = EventNumber;
  MakeWaveform(runnumber, filenumber, eventnumber);
  Tracking(true,false,false);
  //MakeBox();
  //Clustering();
  DrawPlot(false,true);

  double liqVel04 = LiqVel(0.25);
  cout <<" Drift Vel: " << liqVel04 << endl;
  TGraph*gstopx = new TGraph();  TGraph*gstopy = new TGraph();
  // gstopx -> SetPoint(0,150-MHeaderD[6],(150+MHeaderD[7])*liqVel04);  gstopx -> SetPoint(1,150-MHeaderD[0],(150+MHeaderD[2])*liqVel04);
  // gstopy -> SetPoint(0,150-MHeaderD[8],(150+MHeaderD[7])*liqVel04);  gstopy -> SetPoint(0,150-MHeaderD[1],(150+MHeaderD[2])*liqVel04);
  // cTrack->cd(1);   
  // gstopx -> SetMarkerColor(2); gstopx -> SetMarkerStyle(20);gstopx -> SetMarkerSize(2);gstopx -> Draw("same p");  
  // cTrack->cd(2); 
  // gstopy -> SetMarkerColor(2); gstopy -> SetMarkerStyle(20);gstopy -> SetMarkerSize(2);gstopy -> Draw("same p"); 
  cout << " X in : " << MHeaderD[0]
       << " Y: " << MHeaderD[1]
       << " Z: " << MHeaderD[2] << endl;
  cout << " X Stop : " << MHeaderD[6]
       << " Y: " << MHeaderD[7]
       << " Z: " << MHeaderD[8] << endl;
  // cout << " Particles: " << setw(3) << Header[0]
  //      << " Eta: " << setw(3) << Header[1]
  //      << " Deuteron: " << setw(3) << Header[2]
  //      << " Proton: " << setw(3) << Header[3]
  //      << " Pi+ or Pi-: " << setw(3) << Header[4]
  //      << " Pi0: " << setw(3) << Header[5]
  //      << " Gamma: " << setw(3) << Header[6]
  //      << " e+-: " << setw(3) << Header[7] << endl;
  // cout << " Secondary: " << setw(3) << Header[8]
  //      << " Sec Eta: " << setw(3) << Header[9]
  //      << " Sec Deuteron: " << setw(3) << Header[10]
  //      << " Sec Proton: " << setw(3) << Header[11]
  //      << " Sec Pi+ or Pi-: " << setw(3) << Header[12]
  //      << " Sec Pi0: " << setw(3) << Header[13]
  //      << " Sec Gamma: " << setw(3) << Header[14]
  //      << " Sec e+-: " << setw(3) << Header[15]<< endl;
  // for(int ii=0; ii < Header[0] ; ii++){
  //   cout << " PID     : " << PID->at(ii) << " ParID   : " << ParID->at(ii) << " ParPID   : " << ParPID->at(ii) << " trkID   : " << trkID->at(ii)<< " Length : " << Length->at(ii) << " Pos X : " << -1*Position->at(ii).at(0)+150. << " Y : " << -1*Position->at(ii).at(2)+150. << " Z : " << Position->at(ii).at(1)+150. << " Dir X : " << Direction->at(ii).at(0) << " Y : " << Direction->at(ii).at(2) << " Z : " << Direction->at(ii).at(1)  << endl;
  // }
  // cout << "Stop X : "<< 150-StopPos[0] << " Y : "<< 150-StopPos[2] << " Z : "<< 150+StopPos[1] << endl;
  // cout << "Injec X : "<< 150-InjecPos[0] << " Y : "<< 150-InjecPos[2] << " Z : "<< 150+InjecPos[1] << endl;
}

void Next(){
  Plot(eventnumber+1);
}  

void Previous(){
  Plot(eventnumber-1);
}

void MakeWaveform(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber){
  runnumber = RunNumber;
  filenumber = FileNumber;
  eventnumber = EventNumber;

  //Define Histgram
  DefineHist();
  Int_t nevent=tree->GetEntries();         
  tree->GetEntry(eventnumber);
  for(Int_t ich=0;ich<nhch;++ich){
    int ich0 = 63 - ich;
    if (ich > 63) ich0 = 64 + 127 - ich;
    //****Make Waveform****
    for(Int_t ibin=0;ibin<nbin;++ibin){
      hWaveform[ich]->  SetBinContent(nbin-ibin,trvecADC->at(ich).at(ibin));
    }
  } // ch loop
  // Kinetic_Energy->resize(Header[0]); PID->resize(Header[0]);
  // ParID->resize(Header[0]); trkID->resize(Header[0]);
  // Position->resize(Header[0],vector<double>(3,0)); Direction->resize(Header[0],vector<double>(3, 0));
  // GammaStopPos->resize(Header[6],vector<double>(3,0));
}


void Tracking(Bool_t Rawtrack,Bool_t NoRebintrack,  Bool_t NoiseReductrack){
  if(Rawtrack){
    for(Int_t ich=0;ich<nhch/2;ich++){
      for(Int_t ibin=0;ibin<nbin;ibin++){
	hxz -> SetBinContent(ich+1,ibin+1,hWaveform[ich]->GetBinContent(ibin+1));
	hyz -> SetBinContent(ich+1,ibin+1,hWaveform[nhch/2+ich]->GetBinContent(ibin+1));
      }                                                        
    }
  }
}


void DrawPlot(bool RawDisp,bool TrackDisp){
  bool rawwaveform=RawDisp;
  bool track=TrackDisp;

  if(rawwaveform){
    for(Int_t ic=0;ic<8;++ic){
      delete cWaveform[ic];
      cWaveform[ic] = new TCanvas(Form("cWaveform%d",ic),Form("cWaveform%d",ic),1600,800);
      cWaveform[ic] -> Divide(4,4);
    }
    
    for(Int_t ich=0; ich < nhch; ich++){
      cWaveform[ich/16] -> cd(ich%16+1);
      cWaveform[ich/16] -> SetGrid(2);
      
      // hWaveform[ich] -> GetXaxis() -> SetRangeUser(-50., TEND);
      // hWaveform[ich] -> GetYaxis() -> SetRangeUser(YMIN, YMAX);
      hWaveform[ich] -> SetStats(0);
      hWaveform[ich] -> Draw("hist");
      
    }
  }
   
  //Draw Track
  if(track){
    if(!gROOT->FindObjectAny("cTrack")){
      delete cTrack;
      cTrack = new TCanvas("cTrack","cTrack",800,400);
      cTrack->Divide(2,1);
    }

    cTrack->cd(1);
    //hxz_raw -> Draw("colz");
    hxz -> GetZaxis() -> SetRangeUser(-20.,300);
    gStyle->SetPalette(kRainBow);
    hxz -> Draw("colz");

    cTrack->cd(2);
    //hyz_raw -> Draw("colz");
    hyz -> GetZaxis() -> SetRangeUser(-20.,300);
    gStyle->SetPalette(kRainBow);
    hyz -> Draw("colz");
    
    
    cTrack->cd(1);
    delete latexhre;latexhre=NULL;                         
    latexhre=new TLatex(0.85,0.93,"");
    latexhre->SetNDC(1);
    latexhre->SetTextFont(102);
    latexhre->SetTextColor(kBlue);
    latexhre->SetTextSize(0.05);
    latexhre->SetTextAlign(31); 
    latexhre->SetText(0.85,0.93,Form("Event%d(%d)",eventnumber, trigCounter[0]));
    latexhre->Draw();  
  }
}


//Delete and redefine Histgram
void DefineHist(){
  //****************For Plot(Int_t ievent) *********************
  //Raw Waveform
  for(Int_t ich=0;ich<nhch;ich++){
    delete hWaveform[ich];
    hWaveform[ich] = new TH1D(Form("PedestalWaveform_ch%d",ich),Form("PedestalWaveform_ch%d;Time (#mus);ADC counts",ich),nbin ,-(Double_t)0.4*ptrig,0.4*(Double_t)(nbin-ptrig));
  }

  
  //*****************For Tracking()********************
  delete hxz;
  hxz = new TH2D(Form("xz Run%d File%d Event%d",runnumber,filenumber,eventnumber),"xzl;x (mm);z (mm)",nhch/2,-10,310,nbin,-(Double_t)0.4*ptrig,(nbin-ptrig)*0.4);
 
  delete hyz;
  hyz = new TH2D("yzl","yzl;y (mm);z (mm)",nhch/2,-10,310,nbin,-(Double_t)0.4*ptrig,(Double_t)(nbin-ptrig)*0.4);  
}

