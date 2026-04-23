#include "myheader.h"
#include <fftw3.h>
#include <vector>
#include <TSystem.h>

//const about Plot
const Int_t nch=256;
const Int_t nhch=128;
const Int_t nbin=1000;
const Int_t ptrig=125;
const Double_t timeconv = 0.4;
const Double_t ADCconv = (Double_t)2./pow(2,12);

//const about analysis
const Int_t bgch=16;
const Int_t PED_STARTBIN=0, PED_STOPBIN=125;
Double_t HGRMSbin=120.,HGRMSbin_nr=50., LGRMSbin=5., LGRMSbin_nr=5.;


TString inFileName;
TFile *fin;
TTree *tree;
TBranch *b_trigCounter,*b_trvecADC,*b_CLKCounter;
unsigned int trigCounter[4];
unsigned int CLKCounter[4];
std::vector<std::vector<double> >  *trvecADC;
Int_t Header[20],MHeaderI[20], Process;Double_t StopPos[3]={0}, InjecPos[3]={0}, HeaderD[20],MHeaderD[20];
vector<vector<double>> *Pos, *Dir;
vector<vector<double>> *GammaStopPos;
vector<double> *Kine, *Length;
vector<int> *PID, *trkID, *ParID, *GParID, *ParPID, *GParPID;
TLatex*latexhre;

//Global variables
Int_t runnumber;
Int_t filenumber;
Int_t eventnumber;

//TH1D 
TH1D *hPedN2;
//TH1D 
TH1D *hped,*hPedN,*hWaveform[nch],*hWaveform_nr[nch], *hWaveform_nr2[nch];

//TH2D
TH2D *hxz;
TH2D *hyz;


//Noise Reduction
const Int_t npx=1000;
const Int_t nBoard=4;
const Int_t nRebin=1;
//const Int_t snRebin=5;
const Int_t fitbin=1;
const Int_t nrbin=nbin/fitbin;
const Int_t nfit=nbin/fitbin;
const Double_t Fitregion=300.;
TH1D*hProX;
TH1D*hProY;
TH1D*hDistX[nBoard];
TH1D*hDistY[nBoard];
TH1D*hnrform[nBoard*2];
TF1*fDistX[nBoard];
TF1*fDistY[nBoard];
Int_t fittime=0;


TFile *outRootFile;
TTree *WaveTree; 
std::vector<std::vector<double> > trvecNRADC(nhch,std::vector<double>(nbin+1));

//Pedestal Estimation
TF1*fped;
TF1*fped_nr;

Int_t seed;
Int_t sigma;
Int_t dx=0, dy=0;
int mode =0;

  Double_t pedmean[nhch]={1129.02,1199.2,801.481,548.087,755.391,980.589,1026.92,1195.84,1048.04,908.585,1192.98,681.027,1165.85,850.385,779.237,1170.78,990.576,960.688,863.949,1124.11,1041.61,1349.93,1218.91,1110.95,1017.85,877.086,1246.31,956.836,942.081,1081.23,1008.96,981.264,1129.78,1382.82,1220.61,1077.7,1419.98,1131.8,1255.47,848.97,1118.43,1302.32,1480.59,848.19,1311.55,1152.19,1036.39,1291.89,1006.16,1057.62,1013.04,715.836,1199.07,892.052,1288.58,1105.02,1076.02,755.2,822.814,1143.31,1216.58,1047.7,1054.62,901.356,1038.55,972.701,1049.7,918.545,1216.04,1078.68,1208.05,1002.85,1123.32,1329.81,1219.8,1317.81,956.498,1700.07,1116.32,1066.46,1125.43,1234.29,1095.82,1191.21,1139.12,931.042,1254.02,1080.23,1051.15,1035.37,1340.63,1173.03,926.629,1180.19,1059.97,1173.35,1029.83,943.416,1043.1,867.078,852.112,1157.2,1189.81,1114.83,904.913,1084.47,1071.75,811.183,999.13,1185.61,1305.96,1205.43,1059.98,1089.58,1165.91,1091.92,1020.76,1204.84,1191.7,1240.49,1112.39,968.871,1118.52,883.354,888.35,1021.23,979.18,905.46};

int ddd=0;

void LoadRootFile(Int_t RunNumber, Int_t FileNumber,const char *fname="PION"){
  // seed = RunNumber;
  // sigma = FileNumber;
  runnumber = RunNumber;
  filenumber = FileNumber;
  std::cout << "sigma : " << sigma << std::endl;
  //inFileName = Form("./root/Sim/T98/GiBUU/Chain/G4_%d_sigma%d_x%d_y%d_%d.%s.root",runnumber,filenumber,dx,dy,ddd,fname);
  if(mode){
    inFileName = Form("./output/G4_%d_sigma%d_x%d_y%d_%d.%s.root",runnumber,filenumber,dx,dy,ddd,fname);
  }else{
    inFileName = Form("./output/G4_%d_sigma%d_x%d_y%d_%d.%s.root",fname,runnumber,filenumber,dx,dy,ddd,fname);
  }
  
  fin = TFile::Open(inFileName);
  delete tree;
  tree = (TTree*)fin->Get("WaveTree");
  tree->SetBranchAddress("trigCounter", &trigCounter, &b_trigCounter);
  tree->SetBranchAddress("CLKCounter", &CLKCounter, &b_CLKCounter);
  tree->SetBranchAddress("trvecADC", &trvecADC, &b_trvecADC);
  tree->SetBranchAddress("Header", Header);
  tree->SetBranchAddress("HeaderD", HeaderD);
  tree->SetBranchAddress("MHeaderI", MHeaderI);
  tree->SetBranchAddress("MHeaderD", MHeaderD);
  // tree->SetBranchAddress("Position", &Pos);
  // tree->SetBranchAddress("Direction", &Dir);
  //tree->SetBranchAddress("Kinetic_Energy", &Kine);
  // tree->SetBranchAddress("Length", &Length);
  // tree->SetBranchAddress("trkID", &trkID);
  // tree->SetBranchAddress("PID", &PID);
  // tree->SetBranchAddress("ParID", &ParID);
  // tree->SetBranchAddress("ParPID", &ParPID);
  // tree->SetBranchAddress("GParID", &GParID);
  // tree->SetBranchAddress("GParPID", &GParPID);
  // tree->SetBranchAddress("InjecPos", &InjecPos);
  // tree->SetBranchAddress("StopPos", &StopPos);
  // //tree->SetBranchAddress("GammaPosition", &GammaStopPos);
  // tree->SetBranchAddress("Process", &Process);
  
  std::cout << "Load " << inFileName << std::endl;
}


double nnn(double t, double t0, double w){
  return cos( 2.0 * TMath::Pi() * (t-t0)/w);              
}

// int main(int argc, char *argv[]){
//   int RunNumber=atoi(argv[1]);
//   int StartFileNumber=atoi(argv[2]);
//   int EndFileNumber=atoi(argv[3]);
//   std::cout << RunNumber << " " << StartFileNumber << " " << EndFileNumber << std::endl;
//   for(Int_t ifile=StartFileNumber; ifile <= EndFileNumber ; ifile++){
//     AddNoiseRootFile(RunNumber, ifile);
//   }
//   //std::cout <<"Finished " << std::endl;
//   outRootFile->Close();
//   //std::cout <<"Close" << std::endl;
// }


void AddNoiseGiBUU(Int_t RunNumber,
		   Int_t Sigma,
		   Double_t Amp,
		   Int_t Dx=0,
		   Int_t Dy=0,
		   Int_t ddd_ =0,
		   const char *fname="PION",
		   int mode_=0){
  seed = RunNumber;
  sigma = Sigma;
  dx = Dx;
  dy = Dy;
  ddd = ddd_;
  mode = mode_;

  // if(Amp==0){
  //   if(ddd==2){
  //     outRootFile = new TFile(Form("./root/Sim/%s/Chain/G4_%d_sigma%d_x%d_y%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(Amp),fname),"create") ;
  //   }else{
  //     outRootFile = new TFile(Form("./root/Sim/GiBUU/Chain/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(Amp),ddd,fname),"create") ;
  //     //outRootFile = new TFile(Form("./root/Sim/T98/GiBUU/Chain/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(Amp),ddd,fname),"create") ;
  //   }
  // }else{    
  //   if(ddd==2){
  //     outRootFile = new TFile(Form("./root/Sim/%s/Chain/G4_%d_sigma%d_x%d_y%d_%d.%s.root",seed, sigma,dx,dy,(Int_t)(100*Amp),fname),"create") ;
  //   }else{
  //     if(mode){
  // 	outRootFile = new TFile(Form("./root/Sim/GiBUU/Chain/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(100*Amp),ddd,fname),"create");
  //     }else{
  // 	outRootFile = new TFile(Form("./root/Sim/%s/Chain/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(100*Amp),ddd,fname),"create");
  //     }
  //     //outRootFile = new TFile(Form("./root/Sim/T98/GiBUU/Chain/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root", seed, sigma,dx,dy,(Int_t)(100*Amp),ddd,fname),"create") ;
  //   }
  // }
TString outDir;
if(ddd==2){
  outDir = Form("./output");
}else{
  if(mode){
    outDir = Form("./output");
  }else{
    outDir = Form("./%s", fname);
  }
}
gSystem->mkdir(outDir, true);

if(Amp==0){
  if(ddd==2){
    outRootFile = new TFile(Form("%s/G4_%d_sigma%d_x%d_y%d_%d.%s.root",
                                 outDir.Data(), seed, sigma, dx, dy, (Int_t)(Amp), fname),
                            "create");
  }else{
    outRootFile = new TFile(Form("%s/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root",
                                 outDir.Data(), seed, sigma, dx, dy, (Int_t)(Amp), ddd, fname),
                            "create");
  }
}else{
  if(ddd==2){
    outRootFile = new TFile(Form("%s/G4_%d_sigma%d_x%d_y%d_%d.%s.root",
                                 outDir.Data(), seed, sigma, dx, dy, (Int_t)(100*Amp), fname),
                            "create");
  }else{
    outRootFile = new TFile(Form("%s/G4_%d_sigma%d_x%d_y%d_%d_%d.%s.root",
                                 outDir.Data(), seed, sigma, dx, dy, (Int_t)(100*Amp), ddd, fname),
                            "create");
  }
}

  
  cout <<"A" <<endl; 
  WaveTree = new TTree("WaveTree","WaveTree");
  WaveTree -> Branch("trigCounter", &trigCounter, "trigCounter[4]/I");
  WaveTree -> Branch("CLKCounter",  &CLKCounter,  "CLKCounter[4]/I"); 
  WaveTree -> Branch("trvecADC",  &trvecNRADC);
  WaveTree->Branch("Header", &Header,"Header[20]/I");
  WaveTree->Branch("HeaderD", &HeaderD,"HeaderD[20]/D");
  WaveTree->Branch("MHeaderI", &MHeaderI,"MHeaderI[20]/I");
  WaveTree->Branch("MHeaderD", &MHeaderD,"MHeaderD[20]/D");
  //WaveTree->Branch("Position", &Pos);
  //WaveTree->Branch("Direction", &Dir);
  //WaveTree->Branch("Kinetic_Energy", &Kine);
  //WaveTree->Branch("Length", &Length);
  //WaveTree->Branch("trkID", &trkID);
  // WaveTree->Branch("ParID", &ParID);
  // WaveTree->Branch("GParID", &GParID);
  //WaveTree->Branch("PID", &PID);
  // WaveTree->Branch("ParPID", &ParPID);
  // WaveTree->Branch("GParPID", &GParPID);
  //WaveTree->Branch("StopPos", &StopPos, "StopPos[3]/D");
  //WaveTree->Branch("InjecPos", &InjecPos, "InjecPos[3]/D");
  //WaveTree->Branch("GammaPosition", &GammaStopPos);
  //WaveTree->Branch("Process", &Process);
  
  //std::vector<std::vector<double> >().swap(trvecNRADC);
  trvecNRADC.resize(nch, std::vector<double>(nbin));
  LoadRootFile(seed,0,fname);
  //Int_t nevent = 10;
  Int_t nevent = tree->GetEntries();
  Double_t t0x=gRandom->Uniform(0,10);     
  Double_t wwx=gRandom->Gaus(10,1); 
  Double_t t0y=gRandom->Uniform(0,10);     
  Double_t wwy=gRandom->Gaus(10,1);
  for(Int_t ievent=0; ievent < nevent; ievent++){
    tree->GetEntry(ievent);
    for(Int_t ich=0; ich < nhch; ich++){
      for(Int_t ibin=0; ibin < nbin; ibin++){
	Double_t tt = (Double_t)ibin*timeconv-50.-timeconv/2.;
	if(ich<nhch/2){
	  if(sigma!=0)trvecNRADC[ich*2][ibin] =(trvecADC->at(ich*2).at(ibin)-pedmean[ich])*Amp+gRandom->Gaus(0,(Double_t)sigma)+100.*nnn(tt,t0x,wwx)+pedmean[ich];
	  if(sigma==0)trvecNRADC[ich*2][ibin] =(trvecADC->at(ich*2).at(ibin)-pedmean[ich])*Amp+pedmean[ich];
	  //std::cout <<"Sigma : " << sigma  << " Noise : " << gRandom->Gaus(0,(Double_t)sigma) + 100.*nnn(tt,t0x,wwx) << std::endl;;
	  //std::cout << "Random Noise : " << gRandom->Gaus(0,(Double_t)sigma) << std::endl;;
	}else{
	  if(sigma!=0)trvecNRADC[ich*2][ibin] =(trvecADC->at(ich*2).at(ibin)-pedmean[ich])*Amp+gRandom->Gaus(0,(Double_t)sigma)+100.*nnn(tt,t0x,wwx)+pedmean[ich];
	  if(sigma==0)trvecNRADC[ich*2][ibin] =(trvecADC->at(ich*2).at(ibin)-pedmean[ich])*Amp+pedmean[ich];
	}
      }
    }
    WaveTree->Fill();
    //if(ievent%(nevent/100)==0)std::cout << "\rEvent : " << ievent << " Finished! " << std::flush;
  }
  outRootFile->cd();
  WaveTree->Write();
  if(gROOT->FindObjectAny("WaveTree")){      
    delete  WaveTree;
    //std::cout <<"Delete NRtree" << std::endl;
  }
  if(gROOT->FindObjectAny("WaveTree")){      
    delete  tree;
    //std::cout <<"Delete tree" << std::endl;
  }
  outRootFile->Close();
  fin->Close();
  delete outRootFile;
}

//Delete and redefine Histgram
void DefineHist(){
  //****************For Plot(Int_t ievent) *********************
  //Raw Waveform
  for(Int_t ich=0;ich<nhch;ich++){
    delete hWaveform[ich];
    hWaveform[ich] = new TH1D(Form("PedestalWaveform_ch%d",ich),Form("PedestalWaveform_ch%d;Time (#mus);ADC counts",ich),nbin ,-(Double_t)0.4*ptrig,0.4*(Double_t)(nbin-ptrig));
    delete hWaveform_nr[ich];
    hWaveform_nr[ich] = new TH1D(Form("hWaveform_nr%d",ich+1),Form("hWaveform_noise_reducted_ch%d;Time (#mus);ADC counts",ich+1),nbin/fitbin,-(Double_t)0.4*ptrig,0.4*(Double_t)(nbin-ptrig));
    delete hWaveform_nr2[ich];
    hWaveform_nr2[ich] = new TH1D(Form("hWaveform_nr2_%d",ich+1),Form("hWaveform_noise_reducted2_ch%d;Time (#mus);ADC counts",ich+1),nbin/fitbin,-(Double_t)0.4*ptrig,0.4*(Double_t)(nbin-ptrig));

  }
  //****************For MakeWaveform()***********************
  //Noise Reducted Waveform
  for(Int_t iboard=0;iboard<nBoard*2;iboard++){
    delete hnrform[iboard];
    hnrform[iboard] = new TH1D(Form("hnrform_board%d",iboard),Form("hnrform_Board%d;Time (#mus);ADC counts",iboard),nbin/fitbin,-(Double_t)0.4*ptrig,0.4*(Double_t)(nbin-ptrig));
  }  

  //*****************For Tracking()********************
  delete hxz;
  hxz = new TH2D(Form("xz Run%d File%d Event%d",runnumber,filenumber,eventnumber),"xzl;x (mm);z (#mus)",nhch/2,-0.5,63.5,nbin,-(Double_t)0.4*ptrig,(nbin-ptrig)*0.4);
 
  delete hyz;
  hyz = new TH2D("yzl","yzl;y (mm);z (#mus)",nhch/2,-0.5,63.5,nbin,-(Double_t)0.4*ptrig,(Double_t)(nbin-ptrig)*0.4);  

}

