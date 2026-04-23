#include "myheader.h"
#include <fftw3.h>
#include <vector>

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
Int_t Header[20],MHeaderI[20], Process;Double_t StopPos[3]={0}, InjecPos[3]={0}, HeaderD[20],MHeaderD[20];
vector<vector<double>> *Pos, *Dir;
vector<vector<double>> *GammaStopPos;
vector<double> *Kine, *Length;
vector<int> *PID, *trkID, *ParID, *GParID, *ParPID, *GParPID;


//Pedestal Estimation
TF1*fped;
TF1*fped_nr;

Int_t sigma;
Int_t seed;
Int_t amp;
Int_t dx=0,dy=0;
int mode =0;

void LoadRootFile(Int_t RunNumber, Int_t FileNumber, const char*fname="PION"){
  seed = RunNumber;
  sigma = FileNumber;
  runnumber = RunNumber;
  filenumber = FileNumber;
  
  if(mode){
    inFileName = Form("./output/G4_%d_sigma%d_x%d_y%d_%d_1.%s.root",seed,sigma,dx,dy,amp,fname);
  }else{
    inFileName = Form("./%s/Chain/G4_%d_sigma%d_x%d_y%d_%d_1.%s.root",fname,seed,sigma,dx,dy,amp,fname);
  }
  //dddが固定されていることに注意
  
  //inFileName = Form("./root/Sim/T98/%s/G4_%d_sigma%d_x%d_y%d_%d.%s.root",fname,seed,sigma,dx,dy,amp,fname);
  
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
  //tree->SetBranchAddress("ParID", &ParID);
  //tree->SetBranchAddress("ParPID", &ParPID);
  //tree->SetBranchAddress("GParID", &GParID);
  //tree->SetBranchAddress("GParPID", &GParPID);
  // tree->SetBranchAddress("InjecPos", &InjecPos);
  // tree->SetBranchAddress("StopPos", &StopPos);
  //tree->SetBranchAddress("GammaPosition", &GammaStopPos);
  //tree->SetBranchAddress("Process", &Process);
  
  std::cout << "Load " << inFileName << std::endl;
}

// int main(int argc, char *argv[]){
//   int RunNumber=atoi(argv[1]);
//   int StartFileNumber=atoi(argv[2]);
//   int EndFileNumber=atoi(argv[3]);
//   std::cout << RunNumber << " " << StartFileNumber << " " << EndFileNumber << std::endl;
//   for(Int_t ifile=StartFileNumber; ifile <= EndFileNumber ; ifile++){
//     MakeSingleRootFile(RunNumber, ifile);
//   }
//   //std::cout <<"Finished " << std::endl;
//   outRootFile->Close();
//   //std::cout <<"Close" << std::endl;
// }
void MakeWaveform2(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber);

void MakeNoiseGiBUU(Int_t RunNumber,
		    Int_t FileNumber,
		    Double_t Amp,
		    Int_t Dx =0,
		    Int_t Dy=0,
		    const char*fname="PION",
		    int mode_=0){
  // if(gROOT->FindObjectAny("WaveTree")){
  //   delete ;
  //   std::cout <<"Delete" << std::endl;
  // }
  seed = RunNumber;
  sigma = FileNumber;
  if(Amp==1){
    amp = Amp;
  }else{
    amp=(Int_t)100*Amp;
  }
  std::cout <<"Amp : " << Amp << std::endl;
  dx = Dx;
  dy = Dy;
  mode = mode_;
  
  //outRootFile = new TFile(Form("./root/Sim/T98/%s/G4_%d_sigma%d_x%d_y%d_%d.NR.%s.root",fname, seed, sigma, dx,dy,amp,fname),"create") ;
  if(mode){
    outRootFile = new TFile(Form("./output/G4_%d_sigma%d_x%d_y%d_%d_1.NR.%s.root", seed, sigma, dx,dy,amp,fname),"create") ;
  }else{
    outRootFile = new TFile(Form("./%s/Chain/G4_%d_sigma%d_x%d_y%d_%d_1.NR.%s.root", fname, seed, sigma, dx,dy,amp,fname),"create") ;
  }
  
  WaveTree = new TTree("NRWaveTree","NRWaveTree");
  WaveTree -> Branch("trigCounter", &trigCounter, "trigCounter[4]/I");
  WaveTree -> Branch("CLKCounter",  &CLKCounter,  "CLKCounter[4]/I"); 
  WaveTree -> Branch("trvecADC",  &trvecNRADC);
  WaveTree->Branch("Header", &Header,"Header[20]/I");
  WaveTree->Branch("HeaderD", &HeaderD,"HeaderD[20]/D");
  WaveTree->Branch("MHeaderI", &MHeaderI,"MHeaderI[20]/I");
  WaveTree->Branch("MHeaderD", &MHeaderD,"MHeaderD[20]/D");
  // WaveTree->Branch("Position", &Pos);
  // WaveTree->Branch("Direction", &Dir);
  // WaveTree->Branch("Kinetic_Energy", &Kine);
  // WaveTree->Branch("Length", &Length);
  // WaveTree->Branch("trkID", &trkID);
  // WaveTree->Branch("ParID", &ParID);
  // WaveTree->Branch("GParID", &GParID);
  // WaveTree->Branch("PID", &PID);
  // WaveTree->Branch("ParPID", &ParPID);
  // WaveTree->Branch("GParPID", &GParPID);
  // WaveTree->Branch("StopPos", &StopPos, "StopPos[3]/D");
  // WaveTree->Branch("InjecPos", &InjecPos, "InjecPos[3]/D");
  // WaveTree->Branch("GammaPosition", &GammaStopPos);
  // WaveTree->Branch("Process", &Process);
  
  //std::vector<std::vector<double> >().swap(trvecNRADC);
  trvecNRADC.resize(nhch, std::vector<double>(nbin/fitbin+1));
  LoadRootFile(RunNumber,FileNumber,fname);
  //Int_t nevent = 10;
  Int_t nevent = tree->GetEntries();
  for(Int_t ievent=0; ievent < nevent; ievent++){
    MakeWaveform2(runnumber, filenumber, ievent);
    //std::cout <<"Make Waveform" << std::endl;
    for(Int_t ich=0; ich < nhch; ich++){
      for(Int_t ibin=0; ibin < nbin/fitbin+1; ibin++){
	if(sigma==0)trvecNRADC[ich][ibin] = hWaveform[ich]->GetBinContent(ibin+1);
	if(sigma!=0)trvecNRADC[ich][ibin] = hWaveform_nr2[ich]->GetBinContent(ibin+1);
      }
    }
    WaveTree->Fill();
    //if(ievent%(nevent/100)==0)std::cout << "\rEvent : " << ievent << " Finished! " << std::flush;
  }
  std::cout<<std::endl;
  
  outRootFile->cd();
  WaveTree->Write();
  if(gROOT->FindObjectAny("NRWaveTree")){      
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


void MakeWaveform2(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber){
  runnumber = RunNumber;
  filenumber = FileNumber;
  eventnumber = EventNumber;
  
  Double_t ped;
  Double_t pedRMS;
  
  Double_t pedmean[nhch]={1129.02,1199.2,801.481,548.087,755.391,980.589,1026.92,1195.84,1048.04,908.585,1192.98,681.027,1165.85,850.385,779.237,1170.78,990.576,960.688,863.949,1124.11,1041.61,1349.93,1218.91,1110.95,1017.85,877.086,1246.31,956.836,942.081,1081.23,1008.96,981.264,1129.78,1382.82,1220.61,1077.7,1419.98,1131.8,1255.47,848.97,1118.43,1302.32,1480.59,848.19,1311.55,1152.19,1036.39,1291.89,1006.16,1057.62,1013.04,715.836,1199.07,892.052,1288.58,1105.02,1076.02,755.2,822.814,1143.31,1216.58,1047.7,1054.62,901.356,1038.55,972.701,1049.7,918.545,1216.04,1078.68,1208.05,1002.85,1123.32,1329.81,1219.8,1317.81,956.498,1700.07,1116.32,1066.46,1125.43,1234.29,1095.82,1191.21,1139.12,931.042,1254.02,1080.23,1051.15,1035.37,1340.63,1173.03,926.629,1180.19,1059.97,1173.35,1029.83,943.416,1043.1,867.078,852.112,1157.2,1189.81,1114.83,904.913,1084.47,1071.75,811.183,999.13,1185.61,1305.96,1205.43,1059.98,1089.58,1165.91,1091.92,1020.76,1204.84,1191.7,1240.49,1112.39,968.871,1118.52,883.354,888.35,1021.23,979.18,905.46};
  double scale2[16] = {0,0, -0.6,0, 0.30, 0.50, 0.64, 0.76, 0.85, 0.93, 1.00, 1.06,1.11, 1.16, 1.21, 1.25};
  
  //Define Histgram
  DefineHist();
  Int_t nevent=tree->GetEntries();         
  tree->GetEntry(eventnumber);
  for(Int_t ich=0;ich<nhch;++ich){
    int ich0 = 63 - ich;
    if (ich > 63) ich0 = 64 + 127 - ich;     
    //****Make Waveform****
    for(Int_t ibin=0;ibin<nbin;++ibin){
      //hWaveform[ich]->  SetBinContent(ibin+1,trvecADC->at(ich0*2).at(ibin));
      hWaveform[ich]->  SetBinContent(ibin+1,trvecADC->at(ich0*2).at(ibin)-pedmean[ich0]);
    }
    delete hped;
    hped=NULL;
  } // ch loop
  
  if(sigma!=0){
    //*****************************Noise Reduction***********************************
    //****Make 2DTracking*****
    Tracking(true,false,false);

    //***************Bin Loop************ 
    for(Int_t ifit=0; ifit<nfit; ifit++){
      //Each Bin Distribution
      for(Int_t iboard=0; iboard < 4; iboard++){
	delete hDistX[iboard];
	hDistX[iboard] = new TH1D(Form("hDistX%d",iboard),Form("X ch distribution Board %d",iboard),(Int_t)2*Fitregion/1,-1.*Fitregion,Fitregion);
	delete hDistY[iboard];
	hDistY[iboard] = new TH1D(Form("hDistY%d",iboard),Form("Y ch distribution Board %d",iboard),(Int_t)2*Fitregion/1,-1.*Fitregion,Fitregion);
	delete fDistX[iboard];
	fDistX[iboard] = new TF1(Form("fDistX%d",iboard),"gaus",-1.*Fitregion,Fitregion);
	delete fDistY[iboard];
	fDistY[iboard] = new TF1(Form("fDistY%d",iboard),"gaus",-1.*Fitregion,Fitregion);
      }

      for(Int_t ibin=0;ibin<fitbin;ibin++){
	//Projected Waveform
	delete hProX;
	hProX = hxz -> ProjectionX(Form("hProjectionX%d",ibin+ifit*fitbin), 1+ibin+ifit*fitbin, 1+ibin+ifit*fitbin);
	delete hProY;
	hProY = hyz -> ProjectionX(Form("hProjectionY%d",ibin+ifit*fitbin), 1+ibin+ifit*fitbin, 1+ibin+ifit*fitbin);
      
	for(Int_t ich=2; ich<62; ich++){
	  //Odd or Even + 1,3 or 2,4 	
	  if(ich<32){
	    if(ich%2==0)  hDistX[0] -> Fill(hProX->GetBinContent(1+ich));
	    if(ich%2==1)  hDistX[1] -> Fill(hProX->GetBinContent(1+ich));
	    if(ich%2==0)  hDistY[0] -> Fill(hProY->GetBinContent(1+ich));
	    if(ich%2==1)  hDistY[1] -> Fill(hProY->GetBinContent(1+ich));	
	  }else{
	    if(ich%2==0)  hDistX[2] -> Fill(hProX->GetBinContent(1+ich));
	    if(ich%2==1)  hDistX[3] -> Fill(hProX->GetBinContent(1+ich));
	    if(ich%2==0)  hDistY[2] -> Fill(hProY->GetBinContent(1+ich));
	    if(ich%2==1)  hDistY[3] -> Fill(hProY->GetBinContent(1+ich));
	  }
	}
      }
      for(Int_t iboard=0; iboard < 4; iboard++){

	// X channel common noise
	int ibin = 1;
	int nsum = 0;
	int isum = 0;
	double asum = 999.;
	double ss   = sigma;

	while (hDistX[iboard]->GetBinCenter(ibin) < asum + 3.*ss){
	  nsum = nsum + hDistX[iboard]->GetBinContent(ibin);
	  if (isum == 0 && nsum >= 2) {
	    isum = ibin;
	    asum = hDistX[iboard]->GetBinCenter(ibin);
	  }
	  ibin++;
	}
	double amaxx = asum + scale2[nsum]*ss;

	// Y channel common noise
	ibin = 1;
	nsum = 0;
	isum = 0;
	asum = 999.;
	while (hDistY[iboard]->GetBinCenter(ibin) < asum + 3.*ss){
	  nsum = nsum + hDistY[iboard]->GetBinContent(ibin);
	  if (isum == 0 && nsum >= 2) {
	    isum = ibin;
	    asum = hDistY[iboard]->GetBinCenter(ibin);
	  }
	  ibin++;
	}
	double amaxy = asum + scale2[nsum]*ss;
	hnrform[iboard] -> SetBinContent(1+ifit, amaxx);
	hnrform[iboard+4] -> SetBinContent(1+ifit, amaxy);
      }
    }
  
  
    for(Int_t ich=0; ich < nhch; ich++){
      hWaveform_nr[ich] -> Add(hWaveform[ich],+1.0);
      hWaveform_nr[ich] -> Add(hnrform[ich%2+(ich/32)*2],-1.0);

      //Pedestal Re-calculation
      delete hPedN;     
      hPedN = new TH1D(Form("hPedN_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),10*nbin,-nbin*5.-0.5,nbin*5.-0.5);  
      for(Int_t ibin=PED_STARTBIN;ibin<=PED_STOPBIN;++ibin){
	hPedN -> Fill(hWaveform_nr[ich]->GetBinContent(ibin+1));
      } 
      Double_t sigmafit = 10.;
      Double_t mean = hPedN->GetMean();   
      delete fped;
      fped = new TF1("fPedestal", "gaus", mean - 5 * sigmafit,  mean + 2 * sigmafit);
      fped->SetParameters(hPedN->GetMaximum(), mean, sigmafit);
      //fped->SetParameters(hPedN->GetMaximum(), 0, sigmafit);
      // fped = new TF1("fPedestal", "gaus", mean - 5 * sigmafit, mean + 2 * sigmafit);
      // fped->SetParameters(hPedN->GetMaximum(), mean, sigmafit);
    
      hPedN->Fit(fped, "0LRQ");
      Double_t fitped = fped->GetParameter(1);
      for(Int_t ibin=0;ibin<nbin;++ibin){
	hWaveform_nr2[ich]->  SetBinContent(ibin+1,hWaveform_nr[ich]->GetBinContent(ibin+1)-fitped);
      }
      delete hPedN;
      hPedN = new TH1D(Form("hPedN_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),nbin*10.,-nbin*5.-0.5,nbin*5.-0.5);
      for(Int_t ibin=PED_STARTBIN/fitbin;ibin<=PED_STOPBIN/fitbin;ibin++){
	hPedN -> Fill(hWaveform_nr2[ich]->GetBinContent(ibin+1));
      }
      hWaveform_nr2[ich]->SetBinContent(nbin/fitbin+1,hPedN->GetRMS());
      hWaveform_nr2[ich]->SetBinContent(0,hPedN->GetMean());
    }
  }
}




void MakeWaveform(Int_t RunNumber, Int_t FileNumber,Int_t EventNumber){
  runnumber = RunNumber;
  filenumber = FileNumber;
  eventnumber = EventNumber;

  Double_t ped;
  Double_t pedRMS;

  double scale2[16] = {0,0, -0.6,0, 0.30, 0.50, 0.64, 0.76, 0.85, 0.93, 1.00, 1.06,1.11, 1.16, 1.21, 1.25};

  //Define Histgram
  //std::cout <<"Delete NRtree" << std::endl;
  DefineHist();
  //std::cout <<"Delete NRtree" << std::endl;
  Int_t nevent=tree->GetEntries();         
  tree->GetEntry(eventnumber);
  Double_t sigma = 50.;
  for(Int_t ich=0;ich<nhch;++ich){
    int ich0 = 63 - ich;
    if (ich > 63) ich0 = 64 + 127 - ich;
    //Pedestal Mean Calculation
    delete hPedN;
    hPedN = new TH1D(Form("hPedN_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),nbin,-nbin*5.+0.5,nbin*5.-0.5);
    for(Int_t ibin=PED_STARTBIN;ibin<=PED_STOPBIN;++ibin){
      hPedN -> Fill(trvecADC->at(ich0*2).at(ibin+1));
    }
    
    Double_t mean = hPedN->GetMean();
    delete fped;
    fped = new TF1("fPedestal", "gaus", mean - 6 * sigma, mean + 2 * sigma);
    fped->SetParameters(hPedN->GetMaximum(), mean, sigma);
    
    hPedN->Fit(fped, "0LRQ");
    // if(ich==15){
    //   hcheck[0] -> Add(hPedN,+1.0);
    //   fcheck[0] = new TF1("copy", "gaus", mean - 2 * sigma, mean + 2 * sigma);
    //   for(Int_t ip = 0; ip < 3; ip++) {
    // 	fcheck[0]->SetParameter(ip, fped->GetParameter(ip));
    //   }
    // }
    
    Double_t fitped = fped->GetParameter(1);
    //std::cout << "Ch : "  << ich<< " Mean : " << hPedN->GetMean() << " Fitted Gaussian Mean: " << fitped << std::endl;
    
    ped = hPedN->GetMean();
    pedRMS=hPedN->GetRMS();
    
    hWaveform[ich] -> SetBinContent(0,ped);
    hWaveform[ich] -> SetBinContent(nbin+1,pedRMS);
    
    //****Make Waveform****
    for(Int_t ibin=0;ibin<nbin;++ibin){
      //hWaveform[ich]->  SetBinContent(ibin+1,trvecADC->at(ich0*2).at(ibin)-ped);
      hWaveform[ich]->  SetBinContent(ibin+1,trvecADC->at(ich0*2).at(ibin)-fitped);
    }

    delete hPedN;
    hPedN = new TH1D(Form("hPedN_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),nbin*10.,-nbin*5.+0.5,nbin*5.-0.5);
    // delete hPedN2;
    // hPedN2= new TH1D(Form("hPedN2_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),nbin*10.,-nbin*5.,nbin*5.);
    
    for(Int_t ibin=PED_STARTBIN;ibin<=PED_STOPBIN;++ibin){
      hPedN -> Fill(trvecADC->at(ich0*2).at(ibin+1)-fitped);
      //hPedN -> Fill(trvecADC->at(ich0*2).at(ibin+1)-ped);
      //hPedN2 -> Fill(trvecADC->at(ich0*2).at(ibin+1)-ped);
    }
    
    ped = hPedN->GetMean();
    pedRMS=hPedN->GetRMS();
    
    // tcheck[0]->SetPoint(ich,ich+1,pedRMS);
    // tcheck[1]->SetPoint(ich,ich+1,ped);
    
    // ped = hPedN2->GetMean();
    // pedRMS = hPedN2->GetRMS();
    
    
    delete hped;
    hped=NULL;
  } // ch loop

  //*****************************Noise Reduction***********************************
  //****Make 2DTracking*****
  Tracking(true,false,false);
  
  //TCanvas *c1 = new TCanvas("c1","c1",1600,800);
  //c1->Divide(4,2);
  //c1->Print("test.pdf[");

  //***************Bin Loop************ 
  for(Int_t ifit=0; ifit<nfit; ifit++){
    //Each Bin Distribution
    for(Int_t iboard=0; iboard < 4; iboard++){
      delete hDistX[iboard];
      hDistX[iboard] = new TH1D(Form("hDistX%d",iboard),Form("X ch distribution Board %d",iboard),(Int_t)2*Fitregion/1,-1.*Fitregion,Fitregion);
      delete hDistY[iboard];
      hDistY[iboard] = new TH1D(Form("hDistY%d",iboard),Form("Y ch distribution Board %d",iboard),(Int_t)2*Fitregion/1,-1.*Fitregion,Fitregion);
      delete fDistX[iboard];
      fDistX[iboard] = new TF1(Form("fDistX%d",iboard),"gaus",-1.*Fitregion,Fitregion);
      delete fDistY[iboard];
      fDistY[iboard] = new TF1(Form("fDistY%d",iboard),"gaus",-1.*Fitregion,Fitregion);
    }

    for(Int_t ibin=0;ibin<fitbin;ibin++){
      //Projected Waveform
      delete hProX;
      hProX = hxz -> ProjectionX(Form("hProjectionX%d",ibin+ifit*fitbin), 1+ibin+ifit*fitbin, 1+ibin+ifit*fitbin);
      delete hProY;
      hProY = hyz -> ProjectionX(Form("hProjectionY%d",ibin+ifit*fitbin), 1+ibin+ifit*fitbin, 1+ibin+ifit*fitbin);
      
      for(Int_t ich=2; ich<62; ich++){
	//Odd or Even + 1,3 or 2,4 	
	if(ich<32){
	  if(ich%2==0)  hDistX[0] -> Fill(hProX->GetBinContent(1+ich));
	  if(ich%2==1)  hDistX[1] -> Fill(hProX->GetBinContent(1+ich));
	  if(ich%2==0)  hDistY[0] -> Fill(hProY->GetBinContent(1+ich));
	  if(ich%2==1)  hDistY[1] -> Fill(hProY->GetBinContent(1+ich));	
	}else{
	  if(ich%2==0)  hDistX[2] -> Fill(hProX->GetBinContent(1+ich));
	  if(ich%2==1)  hDistX[3] -> Fill(hProX->GetBinContent(1+ich));
	  if(ich%2==0)  hDistY[2] -> Fill(hProY->GetBinContent(1+ich));
	  if(ich%2==1)  hDistY[3] -> Fill(hProY->GetBinContent(1+ich));
	}
      }
    }
    for(Int_t iboard=0; iboard < 4; iboard++){
      // X channel common noise
      int ibin = 1;
      int nsum = 0;
      int isum = 0;
      double asum = 999.;
      double ss   = 6.;
      while (hDistX[iboard]->GetBinCenter(ibin) < asum + 3.*ss){
	nsum = nsum + hDistX[iboard]->GetBinContent(ibin);
	if (isum == 0 && nsum >= 2) {
	  isum = ibin;
	  asum = hDistX[iboard]->GetBinCenter(ibin);
	}
	ibin++;
      }
      double amaxx = asum + scale2[nsum]*ss;

      // Y channel common noise
      ibin = 1;
      nsum = 0;
      isum = 0;
      asum = 999.;
      while (hDistY[iboard]->GetBinCenter(ibin) < asum + 3.*ss){
	nsum = nsum + hDistY[iboard]->GetBinContent(ibin);
	if (isum == 0 && nsum >= 2) {
	  isum = ibin;
	  asum = hDistY[iboard]->GetBinCenter(ibin);
	}
	ibin++;
      }
      double amaxy = asum + scale2[nsum]*ss;
      hnrform[iboard] -> SetBinContent(1+ifit, amaxx);
      hnrform[iboard+4] -> SetBinContent(1+ifit, amaxy);
    }
    //c1->Update();
    //c1->Print("test.pdf");
  }
  //c1->Print("test.pdf]");
  
  for(Int_t ich=0; ich < nhch; ich++){
    //Rebin and Scaling
    //hWaveform_raw[ich] -> Add(hWaveform[ich],+1.0);
    hWaveform[ich] -> Rebin(fitbin);
    hWaveform[ich] -> Scale(1./(Double_t)fitbin);
  }  

  for(Int_t ich=0; ich < nhch; ich++){
    hWaveform_nr[ich] -> Add(hWaveform[ich],+1.0);
    // if (ich<32&&ich%2==0) hWaveform_nr[ich] -> Add(hWaveform[4],-1.0);
    // else if (ich<32&&ich%2==1) hWaveform_nr[ich] -> Add(hWaveform[3],-1.0);
    // else if (ich<64&&ich%2==0) hWaveform_nr[ich] -> Add(hWaveform[60],-0.95);
    // else if (ich<64&&ich%2==1) hWaveform_nr[ich] -> Add(hWaveform[61],-0.9);
    // else hWaveform_nr[ich] -> Add(hnrform[ich%2+(ich/32)*2],-1.0);
    hWaveform_nr[ich] -> Add(hnrform[ich%2+(ich/32)*2],-1.0);
    delete hPedN;
    hPedN = new TH1D(Form("hPedN_ch%d",ich),Form("hPedN_ch%d;ADC counts (#mus);Number of Bin",ich),nbin*10.,-nbin*5.,nbin*5.);
    for(Int_t ibin=PED_STARTBIN/fitbin;ibin<=PED_STOPBIN/fitbin;ibin++){
      hPedN -> Fill(hWaveform_nr[ich]->GetBinContent(ibin+1));
    }
    pedRMS=hPedN->GetRMS();
    ped = hPedN->GetMean();
    // hWaveform_nr[ich]->SetBinContent(nbin/fitbin+1,pedRMS);
    // hWaveform_nr[ich] -> SetBinContent(0,ped);
    // tcheck[2]->SetPoint(ich,ich+1,pedRMS);
    // tcheck[3]->SetPoint(ich,ich+1,ped);
    //tcheck[2]->SetPoint(ich,ich+1,5*pedRMS);

    delete fped_nr;
    fped_nr = new TF1("fPedestal", "gaus", -sigma, sigma);
    fped_nr->SetParameters(hPedN->GetMaximum(), 0., sigma);
    
    hPedN->Fit(fped_nr, "0LRQ");
    
    // if(ich==15){
    //   hcheck[2] -> Add(hPedN,+1.0);
    //   fcheck[2] = new TF1("copy", "gaus", -2*sigma, 2*sigma);
    //   for(Int_t ip = 0; ip < 3; ip++) {
    // 	fcheck[2]->SetParameter(ip, fped_nr->GetParameter(ip));
    //   }
    // }
    hWaveform_nr[ich]->SetBinContent(nbin/fitbin+1,fped_nr->GetParameter(2));
    hWaveform_nr[ich] -> SetBinContent(0,fped_nr->GetParameter(1));
       
    //Double_t fitped = fped->GetParameter(1);
    //std::cout << "Ch : "  << ich<< " Mean : " << hPedN->GetMean() << " Fitted Gaussian Mean: " << fped_nr->GetParameter(1) << std::endl;
  }
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

  // if(NoRebintrack){
  //   for(Int_t ich=0;ich<nhch/2;ich++){
  //     for(Int_t ibin=0;ibin<nbin;ibin++){
  // 	hxz_raw -> SetBinContent(ich+1,ibin+1,hWaveform_raw[ich]->GetBinContent(ibin+1));
  // 	hyz_raw -> SetBinContent(ich+1,ibin+1,hWaveform_raw[nhch/2+ich]->GetBinContent(ibin+1));
  //     }                                                        
  //   }
  // }
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

