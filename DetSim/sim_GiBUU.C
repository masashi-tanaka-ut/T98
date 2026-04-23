#include "myheader.h"
#include <TMath.h>
#include <TVector3.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <cstdio>
#include <iostream>

static const int kNChan = 64;

double dstep;
TVector3 r0, r, rr;
TH1D* hx[kNChan];
TH1D* hy[kNChan];
TH1D* hx2[kNChan];
TH1D* hy2[kNChan];

TH2D* hhx;
TH2D* hhy;
TNtuple* n1;
int ichx, ichy;
double t0;
TF1* f1;
double ss[50];

TVector3 rx1(0.3, 0.3, 310.);
TVector3 rx2(4.2, 0.3, 310.);
TVector3 rx3(0.3, 4.2, 310.);
TVector3 ry1(4.7, 4.7, 310.);
TVector3 ry2(4.7, 0.8, 310.);
TVector3 ry3(0.8, 4.7, 310.);

  Double_t pedmean[128]={1129.02,1199.2,801.481,548.087,755.391,980.589,1026.92,1195.84,1048.04,908.585,1192.98,681.027,1165.85,850.385,779.237,1170.78,990.576,960.688,863.949,1124.11,1041.61,1349.93,1218.91,1110.95,1017.85,877.086,1246.31,956.836,942.081,1081.23,1008.96,981.264,1129.78,1382.82,1220.61,1077.7,1419.98,1131.8,1255.47,848.97,1118.43,1302.32,1480.59,848.19,1311.55,1152.19,1036.39,1291.89,1006.16,1057.62,1013.04,715.836,1199.07,892.052,1288.58,1105.02,1076.02,755.2,822.814,1143.31,1216.58,1047.7,1054.62,901.356,1038.55,972.701,1049.7,918.545,1216.04,1078.68,1208.05,1002.85,1123.32,1329.81,1219.8,1317.81,956.498,1700.07,1116.32,1066.46,1125.43,1234.29,1095.82,1191.21,1139.12,931.042,1254.02,1080.23,1051.15,1035.37,1340.63,1173.03,926.629,1180.19,1059.97,1173.35,1029.83,943.416,1043.1,867.078,852.112,1157.2,1189.81,1114.83,904.913,1084.47,1071.75,811.183,999.13,1185.61,1305.96,1205.43,1059.98,1089.58,1165.91,1091.92,1020.76,1204.84,1191.7,1240.49,1112.39,968.871,1118.52,883.354,888.35,1021.23,979.18,905.46};

int nx  = 50;
int ny  = 50;
double bx = 5.0/nx;
double by = 5.0/ny;
int nmch = 10;
int ntime = 300;
int  nmapdata = nx*ny*nmch*ntime;
double *mapdata;
std::vector<std::vector<double> > trvecADC;
int ddd=0;
double LiqVel025;
double LiqVel150;

// Diffusion  https://lar.bnl.gov/properties/
const double DiffT=9.5;  // cm2/s
const double DiffL=6.3;  // cm2/s

const double SigmaT=sqrt(1.0e-4 * 2.0 * DiffT);   // mm/sqrt(us);
const double SigmaL=sqrt(1.0e-4 * 2.0 * DiffL);   // mm/sqrt(us);

TVector3 dr(double alpha = 0.5)
{
    TVector3 dd;
    if (r.z() < 300.0) {
        dd.SetXYZ(0.0, 0.0, dstep);
    } else {
        dd.SetXYZ(
            -(1.0 - alpha) * (r.x() - 2.5) / 10.0 * dstep,
            -(1.0 - alpha) * (r.y() - 2.5) / 10.0 * dstep,
             dstep
        );
    }
    return dd;
}

double gexp(double t, double delta, double tau, double sigma)
{
    double lambda = 1.0 / tau;
    t -= delta;
    if (t < -sigma * 20.0) return 0;
    double val = (lambda / 2.0)
        * std::exp(lambda / 2.0 * (lambda * sigma * sigma - 2.0 * t))
        * TMath::Erfc((lambda * sigma * sigma - t)/(std::sqrt(2.0)*sigma));
    return val;
}

Double_t wfun(Double_t* x, Double_t* par)
{
    // par[0] = tau, par[1] = sigma
    return gexp(x[0], 0.0, par[0], par[1]);
}

double BirksLaw(double EField, double dEdx){
  double paraA   = 0.8;
  double parak   = 0.0486;
  return  paraA/(1+(parak/EField*dEdx));         
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

double DriftWeight(double driftTime, int mode) {
    if (mode == 0) return 1.0;
    if (mode == 1) {
        double tau = 200.0;
        return std::exp(-driftTime / tau);
    }
    return 1.0;
}


double nnn(double t, double t0, double w){
  
  return cos( 2.0 * TMath::Pi() * (t-t0)/w);
  
}

//static double gEnergyWeight = 20*1e19; // energy Weight
const double rho = 1.4; // g/cm3
const double timeconv = 0.4; //us


TChain *fChain, *fheader,*fbeam, *fpar;
TTree*tree;
TString InputRootFile, OutputRootFile;
//struct stat statBuf;
// Declaration of leaf types     
Int_t event, pid;
Double_t mass; 
Int_t mat;
Double_t t, x, y, z, e, l, k;   
Int_t  trkid;
//G4 Header
Int_t  eventID, CaptureFlag;   
Double_t injecposx, injecposy, injecposz, injecdirx, injecdiry, injecdirz;   
Double_t stopposx, stopposy, stopposz, E_in;
Int_t  Num_pi0,  Num_secpi0;
//Beam
Int_t eveID;
Double_t x_bpc2,y_bpc2,z_bpc2,dxdz_bpc2, dydz_bpc2, beam_Mom;
//Par
Int_t peventID, pProcess, pmat, ptrkid, ppid, pparid;
Double_t pMom_x, pMom_y, pMom_z, pMom, pk;
// List of branches
TBranch *b_event,*b_pid,*b_mass,*b_mat,*b_t,*b_x,*b_y,*b_z,*b_e,*b_l;
TBranch *b_k,*b_trkid;
TBranch*b_eventID,*b_CaptureFlag,*b_InjecPos,*b_StopPos,*b_NumPi0,*b_NumSecPi0;   //!
Int_t trigCounter[4]={0};
Int_t CLKCounter[4]={0};
Int_t EntryNumber=0;

Int_t Header[20],MHeaderI[20];
Double_t Edep=0, HeaderD[20],MHeaderD[20];

Int_t Process;Double_t StopPos[3]={0}, InjecPos[3]={0};

int  num_ppi =0, num_npi =0, num_p =0, num_pi0 =0, num_ee=0,
  num_secppi =0, num_secnpi =0, num_secp =0, num_secpi0 =0,
  num_oth =0, num_secoth =0;
int pentry = 0;
int pbarflag = 0;
int mode_ =0;

void LoadG4RootFile(Int_t FileNumber,const char *fname="sim"){
  delete fChain;
  fChain = new TChain("trk");
  delete fheader;
  fheader = new TChain("Header");
  delete fbeam;
  fbeam = new TChain("Beam");
  delete fpar;
  fpar = new TChain("Par");
  if(strcmp(fname, "PBAR") == 0){
    cout <<"Particle Name: " << fname << endl;
    if(mode_){
      //InputRootFile = Form("./root/Geant4/GiBUU/B4_%d.root",FileNumber);
      InputRootFile = Form("./input/B4_%d.root",FileNumber);
    }else{
      InputRootFile = Form("./input/B4_%d.root",FileNumber);
    }
  }else if(strcmp(fname, "DBAR") == 0){
    cout <<"Particle Name: " << fname << endl;
    InputRootFile = Form("./%s/B4_%d.root",fname,FileNumber);
  }
  //InputRootFile = Form("root/Geant4/T98/%s/B4_%d.root",fname,FileNumber);
  cout << " " << InputRootFile << endl;
  fChain->Add(InputRootFile);
  fheader->Add(InputRootFile);
  fbeam->Add(InputRootFile);
  fpar->Add(InputRootFile);
 
  fChain->SetBranchAddress("event", &event, &b_event);
  fChain->SetBranchAddress("pid", &pid, &b_pid); 
  fChain->SetBranchAddress("mass", &mass, &b_mass);
  fChain->SetBranchAddress("trkid", &trkid, &b_trkid);   
  fChain->SetBranchAddress("mat", &mat, &b_mat); 
  fChain->SetBranchAddress("t", &t, &b_t);  
  fChain->SetBranchAddress("x", &x, &b_x);  
  fChain->SetBranchAddress("y", &y, &b_y);  
  fChain->SetBranchAddress("z", &z, &b_z);  
  fChain->SetBranchAddress("e", &e, &b_e);  
  fChain->SetBranchAddress("l", &l, &b_l);
  
  fheader->SetBranchAddress("EventID", &eventID);
  fheader->SetBranchAddress("CaptureFlag", &CaptureFlag);
  fheader->SetBranchAddress("InjecPosx", &injecposx);
  fheader->SetBranchAddress("InjecPosy", &injecposy);
  fheader->SetBranchAddress("InjecPosz", &injecposz);
  fheader->SetBranchAddress("StopPosx", &stopposx);
  fheader->SetBranchAddress("StopPosy", &stopposy);
  fheader->SetBranchAddress("StopPosz", &stopposz);
  fheader->SetBranchAddress("InjecDirx", &injecdirx);
  fheader->SetBranchAddress("InjecDiry", &injecdiry);
  fheader->SetBranchAddress("InjecDirz", &injecdirz);
  fheader->SetBranchAddress("E_in", &E_in);

  //Beam
  fbeam->SetBranchAddress("EventID", &eveID);
  fbeam->SetBranchAddress("x_bpc2", &x_bpc2);
  fbeam->SetBranchAddress("y_bpc2", &y_bpc2);
  fbeam->SetBranchAddress("z_bpc2", &z_bpc2);
  fbeam->SetBranchAddress("dxdz_bpc2", &dxdz_bpc2);
  fbeam->SetBranchAddress("dydz_bpc2", &dydz_bpc2);
  fbeam->SetBranchAddress("beam_Mom", &beam_Mom);
  
  //paricle
  fpar->SetBranchAddress("EventID", &peventID);
  fpar->SetBranchAddress("pProcess", &pProcess);
  fpar->SetBranchAddress("pmat", &pmat);
  fpar->SetBranchAddress("ptrkid", &ptrkid);
  fpar->SetBranchAddress("ppid", &ppid);
  fpar->SetBranchAddress("pparid", &pparid);
  fpar->SetBranchAddress("pMom_x", &pMom_x);
  fpar->SetBranchAddress("pMom_y", &pMom_y);
  fpar->SetBranchAddress("pMom_z", &pMom_z);
  fpar->SetBranchAddress("pMom", &pMom);
  fpar->SetBranchAddress("pk", &pk);
}


void diff(double x0, double y0, double z0, int icx0, int icy0, double alpha=0.5,  double scut=3, int weightMode = 0){
  double t0 = (300.-z0)/LiqVel025;
  double tdiff = SigmaT*sqrt(t0);
  double xs = x0 - scut*tdiff;
  double xe = x0 + scut*tdiff;
  double ys = y0 - scut*tdiff;
  double ye = y0 + scut*tdiff;
  double pi = 3.14159265358979;
  
  double Energy = e*1e6; // energy (MeV->eV)
  double dedx = e/(l*0.1*rho); //MeV/cm
  double R = BirksLaw(0.25, dedx); // Recombination factor
  //double R = 0.68;  // Recombination	    
  double W = 23.7;  // eV/e-
  double Q_e = 1.6e-19;  // C/e-
  double gain = 8000./50.*1e15; //ADC*us/C
  double purity  = DriftWeight(t0, weightMode);
  double scale = 10e15*3.;
  //printf("Here %f %f %f %f %f %f %f %f %f\n",x0,y0,xs,xe,ys,ye,tdiff,t0,SigmaT);
  
  for (double xx=xs;xx<xe;xx=xx+bx){
    for (double yy=ys;yy<ye;yy=yy+by){
      double r2 = (xx-x0)*(xx-x0) +(yy-y0)*(yy-y0);
      if (r2>scut*scut*tdiff*tdiff) continue;
      double s = 1.0/(2*pi*tdiff*tdiff) * exp(-0.5 * r2/(tdiff*tdiff)) *bx*by * 3;
      
      int ix = xx/bx;
      int iy = yy/by;
    
      int xshift = 0;
      int yshift = 0;
      if(ix >= nx) {
	ix = ix - nx;
	xshift=1;
      }
      if(ix < 0) {
	ix = ix + nx;
	xshift=-1;
      }
      if(iy >= ny) {
	iy = iy - ny;
	yshift=1;
      }
      if(iy < 0) {
	iy = iy + ny;
	yshift=-1;
      }
      
      for (int ich=0; ich<nmch;ich++){
	int icx = icx0 + (ich-2) + xshift;
	int icy = icy0 + (ich-2-5) + yshift;
	if (ich <5){
	  if (icx < 0) continue;
	  if (icx >= 64) continue;
	}else{
	  if(icy < 0) continue;
	  if(icy >= 64) continue;
	}
	
	int imapdata = ntime*(ich + nmch*(iy + ny*ix));
	double tt0 = t0 + -5.0 + 30.*0.5/300.;
	int  ibt0 = (tt0+50.)/400.*1000.+1;
	  
	// for (int itime=0; itime<ntime;itime++){
	//   int ibt=ibt0+itime/4;
	//   //double purity  = DriftWeight(ibt, weightMode);
	//   double purity  = 1;
	//   double q = mapdata[imapdata+itime]*s * timeconv / 0.4 * purity * Energy * (0.79 / W) * Q_e * R*scale;
	//   //double q = mapdata[imapdata+itime]*s * timeconv / dt2Step * purity * Energy * (0.79 / W) * Q_e * R;
	//   if (ich<5) hx2[icx]->AddBinContent(ibt,q);
	//   else hy2[icy]->AddBinContent(ibt,q);
	//   //std::cout << "Q : " << q << std::endl;
	//}
	for (int itime=0; itime<ntime;itime++){
  int ibt=ibt0+itime/4;
  if (ibt < 1 || ibt > 1000) continue;

  //double purity  = DriftWeight(ibt, weightMode);
  double purity  = 1;
  double q = mapdata[imapdata+itime]*s * timeconv / 0.4 * purity * Energy * (0.79 / W) * Q_e * R*scale;
  //double q = mapdata[imapdata+itime]*s * timeconv / dt2Step * purity * Energy * (0.79 / W) * Q_e * R;
  if (ich<5) hx2[icx]->AddBinContent(ibt,q);
  else hy2[icy]->AddBinContent(ibt,q);
}
    }
    }
  }
  return;
}

void track2(double startX=0.,
            double startY=150.,
            double startZ=150.,
            double endX=300.,
            double endY=160.,
            double endZ=160.,
            double alpha      = 0.5,
            int    ddd_        = 0,
            int    weightMode = 0)
{
  const double cellSize = 5.0;

  // f1 = new TF1("f1", wfun, -10, 20, 2);
  // f1->SetParameters(3.0, 1.0);

  double t0x=gRandom->Uniform(0,10);
  double wwx=gRandom->Gaus(10,1);
  double t0y=gRandom->Uniform(0,10);
  double wwy=gRandom->Gaus(10,1);
    
  for (int i = 0; i < kNChan; ++i) {
    hx[i]->Reset();
    hy[i]->Reset();
    hx2[i]->Reset();
    hy2[i]->Reset();
    
    hx[i]->SetLineColor(51 + (i - 35) * 3);
    hy[i]->SetLineColor(51 + (i - 35) * 3);
    hy[i]->SetLineStyle(2);
  }
  n1  = new TNtuple("n1", "n1", "x:y:z:ox:oy:qx:qy");
  double liqVel04 = LiqVel(0.25);
  //double liqVel20 = LiqVel(2.);
  
  Long64_t nbytes = 0, nb = 0;
  Int_t nentries = fChain->GetEntries();
  
  fpar -> GetEntry(pentry);
  for (Long64_t jentry=EntryNumber; jentry<nentries;jentry++) {
    if(jentry%100==0)cout <<"Entry: " << jentry << endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(!mode_&&pid==-2212)continue;
    if( (event==0&&trigCounter[0]!=0 ) || (mat==0&&event==0&&pid==0)){
      cout << "Entry: " << jentry << endl;
      continue;
    }
    if((int)trigCounter[0]==event){
      if(mat==1 & e > 0.001){
      	Edep += e;
      	int tmpX = static_cast<int>(std::floor((x+160.) / cellSize));
      	int tmpY = static_cast<int>(std::floor((z+160.) / cellSize));
      	if (tmpX < 2 || tmpX >= kNChan-2) continue;
      	if (tmpY < 2 || tmpY >= kNChan-2) continue;
      	ichx = tmpX;
      	ichy = tmpY;
      	double locx = (x+160.) - tmpX * cellSize;
      	double locy = (z+160.) - tmpY * cellSize;
      	t0 = (300.0 - (y+150.)) / liqVel04;
      	//printf("%i %f %f %f %f\n",pid,e,l,dedx,b);
      	//printf("%f %f %i %i\n",x,z,xch,zch);
      	diff(locx, locy, y+160., ichx, ichy);
      	//diff(locx, locy, rt.Z(), ichx, ichy);
      }
    }
    if ((int)trigCounter[0]!=event || jentry==nentries-1){
      if(pbarflag==1){
	// while(peventID==trigCounter[0]){
	//   cout <<"A" << endl;
	//   if(pk>10){
	//     if(ppid==2212){num_p++;//proton
	//     }else if(ppid==211){num_ppi++;//pi+ 
	//     }else if(ppid==-211){num_npi++;// pi-
	//     }else{num_oth++;}
	//     if(pparid==1){
	//       if(ppid==2212){num_secp++;//proton
	//       }else if(ppid==211){num_secppi++;//pi+
	//       }else if(ppid==-211){num_secnpi++;// pi-
	//       }else{num_secoth++;}
	//     }
	//   }
	//   pentry++;
	//   fpar -> GetEntry(pentry);
	// }
      }
      
      Header[0] = num_ppi;//RunNumber
      Header[1] = num_npi;//FileNumber
      Header[2] = num_p;//EventNumber
      Header[3] = num_pi0;//TIGArsync
      Header[4] = num_oth;//Trigtype
      Header[5] = num_secppi;//FADCFileNum
      Header[6] = num_secnpi;//FADCEventNum
      Header[7] = num_secp;//Number of Pileup
      
      double pmass = 938.272;
      
      HeaderD[0] = z_bpc2;//X pos @BPC2
      HeaderD[1] = y_bpc2;//Y pos @BPC2
      HeaderD[2] = dxdz_bpc2;//dy/dz_bpc @ BPC2
      HeaderD[3] = dydz_bpc2;//dy/dz_bpc @ BPC2
      HeaderD[4] = std::sqrt(beam_Mom*beam_Mom+pmass*pmass)-pmass;//Momentum->Kinetic Energy
      // HeaderD[5] = num_secppi;//Reconstructed X
      // HeaderD[6] = num_secnpi;//Reconstructed Y
      // HeaderD[7] = num_secp;  //Reconstructed Z
      
      MHeaderI[0] = num_secppi;
      MHeaderI[1] = num_secnpi;
      MHeaderI[2] = num_secp;
      MHeaderI[3] = num_secpi0;
      MHeaderI[4] = num_secoth;
      MHeaderI[5] = num_ppi;
      MHeaderI[6] = num_npi;
      MHeaderI[7] = num_p;
      MHeaderI[8] = num_pi0;
      MHeaderI[9] = num_oth;
      MHeaderI[10] = CaptureFlag;
      
      MHeaderD[0] = 150.-injecposx;
      MHeaderD[1] = 150.+injecposy;
      MHeaderD[2] = 150.-injecposz;
      MHeaderD[3] = injecdirz/injecdirx;//dZ/dx_i
      MHeaderD[4] = injecdiry/injecdirx;//dZ/dy_i
      MHeaderD[5] = E_in;//E_i
      MHeaderD[6] = 150.-stopposx;
      MHeaderD[7] = 150.+stopposy;
      MHeaderD[8] = 150.-stopposz;
      MHeaderD[9] = Edep;//Totsl Edep
      
      cout << "Stop X: " << stopposx << " Y: " << stopposy <<" Z: " << stopposz << endl;
      EntryNumber = jentry;
      std::cout << "Finished jentry : " << jentry << std::endl;
      std::cout << "Trig Counter : " << (int)trigCounter[0] << " Event : " << event << std::endl;
      break;
    }
  }
}

void sim_GiBUU(int nev, int nmap=100, int file_=0,  int seed=0, const char *fname="sim",int dx=0, int dy=0, int ddd_=1,int mode =0){
  ddd = ddd_;
  gRandom->SetSeed(seed);
  mode_ = mode; // mode=0:FTF, mode=1: GiBUU
  // if(mode_){
  //   OutputRootFile= Form("./root/sim/G4_%i_sigma0_x%d_y%d_%d.%s.root",file_,dx,dy,ddd,fname);
  // }else{
  //   OutputRootFile= Form("./root/sim/G4_%i_sigma0_x%d_y%d_%d.%s.root",file_,dx,dy,ddd,fname);
  // }
  // TFile *ofile = new TFile(OutputRootFile,"recreate");
  TString outDir;
  if(mode_){
    outDir = Form("./output");
  }else{
  outDir = Form("./%s", fname);
  }
  gSystem->mkdir(outDir, true);

  OutputRootFile = Form("%s/G4_%i_sigma0_x%d_y%d_%d.%s.root",
                      outDir.Data(), file_, dx, dy, ddd, fname);

TFile *ofile = new TFile(OutputRootFile,"recreate");
 
  nx  = 100;
  ny  = 100;
  bx = 5.0/nx;
  by = 5.0/ny;
  nmch = 10;
  ntime = 300;
  nmapdata = nx*ny*nmch*ntime;
  if(strcmp(fname, "PBAR") == 0)pbarflag = 1;
  
  mapdata = (double *)malloc(sizeof(double) * nmapdata);
  
  // FILE *fdata = fopen(Form("mapdata/map_x%d_y%d_%d.dat",dx,dy,ddd),"rb");
  // //FILE *fdata = fopen(Form("mapdata/map_%i.dat",nmap),"rb");
  // cout << fread(mapdata,(sizeof *mapdata),nmapdata,fdata)<<endl;;
  // fclose(fdata);
  // printf("Map data ok\n");
 mapdata = (double *)malloc(sizeof(double) * nmapdata);

FILE *fdata = popen(Form("bzcat mapdata/map_x%d_y%d_%d.dat.bz2",dx,dy,ddd),"r");
//FILE *fdata = fopen(Form("mapdata/map_x%d_y%d_%d.dat",dx,dy,ddd),"rb");
//FILE *fdata = fopen(Form("mapdata/map_%i.dat",nmap),"rb");
if (!fdata) {
  printf("ERROR: cannot open mapdata/map_x%d_y%d_%d.dat.bz2\n", dx, dy, ddd);
  return;
}

size_t nread = fread(mapdata, sizeof(*mapdata), nmapdata, fdata);
cout << nread << endl;
fclose(fdata);

if (nread != (size_t)nmapdata) {
  printf("ERROR: mapdata size mismatch (%zu / %d)\n", nread, nmapdata);
  return;
}

printf("Map data ok\n");

 
  f1 = new TF1("f1",wfun,-10,20,2);
  f1->SetParameters(3.0,1.0);
  for (int ii =0;ii<50;ii++){
    double dt2 = -2.0 + ii*0.4;
    double s = f1->Eval(dt2);
    ss[ii]=s;
  }  
  int nch=128;
  int nhch=64;
  int nbin=1000;
  LiqVel025 = LiqVel(0.25);
  LiqVel150 = LiqVel(1.5);  
  
  std::vector<std::vector<double> > trvecADC;
  trvecADC.resize(nch*2, std::vector<double>(nbin+1));
  TTree *WaveTree = new TTree("WaveTree","WaveTree");
  WaveTree -> Branch("trigCounter", &trigCounter, "trigCounter[4]/I");
  WaveTree -> Branch("CLKCounter",  &CLKCounter,  "CLKCounter[4]/I"); 
  WaveTree -> Branch("trvecADC",  &trvecADC);
  WaveTree->Branch("Header", Header,"Header[20]/I");//0=Sec Num, 1=pi+-, 2=pi0, 3= proton, 4=gamma(in TPC)
  WaveTree->Branch("HeaderD", HeaderD,"HeaderD[20]/D");//0=Sec Num, 1=pi+-, 2=pi0, 3= proton, 4=gamma(in TPC)
  WaveTree->Branch("MHeaderI", MHeaderI,"MHeaderI[20]/I");//0=Sec Num, 1=pi+-, 2=pi0, 3= proton, 4=gamma(in TPC)
  WaveTree->Branch("MHeaderD", MHeaderD,"MHeaderD[20]/D");//0=Sec Num, 1=pi+-, 2=pi0, 3= proton, 4=gamma(in TPC)

  for(Int_t i=0; i < nhch; i++){
    hx[i] = new TH1D(Form("hx3_%d", i), Form("hx3_%d", i), 1000, -50, 350);
    hy[i] = new TH1D(Form("hy3_%d", i), Form("hy3_%d", i), 1000, -50, 350);
    hx2[i] = new TH1D(Form("hx2_%d", i), Form("hx2_%d", i), 1000, -50, 350);
    hy2[i] = new TH1D(Form("hy2_%d", i), Form("hy2_%d", i), 1000, -50, 350);
  }

  LoadG4RootFile(file_,fname);
  fheader->GetEntry(0);fbeam->GetEntry(0);
  
  for (int i=0;i<nev;i++) {
    std::cout << "\rEvent " <<i  << "/" << nev << " start" << std::flush;;
    cout << "" << endl;
    // double delta1=gRandom->Gaus(0,10);
    // double delta2=gRandom->Gaus(0,10);
    // double delta3=gRandom->Gaus(0,10);
    // double delta4=gRandom->Gaus(0,10);
    track2(0.,150.,150., 300.,150.,150.);

    for (int ich=0; ich<64;ich++){
      int nBins = hx[ich]->GetNbinsX();
      for (int b = 1; b <= nBins; ++b) {
	double noiseX = pedmean[ich];
	double noiseY = pedmean[ich+64];
	hx[ich]->Fill((Double_t)b*0.4-50.-0.2, noiseX);
	hy[ich]->Fill((Double_t)b*0.4-50.-0.2, noiseY);
	// hx[i]->AddBinContent(b, noiseX);
	// hy[i]->AddBinContent(b, noiseY);
      }
    }
    
    double pi = 3.14159265358979;  
    for (int ich=0; ich<64;ich++){
      for (int itime=0; itime<1000;itime++){
	double qx = hx2[ich]->GetBinContent(itime+1);
	double qy = hy2[ich]->GetBinContent(itime+1);
	double t0 = hx2[ich]->GetBinCenter(itime+1);
	double sigt = 0.;
	if (t0>0) sigt= SigmaL*sqrt(t0);
	for (double tt = -sigt*3;tt<sigt*3;tt=tt+0.4){
	  double s = 1.0/sqrt(2*pi)/sigt * exp(-0.5 * tt*tt/sigt/sigt)*0.4;
	  hx[ich]->Fill(t0+tt,qx*s);
	  hy[ich]->Fill(t0+tt,qy*s);
	}
      }
    }
    
    for(Int_t ich=0; ich < nhch; ich++){
      for(Int_t ibin=0; ibin < nbin; ibin++){
	trvecADC[ich*2][ibin]    = hx[ich]->GetBinContent(ibin+1);
	trvecADC[128+ich*2][ibin] = hy[ich]->GetBinContent(ibin+1);
      }
    }
    trigCounter[0]=i;
    //if(CaptureFlag==1&&Num_secpi0==0)WaveTree->Fill();
    cout << "Event: " << event << " CaptureFlag: " << CaptureFlag << endl;
    WaveTree->Fill();
    trigCounter[0]++;
    
    num_ppi=0; num_npi=0; num_p =0; num_pi0=0; num_oth=0;
    num_secppi=0;  num_secnpi=0;num_secp =0; num_secpi0=0; num_secoth=0;
    
    cout << "TrigCounter: " << trigCounter[0] << " Event: " << event << endl;
    //fheader->GetEntry(trigCounter[0]);
    fheader->GetEntry(event); fbeam->GetEntry(event);
    cout << "Event: " << event << " CaptureFlag: " << CaptureFlag << endl;
  }
  
  WaveTree->Print();
  ofile->cd();    
  WaveTree->Write();
  ofile->Close();

  return;
  
}
