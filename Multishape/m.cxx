/****************************************************************
Readme!!
This code takes the total multiplicity (n), no. of pairs exhibiting femtoscopy (a) and no. of pairs with elliptical flow (b) as arguments 
and simulates an ideal coleration function. Ideally the pairs generate from femtoscopy are very close to each other and 
so the relative rapidity and relative azimuthal angle is nearly zero - resulting a peak near zero in the correlation funstion R2.
The various patterns of the anisotropic flow of ultra-RHIC can be characterized by the Fourier expansion of the invariant 
triple distribution of particle pairs. See eqn. 2 and 3 of the reference below:
https://arxiv.org/abs/1102.3010 
The conservation of momentum (represented by the first term i.e. v1 in the Fourier expansion and is sinusoidal in phi - psi i.e.cos(phi - psi)
where phi is the azimuthal angle and psi is the reaction plane angle.
Similarly the eliptical flow v2 due to the Lorentz boost is quantified by v2 i.e. cos2(phi - psi) which is evident in the ridge of 
the correlation function R2.

NOTE: The total number of particles N < 500. (N = S +2*a +2*b) where S is the # of single particles.

Code Usage:
make clean
make
./m -n 40 -a 10 -b 1

Output: Running like the command above creates a pdf called 'm_40_10_1_0.100.pdf' in the ps directory of the current working directory
******************************************************************/
#include "MultCum.h"
//
#include <TROOT.h>
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFunc.h"
#include "Math/Math.h"
#include "Math/SpecFunc.h"
#include "Math/SpecFuncMathCore.h"
#include <cmath>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstring>
//
#include <TRatioPlot.h>
#include <TSystem.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TRandom.h>
#include <TLatex.h>

double YL,YU;	// global scope for just these so i can dy-weight in the C2 fitting function
void DisplayHelp();

//---- prototypes...
void m(int N_target, int Pa_K1, int Pb_K1);
double binomial_pdf(unsigned int k, double p, unsigned int n);
double negative_binomial_pdf(unsigned int k, double p, double n);
double poisson_pdf(unsigned int n, double mu);
Double_t fBD(Double_t *x, Double_t *par);
Double_t fPoi(Double_t *x, Double_t *par);
Double_t fNBD(Double_t *x, Double_t *par);
Double_t fC2fit(Double_t *x, Double_t *par);
//
double xlimits[2];
double ylimits[2];
void RangeFinderH(int, TH1D*);
void RangeFinderG(int, TGraph*);
void RangeFinderGE(int, TGraphErrors*);
void RangeCheck();
//void XRangeFinder2D(TH2* h=0);
//void RangeFinder2D(TH2* h=0);
void XRangeFinder(TH1* h=0);
void XRangeCheck();
void GenericProject(TH2D* h2,TString steerav,TString steerax,TH1D* h1);

using namespace std;

//-------------------------------------------------------------
int main(int argc, char **argv){
  int N_target;// = 40.0;         // desired ensemble-average total multiplicity
  int Pa_K1;// =  1.0;         // avg number of pairs per evt from correlation source "a" (FEMTOSCOPY or HBT effect), PDF is Poisson
  int Pb_K1;//  =  10.0;        // avg number of pairs per evt from correlation source "b" (ELLIPTICAL FLOW), PDF is Poisson

  for (int i=1; i<argc; i++) {
    if (argc>1){
      if (argv[i][0] == '-') {
        switch (argv[i][1]) {
        case 'n':
          N_target         = atoi(argv[++i]);
          break;
        case 'a':
          Pa_K1          = atoi(argv[++i]);
          break;
	case 'b':
          Pb_K1          = atoi(argv[++i]);
          break;
	}
      }
    }
    else DisplayHelp(); 
  }
  m(N_target,Pa_K1,Pb_K1);
}
  

void DisplayHelp(){
  cout<<endl;
  cout<<"USAGE: ./m [OPTIONS]"<<endl;
  cout<<"                -h          ... Display this help"<<endl;
  cout<<"                -n  ... Total multiplicity"<<endl;
  cout<<"                -a  ... Average no. of pairs from FEMTOSCOPY or HBT effect < n"<<endl;
  cout<<"                -b  ... Average no. of pairs from ELLIPTICAL FLOW < n"<<endl;
}
//-------------------------------------------------------------
 void m(int N_target, int Pa_K1, int Pb_K1){

  //--------------------------------------
  //--------------------------------------
  //
  int NEVT	= 100000;
  double Pa_sigdy		=  0.20;	// sigma of 2D Gaussian in dy direction
  double Pa_sigdphi	=  10.0;	// sigma of 2D Gaussian in dphi direction
  //
  //--------------------------------------
  //--------------------------------------
  cout<<"Multiplicity = "<<N_target<<endl;
  int S_K1			= N_target - 2*Pa_K1 - 2*Pb_K1;		// avg number of singles per event, PDF is Poisson
  cout<<"1:Multiplicity = "<<N_target<<endl;
  if (S_K1<1){ cout<<"not enough singles = ..."<<S_K1<<endl; exit(0); }

  //double lam		= (2.*Pb_K1)/(N_target+2.*Pb_K1);
  //double v2check	= 0.5*pow(lam,2.);
  //cout<<"lambda = "<<lam<<" v2check = "<<v2check<<endl;

  if (S_K1<0){ cout<<"S_K1<0... are you sure?"<<endl; exit(0); }
  double nevtM		= (double)NEVT/1000000;
  TString dir			= TString("~/RHIC/multshape/ps/");
  TString base		= TString(Form("m_%d_%d_%d_%.3f",N_target,Pa_K1,Pb_K1,nevtM));
  cout<<"base = "<<base.Data()<<endl;
  TString outfile		= dir + base + TString(".ps");
  TString outfileO	= dir + base + TString(".ps(");
  TString outfileC	= dir + base + TString(".ps]");
  TString outfileP	= dir + base + TString(".pdf");

  //---- define poisson parent for SINGLES and PAIRS sources
  TF1* fparent_S;		// parent Poisson PDF for singles multiplicity per event
  fparent_S	= new TF1("fparent_S",fPoi,0.,100.,2);
  fparent_S	->SetParameter(0,NEVT );
  fparent_S	->SetParameter(1,S_K1 );
  TF1* fparent_Pa;		// parent Poisson PDF for pairs multiplicity per event, correlation source "a"
  fparent_Pa	= new TF1("fparent_Pa",fPoi,0.,100.,2);
  fparent_Pa	->SetParameter(0,NEVT );
  fparent_Pa	->SetParameter(1,Pa_K1);
  TF1* fparent_Pb;		// parent Poisson PDF for pairs multiplicity per event, correlation source "b"
  fparent_Pb	= new TF1("fparent_Pb",fPoi,0.,100.,2);
  fparent_Pb	->SetParameter(0,NEVT );
  fparent_Pb	->SetParameter(1,Pb_K1);
  //
  //---- there are used to describe the kinematic correlation (in dy,dphi) amongst pair particles...
  TF2* fparent_Pa_src;		// parent kinematic shape of sharp SRC correlated pair peak 
  fparent_Pa_src	= new TF2("fparent_Pa_src","[0]*TMath::Gaus(x,[1],[2])*TMath::Gaus(y,[3],[4])",-1,1,-90,90);
  fparent_Pa_src	->SetParameter(0,  10000.   );	// arbitrary
  fparent_Pa_src	->SetParameter(1,      0.   );	// SRC dy mean
  fparent_Pa_src	->SetParameter(2, Pa_sigdy  );	// SRC dy sig
  fparent_Pa_src	->SetParameter(3,      0.   );	// SRC dphi mean
  fparent_Pa_src	->SetParameter(4, Pa_sigdphi);	// SRC dphi sig
  TF1* fparent_Pb_v2;			// parent kinematic shape of ReactionPlane correlated pairs
  fparent_Pb_v2	= new TF1("fparent_Pb_v2","1.0+TMath::Cos(2.*x/57.2958)",0.,360.);
  //		fparent_Pb_v2	->SetParameter(0, 10000. );		// arbitrary
  //		fparent_Pb_v2	->SetParameter(2, Pb_v2  );		// v2 dy sig
  //
  //----- histograms for generated multiplicities per event...
  TH1D *hm0	= new TH1D("hm0","Multiplicity/Event",100,-0.5,99.5);
  hm0->SetMarkerStyle(20);
  hm0->SetMarkerSize(0.8);
  hm0->SetMarkerColor(1);
  hm0->SetLineColor(1);
  hm0->SetLineWidth(2);
  TH1D *hS	= new TH1D("hS","S Multiplicity/Event",100,-0.5,99.5);
  TH1D *hPa	= new TH1D("hPa","Pa Multiplicity/Event",100,-0.5,99.5);
  TH1D *hPb	= new TH1D("hPb","Pb Multiplicity/Event",100,-0.5,99.5);
  TH1D *hN	= new TH1D("hN","N Multiplicity/Event",100,-0.5,99.5);
  hS->SetMarkerStyle(20);
  hPa->SetMarkerStyle(20);
  hPb->SetMarkerStyle(20);
  hN->SetMarkerStyle(20);
  hS->SetMarkerSize(0.7);
  hPa->SetMarkerSize(0.7);
  hPb->SetMarkerSize(0.7);
  hN->SetMarkerSize(0.7);
  //
  TH1D *hv2dphi	= new TH1D("hv2dphi","v2 dphi",90,-90,270);
  //
  double YNB		=   20;
  YL		= -1.0;			// global scope
  YU		=  1.0;			// global scope
  double YBW		= (YU-YL)/YNB;
  double DYL		= YL-YU;			// this is the binning rho2 needs to be hermetic...	
  double DYU		= YU-YL;			// this is the binning rho2 needs to be hermetic...
  double DYNB		= (DYU-DYL)/YBW;	// so Y and DY have the same bin width... 
  cout<<" Y binning: "<<YNB<<" "<<YL<<" "<<YU<<endl;
  cout<<"DY binning: "<<DYNB<<" "<<DYL<<" "<<DYU<<" \t "<<YBW<<" \t "<<(DYU-DYL)/YBW<<endl;
  double PHINB	=   24;
  double PHIL		=    0;
  double PHIU		=  360;
  double PHIBW	= (PHIU-PHIL)/PHINB;
  double DPHINB	=   24;
  double DPHIL	=  -90;
  double DPHIU	=  270;
  //------------------------
  int		rho1ResFacY		= 1;
  int		rho1ResFacPhi	= 1;
  //------------------------
  TH1D *hmult	= new TH1D("hmult","hmult"  	   ,100,-0.5,99.5);
  TH2D *hrho1	= new TH2D("hrho1","#rho_{1}(y,#varphi);y;#varphi"    ,rho1ResFacY*YNB,YL,YU,rho1ResFacPhi*PHINB,PHIL,PHIU);
  TH2D *hrho2	= new TH2D("hrho2","#rho_{2}(dy,d#varphi);dy;d#varphi",DYNB,DYL,DYU,DPHINB,DPHIL,DPHIU);
	
  //---- Event Loop
  //
  cout<<"starting event loop..."<<endl;
  for (int ievt=0;ievt<NEVT;ievt++){
    //
    //---- GENERATE EVENT
    //---- generate number of singles and pairs in this event...
    int S,Pa,Pb,N;
    S	= fparent_S ->GetRandom();	// # singles this event
    Pa	= fparent_Pa->GetRandom();	// # pairs this event, correlation source "a"
    Pb	= fparent_Pb->GetRandom();	// # pairs this event, correlation source "b"
    if (Pa_K1==0.&&Pa!=0.) Pa=0.;	// protect (poisson technically not defined if mean is zero)
    if (Pb_K1==0.&&Pb!=0.) Pb=0.;	// protect (poisson technically not defined if mean is zero)
    N	= S + 2*Pa + 2*Pb;
    hm0->Fill(N);
    if (N<2) continue; 		// should monitor leakage here
    if (N>=500){ cout<<"array"<<endl; exit(0); } 
    hS	->Fill(S);
    hPa	->Fill(Pa);
    hPb	->Fill(Pb);
    hN	->Fill(N);
    //---- generate the kinematics of each track, including those in correlated pairs
    int k			=  0;
    double y[200]	= {0};
    double phi[200]	= {0};
    for (int j=0;j<S;j++){			//----- SINGLES
      y[k]	=   YL +     (YU-YL)*gRandom->Rndm();	// [YL,YU]
      phi[k]	= PHIL + (PHIU-PHIL)*gRandom->Rndm();	// [0,360]
      if (y[k]>=YU){ cout<<"yu... "<<y[k]<<endl; }
      if (y[k]< YL){ cout<<"yl... "<<y[k]<<endl; }
      if (phi[k]>=PHIU){ cout<<"fu... "<<phi[k]<<endl; }
      if (phi[k]< PHIL){ cout<<"fl... "<<phi[k]<<endl; }
      ++k;
    }
    for (int j=0;j<Pa;j++){			//----- PAIRS, correlation source "a" (FEMTO)
      y[k]	=   YL +     (YU-YL)*gRandom->Rndm();	// [YL,YU]
      phi[k]	= PHIL + (PHIU-PHIL)*gRandom->Rndm();	// [0,360]
      ++k;
      double srcdevdy,srcdevdphi;
      bool finding	= true;
      double ya		= y[k-1];
      double phia		= phi[k-1];
      while (finding){
	fparent_Pa_src->GetRandom2(srcdevdy,srcdevdphi);
	if ( ( (ya  +srcdevdy  )>=YL   && (ya  +srcdevdy  )<YU  )
	     && ( (phia+srcdevdphi)>=PHIL && (phia+srcdevdphi)<PHIU ) ){
	  y[k]	= y[k-1]   + srcdevdy;		// 2nd particle in pair placed w.r.t. 1st particle in pair!
	  phi[k]	= phi[k-1] + srcdevdphi;	// 2nd particle in pair placed w.r.t. 1st particle in pair!
	  ++k;
	  finding	= false;
	}
      }
    }
    for (int j=0;j<Pb;j++){			//----- PAIRS, correlation source "b" (ELLIPTIC FLOW)
      double phiRP	= 360.*gRandom->Rndm();
      y[k]	=    YL + (YU-YL)*gRandom->Rndm();		// [YL,YU]
      phi[k]	= phiRP + fparent_Pb_v2->GetRandom(); 	// [0,360]
      if (phi[k]<    0) phi[k]+=360.;
      if (phi[k]>= 360) phi[k]-=360.;
      ++k;
      y[k]	=    YL + (YU-YL)*gRandom->Rndm();		// [YL,YU]
      phi[k]	= phiRP + fparent_Pb_v2->GetRandom(); 	// [0,360]
      if (phi[k]<    0) phi[k]+=360.;
      if (phi[k]>= 360) phi[k]-=360.;
      ++k;
      double v2dphi	= (phi[k-2]-phi[k-1]);
      if (v2dphi<  -90) v2dphi += 360.; 
      if (v2dphi>= 270) v2dphi -= 360.; 
      hv2dphi->Fill(v2dphi);
    }
    if (k!=N){ cout<<"k="<<k<<" N="<<N<<" mismatch!"<<endl; exit(0); }
    //---- event is now generated!
    //
    //---- CORRELATION INCREMENT
    //---- now increment the correlation densities ...
    hmult		->Fill(N);
    for (int i=0;i<N;i++){
      hrho1	->Fill(y[i],phi[i],1.);
    }
    for (int i=0;i<N;i++){		// note with my new method, there is no need to bin center!!
      for (int j=0;j<N;j++){
	if (j!=i){
	  double dy	=   y[i] -   y[j];
	  double dphi	= phi[i] - phi[j];
	  if (dphi< DPHIL)dphi += 360.;
	  if (dphi>=DPHIU)dphi -= 360.;
	  hrho2	->Fill(dy,dphi,1.);
	}
      }
    }
    // 
  }
  cout<<"event loop done..."<<endl;
  //
  //---- end Event Loop...

  //---- scale the generator monitoring distributions so as to be "per event"
  double nseen	= hmult->GetEntries();
  cout<<"Nseen = "<<nseen<<endl;
  hv2dphi	->Scale(1./nseen);
  hS		->Scale(1./nseen);
  hPa		->Scale(1./nseen);
  hPb		->Scale(1./nseen);
  hN		->Scale(1./nseen);

  //---- calculate the correlation function...
  //
  //---- normalize densities...
  hrho1	->Scale(1./nseen);
  hrho2	->Scale(1./nseen);
  //
  //---- convolution BINGAP ODD/EVEN histogram setup!
  //----		!!!!! note the use of YNB here, this is key to this approach !!!!!
  //----		we have two overlapping&offset dy-histograms covering
  int nby	= hrho1->GetNbinsX();
  if (nby!=YNB){ cout<<"binning"<<endl; exit(0); }
  TH2D *hrho1rho1U  = nullptr;		// "upper" binning, covers whole space... nbins=YNB
  TH2D *hrho1rho1Ue = nullptr;		// "upper" binning, covers whole space... nbins=YNB
  TH2D *hrho1rho1Un = nullptr;		// "upper" binning, covers whole space... nbins=YNB
  TH2D *hrho1rho1L = nullptr;		// "lower" binning, covers all but YBW/2 at extremes... nbins=YNB-1
  TH2D *hrho1rho1Le = nullptr;		// "lower" binning, covers all but YBW/2 at extremes... nbins=YNB-1
  TH2D *hrho1rho1Ln = nullptr;		// "lower" binning, covers all but YBW/2 at extremes... nbins=YNB-1
  bool YNBISEVEN	= false;
  if ((int)nby%2==0){				// rho1 rapidity axis has an EVEN number of bins
    YNBISEVEN	= true;			// ODD bingaps cover the whole y-space
    hrho1rho1U	= new TH2D("hrho1rho1U" ,"Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1Ue	= new TH2D("hrho1rho1Ue","Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1Un	= new TH2D("hrho1rho1Un","Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1L	= new TH2D("hrho1rho1L" ,"Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
    hrho1rho1Le	= new TH2D("hrho1rho1Le","Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
    hrho1rho1Ln	= new TH2D("hrho1rho1Ln","Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
  } else if ((int)nby%2==1){		// rho1 rapidity axis has an ODD number of bins
    YNBISEVEN	= false;		// EVEN bingaps cover the whole y-space
    hrho1rho1U	= new TH2D("hrho1rho1U" ,"Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1Ue	= new TH2D("hrho1rho1Ue","Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1Un	= new TH2D("hrho1rho1Un","Upper #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB  ,DYL       ,DYU       ,DPHINB,DPHIL,DPHIU);
    hrho1rho1L	= new TH2D("hrho1rho1L" ,"Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
    hrho1rho1Le	= new TH2D("hrho1rho1Le","Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
    hrho1rho1Ln	= new TH2D("hrho1rho1Ln","Lower #rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi",YNB-1,DYL+YBW/2.,DYU-YBW/2.,DPHINB,DPHIL,DPHIU);
  }								// end rho1 rapidity axis nbins check
  //
  //---- form the denominator rho1rho1 with convolution...
  //---- 	with the upper and lower hists, there is no longer a need to "bin-center" in Increment! 
  TAxis* axy		= (TAxis*)hrho1->GetXaxis();
  TAxis* axphi	= (TAxis*)hrho1->GetYaxis();
  TAxis* axdy		= (TAxis*)hrho2->GetXaxis();
  TAxis* axdphi	= (TAxis*)hrho2->GetYaxis();
  TAxis* axdyU	= (TAxis*)hrho1rho1U->GetXaxis();
  TAxis* axdyL	= (TAxis*)hrho1rho1L->GetXaxis();
  for (int iy1=1;iy1<=axy->GetNbins();iy1++){
    double y1	= axy->GetBinCenter(iy1);
    for (int iy2=1;iy2<=axy->GetNbins();iy2++){
      double y2	= axy->GetBinCenter(iy2);
      for (int iphi1=1;iphi1<=axphi->GetNbins();iphi1++){
	double phi1	= axphi->GetBinCenter(iphi1);
	for (int iphi2=1;iphi2<=axphi->GetNbins();iphi2++){
	  double phi2	= axphi->GetBinCenter(iphi2);
	  //
	  int kBinGap_dy		= iy1 - iy2;
	  bool BINGAPISEVEN	= false;					// bin gap is odd
	  if (kBinGap_dy%2==0){BINGAPISEVEN = true;  }	// bin gap is even
	  //
	  double n1	= hrho1	->GetBinContent(iy1,iphi1);
	  double n1e	= hrho1	->GetBinError  (iy1,iphi1);
	  double n2	= hrho1	->GetBinContent(iy2,iphi2);
	  double n2e	= hrho1	->GetBinError  (iy2,iphi2);
	  //
	  double dy	=   y1 -   y2;
	  double dphi	= phi1 - phi2;
	  if (dphi< DPHIL)dphi += 360.;
	  if (dphi>=DPHIU)dphi -= 360.;
	  //int kb; double dyc,dphic;
	  if ( ( YNBISEVEN && !BINGAPISEVEN)
	       || (!YNBISEVEN &&  BINGAPISEVEN)) {
	    hrho1rho1U	->Fill(dy,dphi,  n1*n2 );	
	    hrho1rho1Ue	->Fill(dy,dphi, (n1e*n2)*(n1e*n2)+(n2e*n1)*(n2e*n1) );	
	    hrho1rho1Un	->Fill(dy,dphi,  1.0   );	
	  } else {
	    hrho1rho1L	->Fill(dy,dphi,  n1*n2 );	
	    hrho1rho1Le	->Fill(dy,dphi, (n1e*n2)*(n1e*n2)+(n2e*n1)*(n2e*n1) );	
	    hrho1rho1Ln	->Fill(dy,dphi,  1.0   );	
	  }
	  //
	}
      }
    }
  }
  //
  //---- now form C2,R2....
  TH2D* hrho1rho1		= (TH2D*)hrho2->Clone("hrho1rho1");
  hrho1rho1		->SetTitle("#rho_{1}#rho_{1}(dy,d#varphi);dy;d#varphi");
  hrho1rho1		->Reset();
  TH2D* hrho1rho1e	= (TH2D*)hrho1rho1->Clone("hrho1rho1e");
  hrho1rho1e	->Reset();
  TH2D* hC2			= (TH2D*)hrho2->Clone("hC2");
  hC2			->SetTitle("C_{2}(dy,d#varphi);dy;d#varphi");
  hC2			->Reset();
  TH2D* hC2tri		= (TH2D*)hrho2->Clone("hC2");		// C2 divided by the dy "triangle"
  hC2tri		->SetTitle("C_{2}(dy,d#varphi)/#Delta;dy;d#varphi");
  hC2tri		->Reset();
  TH2D* hR2			= (TH2D*)hrho2->Clone("hR2");
  hR2			->SetTitle("R_{2}(dy,d#varphi);dy;d#varphi");
  hR2			->Reset();
  for (int idy=1;idy<=axdy->GetNbins();idy++){				// loop over rho2 binning!
    for (int idphi=1;idphi<=axdphi->GetNbins();idphi++){	// loop over rho2 binning!
      double vdy		=  axdy  ->GetBinCenter (idy  );
      double vdphi	=  axdphi->GetBinCenter (idphi);
      double vrho2	=  hrho2 ->GetBinContent(idy,idphi);
      double vrho2e	=  hrho2 ->GetBinError  (idy,idphi);
      int kbU			=  axdyU ->FindBin(vdy);
      int kbL			=  axdyL ->FindBin(vdy);
      double vrho11,vrho11e,vrho11n,vrho11U,vrho11Ue,vrho11Un,vrho11L,vrho11Le,vrho11Ln;
      //if (idphi==1)cout<<idy<<" "<<vdy<<" "<<kbU<<" "<<kbL<<" \t "
      //				<<axdy->GetNbins()<<" "<<axdy->GetXmin()<<" "<<axdy->GetXmax()<<" "
      //				<<axdyU->GetNbins()<<" "<<axdyU->GetXmin()<<" "<<axdyU->GetXmax()<<" "
      //				<<axdyL->GetNbins()<<" "<<axdyL->GetXmin()<<" "<<axdyL->GetXmax()<<" "
      //				<<endl;
      if (idy==1||idy==axdy->GetNbins()){		// dy is in-range only in the UPPER histogram
	vrho11		=      hrho1rho1U ->GetBinContent(kbU,idphi)/2.;				// half a bin here!
	vrho11n		=      hrho1rho1Un->GetBinContent(kbU,idphi);					//
	vrho11e		= sqrt(hrho1rho1Ue->GetBinContent(kbU,idphi))/vrho11n*sqrt(2.);	// half a bin here!
      } else {								// dy is in-range in both UPPER and LOWER histograms
	vrho11U		=      hrho1rho1U ->GetBinContent(kbU,idphi);
	vrho11Un	=      hrho1rho1Un->GetBinContent(kbU,idphi);
	vrho11Ue	= sqrt(hrho1rho1Ue->GetBinContent(kbU,idphi))/vrho11Un;
	vrho11L		=      hrho1rho1L ->GetBinContent(kbL,idphi);
	vrho11Ln	=      hrho1rho1Ln->GetBinContent(kbL,idphi);
	vrho11Le	= sqrt(hrho1rho1Le->GetBinContent(kbL,idphi))/vrho11Ln;
	vrho11		= ( vrho11U + vrho11L )/2.;	
	vrho11e		= sqrt(vrho11Ue*vrho11Ue + vrho11Le*vrho11Le);
      }
      double vC2=0.,vC2e=0.,vR2=0.,vR2e=0.;
      if (vrho11>0.){
	vC2		 = vrho2 - vrho11;
	vC2e	 = sqrt(vrho2e*vrho2e + vrho11e*vrho11e);
	vR2		 = vrho2 / vrho11; 
	vR2e	 = vR2*sqrt( pow(vrho2e/vrho2,2) + pow(vrho11e/vrho11,2) );
	vR2		-= 1.0; 	// this -1 has no uncertainty (known exactly)
      }
      double valdytri	= (1.0 - fabs(vdy)/(YU-YL));		// dy-triangle height (will be used as a C2 scale factor).
      hrho1rho1	->SetBinContent(idy,idphi,vrho11 );
      hrho1rho1	->SetBinError  (idy,idphi,vrho11e);
      hC2			->SetBinContent(idy,idphi,vC2    );
      hC2			->SetBinError  (idy,idphi,vC2e   );
      hC2tri		->SetBinContent(idy,idphi,vC2 /valdytri);	// C2 scaled by the dy-triangle
      hC2tri		->SetBinError  (idy,idphi,vC2e/valdytri);	// C2 scaled by the dy-triangle
      hR2			->SetBinContent(idy,idphi,vR2    );
      hR2			->SetBinError  (idy,idphi,vR2e   );
    }
  }
  //
  //---- project R2 onto dy and dphi...
  TH1D *hR2dy		= new TH1D("hR2dy"  ,"R2dy"  ,DYNB,DYL,DYU);
  TH1D *hR2dphi	= new TH1D("hR2dphi","R2dphi",DPHINB,DPHIL,DPHIU);
  hR2dy->SetLineWidth(2); 
  hR2dy->SetMarkerStyle(1); 
  hR2dy->SetMarkerColor(1); 
  hR2dy->SetLineColor(1); 
  hR2dphi->SetLineWidth(2); 
  hR2dphi->SetMarkerStyle(1); 
  hR2dphi->SetMarkerColor(1); 
  hR2dphi->SetLineColor(1); 
  GenericProject(hR2,"Avg","X",hR2dy  );
  GenericProject(hR2,"Avg","Y",hR2dphi);
  //		
  //for (int idy=1;idy<axdy->GetNbins();idy++){
  //	for (int idphi=1;idphi<axphi->GetNbins();idphi++){
  //		cout<<axdy->GetBinCenter(idy)<<" "
  //			<<axdphi->GetBinCenter(idphi)<<" \t "
  //			<<hrho1rho1->GetBinContent(idy,idphi)<<" "
  //			<<hrho2->GetBinContent(idy,idphi)<<" "
  //			<<hR2->GetBinContent(idy,idphi)<<" "
  //			<<endl;
  //	}
  //}
  //
  //---- done calculating the correlation function...

  //---- here we'll look at the multiplicity distribution itself. 
  //
  //---- get the cumulants of the mult distribution...
  const int NCUM	= 12;
  static const char* cumnames[NCUM]	= {"K_{1}","K_{2}","K_{3}","K_{4}",
    "K_{2}/K_{1}","K_{3}/K_{1}",
    "K_{3}/K_{2}","K_{4}/K_{2}","K_{4}/K_{3}",
    "f1","f2","f3"};
  double K[NCUM]	= {0};
  MultCum* MC		= new MultCum(hm0);
  for (int iordm1=0;iordm1<NCUM;iordm1++){
    int iord		= iordm1 + 1;
    K[iordm1]		= MC->GetCumulant(iord);
  }
  cout<<"K1="<<K[0]<<" f1="<<K[9]<<" f2="<<K[10]<<" f3="<<K[11]<<endl;
  double f1=K[9],f2=K[10],f3=K[11];
  //---- done with cumulants...
  //
  //---- now lets fit the mult distribution with some basic forms (poisson, NBD, ...)
  TH1D* hm0t		= (TH1D*)hm0->Clone("hm0t");
  double rawK1	= hm0t->GetMean();
  double rawK2	= pow(hm0t->GetRMS(),2);
  //
  TF1* ffit_BD		= new TF1("ffit_BD"     ,fBD ,0.,100.,3);
  TF1* ffit_Poisson	= new TF1("ffit_Poisson",fPoi,0.,100.,2);
  TF1* ffit_NBD		= new TF1("ffit_NBD"    ,fNBD,0.,100.,3);
  TH1D *hfitpoi	= new TH1D("hfitpoi","hfitpoi",100,-0.5,99.5);
  TH1D *hfitnbd	= new TH1D("hfitnbd","hfitnbd",100,-0.5,99.5);
  TH1D *hratnbd	= new TH1D("hratnbd","hratnbd",100,-0.5,99.5);
  //
  double nbd_p	= 1.-(rawK2-rawK1)/rawK2;
  double nbd_r	= rawK1*rawK1/(rawK2-rawK1);
  ffit_NBD->SetParameter(0,NEVT);
  ffit_NBD->SetParameter(1,nbd_p);
  ffit_NBD->SetParameter(2,nbd_r);
  hm0t->Fit("ffit_NBD","Q0R"); 
  double chi2_NBD			= ffit_NBD->GetChisquare();
  double ndof_NBD			= ffit_NBD->GetNDF();
  double redchi2_NBD		= chi2_NBD/ndof_NBD;
  //
  ffit_Poisson->SetParameter(0,NEVT);
  ffit_Poisson->SetParameter(1,rawK1);
  hm0t->Fit("ffit_Poisson","Q0R"); 
  double chi2_Poisson		= ffit_Poisson->GetChisquare();
  double ndof_Poisson		= ffit_Poisson->GetNDF();
  double redchi2_Poisson	= chi2_Poisson/ndof_Poisson;
  //
  for (int ib=1;ib<hfitpoi->GetNbinsX();ib++){
    double x	= hfitpoi->GetBinCenter(ib);
    double ndat	= hm0t->GetBinContent(ib);
    double ndate= hm0t->GetBinError(ib);
    double npoi	= ffit_Poisson->Eval(x);
    double nnbd	= ffit_NBD->Eval(x);
    //double ne = 0.01;
    //hfit->SetBinError  (ib,ne);
    hfitpoi->SetBinContent(ib,npoi);
    hfitnbd->SetBinContent(ib,nnbd);
    if (nnbd>=0.1){
      hratnbd->SetBinContent(ib,ndat/nnbd);
      hratnbd->SetBinError(ib,ndate/nnbd);
    }
  }
  double nbd_p_fit	= ffit_NBD->GetParameter(1);
  double nbd_r_fit	= ffit_NBD->GetParameter(2);
  double nbd_q_fit	= 1.0 - nbd_p_fit;
  double fitK1		= nbd_r_fit*nbd_q_fit/nbd_p_fit;
  double fitK2		= fitK1/nbd_p_fit;
  //
  hfitpoi->SetEntries(NEVT);
  hfitpoi->SetLineWidth(2); 
  hfitpoi->SetLineColor(4); 
  hfitnbd->SetEntries(NEVT);
  hfitnbd->SetLineWidth(2); 
  hfitnbd->SetLineColor(kGreen+2); 
  hratnbd->SetEntries(NEVT);
  hratnbd->SetLineWidth(2); 
  hratnbd->SetLineColor(kGreen+2); 
  TCanvas *hdummy	= new TCanvas();
  TRatioPlot* rp = new TRatioPlot(hm0,hfitpoi);
  //---- done with mult distribution....

  //---- correlator integration...
  double ainte_rho1;
  double ainte_rho2;
  double ainte_rho11U;
  double ainte_rho11L;
  double ainte_rho11;
  double ainte_C2;
  double ainte_R2;
  double aint_rho1	= hrho1->IntegralAndError(1,0,1,0,ainte_rho1,"");
  double aint_rho2	= hrho2->IntegralAndError(1,0,1,0,ainte_rho2,"");
  double aint_rho11U	= hrho1rho1U->IntegralAndError(1,0,1,0,ainte_rho11U,"");
  double aint_rho11L	= hrho1rho1L->IntegralAndError(1,0,1,0,ainte_rho11L,"");
  double aint_rho11	= aint_rho11U + aint_rho11L;
  ainte_rho11	= sqrt(ainte_rho11U*ainte_rho11U + ainte_rho11L*ainte_rho11L);
  double aint_C2		= hC2->IntegralAndError(1,0,1,0,ainte_C2,"");
  double aint_R2		= hR2->IntegralAndError(1,0,1,0,ainte_R2,"");
  cout<<"Integral           rho1 = "<<aint_rho1<<" +- "<<ainte_rho1<<endl;
  cout<<"Integral           rho2 = "<<aint_rho2<<" +- "<<ainte_rho2<<endl;
  cout<<"Integral      rho1rho1U = "<<aint_rho11U<<" +- "<<ainte_rho11U<<endl;
  cout<<"Integral      rho1rho1L = "<<aint_rho11L<<" +- "<<ainte_rho11L<<endl;
  cout<<"Integral       rho1rho1 = "<<aint_rho11<<" +- "<<ainte_rho11<<endl;
  cout<<"Integral             C2 = "<<aint_C2<<" +- "<<ainte_C2<<endl;
  cout<<"Integral             R2 = "<<aint_R2<<" +- "<<ainte_R2<<endl;
  //cout<<"Integral        rho1/f1 = "<<aint_rho1/f1<<endl;
  //cout<<"Integral        rho2/f2 = "<<aint_rho2/f2<<endl;
  //cout<<"Integral rho1rho1/f1/f1 = "<<aint_rho11/f1/f1<<endl;  

  //---- now let's fit the correlator!
  //---- we're trying to extract the correlation strengths from the data to compare to generator pars
  TF2* fC2	= new TF2("fC2",fC2fit,DYL,DYU,DPHIL,DPHIU,5);
  fC2		->SetParameter(0,0.);			// global offset
  fC2		->SetParameter(1,Pa_K1/5.);		// A prefactor
  fC2		->SetParameter(2,Pa_sigdy);		// A sx
  fC2		->SetParameter(3,Pa_sigdphi);	// A sy
  fC2		->SetParameter(4,0.01);			// B prefactor
  //
  if (Pa_K1==0.){ 
    fC2->FixParameter(1,0.);	// clamp to zero
    fC2->FixParameter(2,0.);
    fC2->FixParameter(3,0.);
  }
  if (Pb_K1==0.){
    fC2->FixParameter(4,0.);
  }
  //
  TFitResultPtr fitr	= hC2tri->Fit("fC2","NRMSQ");
  double fitC2_par[5]	= {0};
  double fitC2_pare[5]= {0};
  const char* fitC2_parname[5] = {"base","A scale","A #sigma(dy)","A #sigma(d#varphi)","B scale"};
  for (int i=0;i<5;i++){
    fitC2_par[i]	= fitr->Value(i);
    fitC2_pare[i]	= fitr->FitResult::Error(i);
    cout<<"FitC2: "<<fitC2_parname[i]<<"\t = "<<fitC2_par[i]<<" +- "<<fitC2_pare[i]<<endl;
  }
  //
  int status = int(fitr);
  cout<<"Fit Status = "<<status<<endl;
  TH2D* hC2fit	= (TH2D*)hrho2->Clone("hC2");
  hC2fit	->SetTitle("Fit C_{2}(dy,d#varphi)/#Delta;dy;d#varphi");
  hC2fit	->Reset();
  TH2D* hC2fittri	= (TH2D*)hrho2->Clone("hC2");
  hC2fittri	->SetTitle("Fit C_{2}(dy,d#varphi);dy;d#varphi");
  hC2fittri	->Reset();
  TH2D* hC2fitrat	= (TH2D*)hrho2->Clone("hC2");
  hC2fitrat	->SetTitle("C_{2}(dy,d#varphi)/FIT;dy;d#varphi");
  hC2fitrat	->Reset();
  //
  for (int idy=1;idy<=axdy->GetNbins();idy++){
    for (int idphi=1;idphi<=axdphi->GetNbins();idphi++){
      double dy		= axdy	->GetBinCenter(idy);
      double dphi		= axdphi->GetBinCenter(idphi);
      double vC2		= hC2	->GetBinContent(idy,idphi);
      double vC2e		= hC2	->GetBinError  (idy,idphi);
      double valdytri	= (1.0 - fabs(dy)/(YU-YL));
      double vC2fit	= fC2	->Eval(dy,dphi);
      hC2fit		->SetBinContent(idy,idphi,vC2fit);
      hC2fittri	->SetBinContent(idy,idphi,vC2fit*valdytri);
      //hC2fitrat	->SetBinContent(idy,idphi,vC2fit/(vC2fit*valdytri));
      //if (fabs(dy)<0.1&&fabs(dphi)<30){
      //	cout<<dy<<" "<<dphi<<" "<<vC2<<" "<<vC2fit<<" "<<vC2/vC2fit<<endl;
      //}
    }
  }		  
  double ainte_C2tri		= 0;
  double aint_C2tri		= hC2tri	->IntegralAndError(1,0,1,0,ainte_C2tri   ,"");
  double ainte_C2fit		= 0;
  double aint_C2fit		= hC2fit	->IntegralAndError(1,0,1,0,ainte_C2fit   ,"");
  double ainte_C2fittri	= 0;
  double aint_C2fittri	= hC2fittri	->IntegralAndError(1,0,1,0,ainte_C2fittri,"");
  //
  //---- now for fun let's fit the R2 projections (page 1 lower right)
  TF1* fR2dy		= new TF1("fR2dy","[0]+[1]*exp(-0.5*(x/[2])**2)",DYL,DYU);
  fR2dy		->SetLineColor(2);
  fR2dy		->SetParameter(0,0.012);	// base
  fR2dy		->SetParameter(1,0.01);		// femto scale
  fR2dy		->SetParameter(2,0.2);		// femto sig dy
  hR2dy		->Fit("fR2dy","QMRS");
  TF1* fR2dphi	= new TF1("fR2dphi","[0]+[1]*cos(2.*x/57.2958)+[2]*exp(-0.5*(x/[3])**2)",DPHIL,DPHIU);
  fR2dphi	->SetLineColor(2);
  fR2dphi	->SetParameter(0,-0.01);	// base
  fR2dphi	->SetParameter(1,0.01);		// flow scale
  fR2dphi	->SetParameter(2,0.01);		// femto scale
  fR2dphi	->SetParameter(3,20.);		// femto sig dphi
  hR2dphi	->Fit("fR2dphi","QMRS");
  TF1* fR2dphi_a	= new TF1("fR2dphi_a","[0]*exp(-0.5*(x/[1])**2)",DPHIL,DPHIU);
  fR2dphi_a	->SetLineColor(16);
  fR2dphi_a	->SetParameter(0, fR2dphi->GetParameter(2) );	// femto scale
  fR2dphi_a	->SetParameter(1, fR2dphi->GetParameter(3) );	// femto sig dphi
  TF1* fR2dphi_b	= new TF1("fR2dphi_b","[0]+[1]*cos(2.*x/57.2958)",DPHIL,DPHIU);
  fR2dphi_b	->SetLineColor(16);
  fR2dphi_b	->SetParameter(0, fR2dphi->GetParameter(0) );	// base	
  fR2dphi_b	->SetParameter(1, fR2dphi->GetParameter(1) );	// flow scale

  //---- prepare to paint...
  //
  gROOT->SetStyle("Modern");
  gStyle->SetTitleAlign(13);
  gStyle->SetTitleX(0.10);
  gStyle->SetTitleY(0.995);
  gStyle->SetTitleFontSize(0.057);
  gStyle->SetOptStat(0);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadTopMargin(0.005);
  gStyle->SetPadBottomMargin(0.08);
  gStyle->SetPadLeftMargin(0.1);
  int ican=-1,iline=-1,ivline=-1; 
  int itxt=-1;
  TCanvas *ccan[1000];
  TLatex *txt[1000]; 
  for (int i=0;i<1000;i++){
    txt[i]	= new TLatex();
    txt[i]	->SetTextFont(42);
    txt[i]	->SetTextSize(0.035);
    txt[i]	->SetTextAlign(12);		// left-just, right-just is 32
    txt[i]	->SetNDC();
  }

  //---- 
  ++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,695);
  ccan[ican]->cd(); ccan[ican]->Divide(3,1,0.0001,0.0001);
  ccan[ican]->cd(1);
  XRangeFinder(hm0);
  gStyle->SetPalette(1);
  rp->SetH1DrawOpt("E1 P");
  rp->SetH2DrawOpt("L");
  rp->SetSeparationMargin(0.01);
  rp->Draw();
  rp->GetUpperPad()->SetLogy(1);
  rp->SetLeftMargin(0.1);
  rp->SetUpTopMargin(0.005);
  rp->SetRightMargin(0.005);
  //rp->SetLowBottomMargin(0.55);
  rp->GetUpperRefXaxis()->SetRangeUser(xlimits[0],xlimits[1]+20);
  rp->GetLowerRefXaxis()->SetRangeUser(xlimits[0],xlimits[1]+20);
  rp->GetLowerRefYaxis()->SetRangeUser(0.01,1.99);
  rp->GetLowYaxis()->SetNdivisions(6, 0, 0, kTRUE);
  //rp->SetConfidenceIntervalColors(4,7);
  TGraphErrors* grat	= (TGraphErrors*)rp->GetLowerRefGraph();
  grat->SetMarkerSize(0.8);
  grat->SetMarkerStyle(20);
  grat->SetMarkerColor(4);
  grat->SetLineColor(4);
  grat->SetLineWidth(2);
  rp->GetUpperPad()->cd();
  hfitnbd->Draw("same");
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.95,Form("NEVT(M)=%0.3f",(double)NEVT/1000000.));
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.91,Form("Sing K_{1} = %d",S_K1));
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.87,Form("PairA K_{1} = %d",Pa_K1));
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.83,Form("PairB K_{1} = %d",Pb_K1));
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.79,Form("nbd fit K_{1} = %f",fitK1));
  ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
  txt[itxt]->DrawLatex(0.97,0.74,Form("nbd fit K_{2} = %f",fitK2));
  for (int iordm1=0;iordm1<8;iordm1++){
    ++itxt; txt[itxt]->SetTextAlign(33); txt[itxt]->SetTextSize(0.04); 
    txt[itxt]->DrawLatex(0.97,0.67-0.04*iordm1,Form("%s = %.3f",cumnames[iordm1],K[iordm1]));
  }
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.40,Form("Poisson #chi^{2} = %0.3f",chi2_Poisson));
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.36,Form("Poisson NDF = %d",(int)ndof_Poisson));
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.32,Form("Poisson #chi^{2}/NDF = %0.3f",redchi2_Poisson));
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.28,Form("NBD #chi^{2} = %0.3f",chi2_NBD));
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.24,Form("NBD NDF = %d",(int)ndof_NBD));
  //++itxt; txt[itxt]->SetTextAlign(32); txt[itxt]->DrawLatex(0.65,0.20,Form("NBD #chi^{2}/NDF = %0.3f",redchi2_NBD));
  rp->GetLowerPad()->cd();
  hratnbd->Draw("HIST same");
  hratnbd->Draw("E same");
  //
  ccan[ican]->cd(2);
  TPad* thispad	= (TPad*)gPad;
  thispad->Divide(1,4,0.0001,0.0001);
  thispad->cd(1);
  hrho1->SetMinimum(0);
  hrho1->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.99,Form("INT=%.3f#pm%.3f",aint_rho1,ainte_rho1));
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.91,Form("f_{1}=%.3f",f1));
  thispad->cd(2);
  hrho2->SetMinimum(0);
  hrho2->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.99,Form("INT=%.3f#pm%.3f",aint_rho2,ainte_rho2));
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.91,Form("f_{2}=%.3f",f2));
  thispad->cd(3);
  hrho1rho1->SetMinimum(0);
  hrho1rho1->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.99,Form("INT=%.3f#pm%.3f",aint_rho11,ainte_rho11));
  ++itxt; txt[itxt]->SetTextSize(0.07); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.97,0.91,Form("f_{1}^{2}=%.3f",f1*f1));
  thispad->cd(4);
  hR2->Draw("lego2");
  ccan[ican]->cd(3);
  TPad* thispad3	= (TPad*)gPad;
  thispad3->Divide(1,4,0.0001,0.0001);
  thispad3->cd(1);
  hR2->Draw("lego2");
  thispad3->cd(2);
  hR2->Draw("colz");
  thispad3->cd(3);
  RangeFinderH(0,hR2dy);
  RangeCheck();
  hR2dy->SetMinimum(ylimits[0]);
  hR2dy->SetMaximum(ylimits[1]);
  hR2dy->Draw("HIST");
  fR2dy->Draw("same");
  hR2dy->Draw("HIST same");
  hR2dy->Draw("e same");
  thispad3->cd(4);
  RangeFinderH(0,hR2dphi);
  RangeCheck();
  if (ylimits[0]>-0.0025) ylimits[0] = -0.0025;
  hR2dphi->SetMinimum(ylimits[0]);
  hR2dphi->SetMaximum(ylimits[1]);
  hR2dphi->Draw("HIST");
  fR2dphi_a->Draw("same");
  fR2dphi_b->Draw("same");
  fR2dphi->Draw("same");
  hR2dphi->Draw("HIST same");
  hR2dphi->Draw("e same");
  //
  ccan[ican]->cd(); ccan[ican]->Update();
  ccan[ican]->Print(outfileO.Data());	

  //---- generator monitoring histograms...
  // 	++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,650);
  // 	ccan[ican]->cd(); ccan[ican]->Divide(3,2,0.0001,0.0001);
  // 		ccan[ican]->cd(1);
  // 			gPad->SetLogy(1);
  // 			XRangeFinder(hS);
  // 			hS->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hS->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hS->Draw("PE");
  // 		ccan[ican]->cd(2);
  // 			gPad->SetLogy(1);
  // 			XRangeFinder(hN);
  // 			hN->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hN->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hN->Draw("PE");
  // 		ccan[ican]->cd(3);
  // 			hv2dphi->SetMinimum(0);
  // 			hv2dphi->SetLineWidth(2);
  // 			hv2dphi->SetFillColor(18);
  // 			hv2dphi->Draw("hist");
  // 		ccan[ican]->cd(4);
  // 			gPad->SetLogy(1);
  // 			XRangeFinder(hPa);
  // 			hPa->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hPa->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hPa->Draw("PE");
  // 		ccan[ican]->cd(5);
  // 			gPad->SetLogy(1);
  // 			XRangeFinder(hPb);
  // 			hPb->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hPb->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  // 			hPb->Draw("PE");
  // 	ccan[ican]->cd(); ccan[ican]->Update();
  // 	ccan[ican]->Print(outfile.Data());	

  //---- correlator integration...
  ++ican; ccan[ican]	= new TCanvas(Form("ccan%d",ican),Form("ccan%d",ican),20*ican,30+20*ican,900,695);
  ccan[ican]->cd(); ccan[ican]->Divide(4,3,0.0001,0.0001);
  //---------- row 1
  ccan[ican]->cd(1);
  gPad->SetLogy(1);
  XRangeFinder(hS);
  hS->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hS->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hS->Draw("PE");
  ++itxt; txt[itxt]->SetTextSize(0.06); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.82,"Singles:");
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.72,Form("Gen K_{1} = %d",S_K1));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.66,Form("K_{1} = %.3f",hS->GetMean()));
  ccan[ican]->cd(2);
  gPad->SetLogy(1);
  XRangeFinder(hPa);
  hPa->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hPa->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hPa->Draw("PE");
  ++itxt; txt[itxt]->SetTextSize(0.06); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.82,"Pairs (CorrSource A):");
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.72,Form("Gen K_{1} = %d",Pa_K1));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.66,Form("K_{1} = %.3f",hPa->GetMean()));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.50,Form("#sigma(dy) = %.3f",Pa_sigdy));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.44,Form("#sigma(d#varphi) = %.3f",Pa_sigdphi));
  ccan[ican]->cd(3);
  gPad->SetLogy(1);
  XRangeFinder(hPb);
  hPb->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hPb->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hPb->Draw("PE");
  ++itxt; txt[itxt]->SetTextSize(0.06); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.82,"Pairs (CorrSource B):");
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.72,Form("Gen K_{1} = %d",Pb_K1));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.66,Form("K_{1} = %.3f",hPb->GetMean()));
  ccan[ican]->cd(4);
  gPad->SetLogy(1);
  XRangeFinder(hN);
  hN->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hN->GetXaxis()->SetRangeUser(xlimits[0],xlimits[1]);
  hN->Draw("PE");
  ++itxt; txt[itxt]->SetTextSize(0.06); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.82,"Total:");
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.72,Form("Gen K_{1} = %d",N_target));
  ++itxt; txt[itxt]->SetTextSize(0.05); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.82,0.66,Form("K_{1} = %.3f",hN->GetMean()));
  //---------- row 2
  ccan[ican]->cd(5);
  hC2->Draw("colz");
  ccan[ican]->cd(6);
  hC2->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.99,Form("2(Pa_K1+Pb_K1) = %d",2*(Pa_K1+Pb_K1)));
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.95,Form("INT C_{2} = %.3f #pm %.4f",aint_C2,ainte_C2));
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.91,Form("K_{2}-K_{1} = %.3f",K[1]-K[0]));
  ccan[ican]->cd(7);
  hC2tri->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.99,Form("INT = %.3f #pm %.4f",aint_C2tri,ainte_C2tri));
  ccan[ican]->cd(8);
  hC2tri->Draw("colz");
  //---------- row 3
  ccan[ican]->cd(9);
  //hC2fitrat->SetMinimum(0.5);
  //hC2fitrat->SetMaximum(1.5);
  //hC2fitrat->Draw("lego2");
  hC2fittri->Draw("colz");
  ccan[ican]->cd(10);
  hC2fittri->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.99,Form("INT = %.3f #pm %.4f",aint_C2fittri,ainte_C2fittri));
  ccan[ican]->cd(11);
  hC2fit->Draw("lego2");
  ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
  txt[itxt]->DrawLatex(0.98,0.99,Form("INT = %.3f #pm %.4f",aint_C2fit,ainte_C2fit));
  for (int i=0;i<5;i++){
    ++itxt; txt[itxt]->SetTextSize(0.04); txt[itxt]->SetTextAlign(33);  
    txt[itxt]->DrawLatex(0.98,0.95-i*0.04,Form("%s = %.3f #pm %.4f",fitC2_parname[i],fitC2_par[i],fitC2_pare[i]));
  }
  ccan[ican]->cd(12);
  hC2fit->Draw("colz");
  ccan[ican]->cd(); ccan[ican]->Update();
  ccan[ican]->Print(outfile.Data());	
  //---- build and scp the pdf file...
	
  if (ican>=0){
    ccan[ican]->Print(outfileC.Data());
    TString sexec;
    sexec	 = TString(Form("pstopdf %s %s",outfile.Data(),outfileP.Data()));
    cout<<sexec.Data()<<endl;
    gSystem->Exec(sexec.Data());
    sexec	 = TString(Form("/bin/rm %s",outfile.Data()));
    cout<<sexec.Data()<<endl;
    gSystem->Exec(sexec.Data());
    //sexec	 = TString(Form("/usr/bin/scp %s wjllope@rhic22.physics.wayne.edu:/Library/WebServer/WebPages/files/",outfileP.Data()));
    //cout<<sexec.Data()<<endl;
    //gSystem->Exec(sexec.Data());
  }

}

void XRangeFinder(TH1* h){
  int nbx		= h->GetNbinsX();
  xlimits[0]	= 9999999;
  xlimits[1]	=       0;
  for (int ibx=1;ibx<=nbx;ibx++){
    int 	thisent	= h->GetBinContent(ibx);
    float	thisx	= h->GetXaxis()->GetBinCenter(ibx);
    if (xlimits[0]>9999990 && thisent>0){ xlimits[0]=thisx; }
    if (                      thisent>0){ xlimits[1]=thisx; }
  }
  XRangeCheck();
}
// void XRangeFinder2D(TH2* h){
// 	int nbx		= h->GetNbinsX();
// 	xlimits[0]	= 9999999;
// 	xlimits[1]	=       0;
// 	for (int ibx=1;ibx<=nbx;ibx++){
// 		TH1D*	h1		= (TH1D*)h->ProjectionY("h1",ibx,ibx);
// 		int 	thisent	= h1->GetEntries();
// 		float	thisx	= h->GetXaxis()->GetBinCenter(ibx);
// 		if (xlimits[0]>9999990 && thisent>0){ xlimits[0]=thisx; }
// 		if (                      thisent>0){ xlimits[1]=thisx; }
// 		delete h1;
// 	}
// 	//
// 	XRangeCheck();
// 	//
// }
void XRangeCheck(){
  if (xlimits[0]>0){ 			// lower limit
    xlimits[0]	*= 0.5;
  } else if (xlimits[0]==0){
    xlimits[0]	= -1;
  } else if (xlimits[0]<0){
    xlimits[0]	*= 1.1;
  }
  if (xlimits[1]>0){ 			// upper limit
    xlimits[1]	*= 1.15;
  } else if (xlimits[1]==0){
    xlimits[1]	=    1;
  } else if (xlimits[1]<0){
    xlimits[1]	*= 0.9;
  }
}
// void RangeFinder2D(TH2* h){
// 	int nbx	= h->GetNbinsX();
// 	int nby	= h->GetNbinsY();
// 	ylimits[0]	= 999999;
// 	ylimits[1]	=      0;
// 	for (int ibx=1;ibx<=nbx;ibx++){
// 		TH1D *h1	= (TH1D*)h->ProjectionY("h1",ibx,ibx);
// 		if (h1->GetEntries()>0){
// 			for (int iby=1;iby<=nby;iby++){
// 				int nent	= h1->GetBinContent(iby);
// 				if (nent==0) continue;
// 				float ybin	= h1->GetXaxis()->GetBinCenter(iby);
// 				if (ybin<ylimits[0]){ ylimits[0] = ybin; }
// 				if (ybin>ylimits[1]){ ylimits[1] = ybin; }
// 			}
// 		}
// 		delete h1;
// 	}
// 	//
// 	RangeCheck();
// 	//
// }
// void RangeCheck(){
// 	if (ylimits[0]>0){ 
// 		ylimits[0]	*= 0.9;
// 	} else if (ylimits[0]==0){
// 		ylimits[0]	=    -1;
// 	} else if (ylimits[0]<0){
// 		ylimits[0]	*= 1.1;
// 	}
// 	if (ylimits[1]>0){ 
// 		ylimits[1]	*= 1.1;
// 	} else if (ylimits[1]==0){
// 		ylimits[1]	=     1;
// 	} else if (ylimits[1]<0){
// 		ylimits[1]	*= 0.9;
// 	}
// }
//----------------------------------------------------------------------------
void RangeFinderH(int isteer, TH1D* h){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  if (!h) return;
  if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  if (debug) cout<<h->GetName()<<" "<<h->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= h->GetNbinsX();
  for (int ibx=1;ibx<=nbx;ibx++){
    double val	= h->GetBinContent(ibx);
    double vale	= h->GetBinError(ibx);
    //if (val==0) continue;
    if (val-vale<ylimits[0]){ ylimits[0] = val-vale; }	// update lower limit...
    if (val+vale>ylimits[1]){ ylimits[1] = val+vale; }	// update upper limit...
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeFinderG(int isteer, TGraph* g){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  //if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  //if (debug) cout<<g->GetName()<<" "<<g->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= g->GetN();
  for (int ibx=0;ibx<nbx;ibx++){
    double x,val;
    g->GetPoint(ibx,x,val);
    //if (val==0) continue;
    if (val<ylimits[0]){ ylimits[0] = val; }
    if (val>ylimits[1]){ ylimits[1] = val; }
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeFinderGE(int isteer, TGraphErrors* g){
  bool debug = false;
  if (isteer>=100){
    debug	 = true;
    isteer	-= 100;
    cout<<"debug on... "<<isteer<<endl;
  }
  if (isteer==0){
    ylimits[0]	=  999999;
    ylimits[1]	= -999999;
  }
  //if (debug) cout<<"start... "<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  //if (debug) cout<<g->GetName()<<" "<<g->GetEntries()<<endl;
  //if (h->GetEntries()<1){
  //	ylimits[0]	=  0;
  //	ylimits[1]	=  1;
  //	return;
  //}
  int nbx	= g->GetN();
  for (int ibx=0;ibx<nbx;ibx++){
    double x,val,vale;
    g->GetPoint(ibx,x,val);
    vale	= g->GetErrorY(ibx);
    //if (val==0) continue;
    if (val-vale<ylimits[0]){ ylimits[0] = val-vale; }
    if (val+vale>ylimits[1]){ ylimits[1] = val+vale; }
    if (debug) cout<<"check... "<<ibx<<" val="<<val<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
  }
  if (debug) cout<<"FINAL... "<<isteer<<" min="<<ylimits[0]<<" max="<<ylimits[1]<<endl;
}
//----------------------------------------------------------------------------
void RangeCheck(){
  if (ylimits[0]>0){ 
    ylimits[0]	*= 0.9;
  } else if (ylimits[0]==0){
    //		ylimits[0]	= -1;
  } else if (ylimits[0]<0){
    ylimits[0]	*= 1.1;
  }
  if (ylimits[1]>0){ 
    ylimits[1]	*= 1.1;
  } else if (ylimits[1]==0){
    //		ylimits[1]	=    1;
  } else if (ylimits[1]<0){
    ylimits[1]	*= 0.9;
  }
}



//-------------------------------------------------------------------------------------
//
//	Service routine to project any 2D R2 onto a 1D R2 projection
//	Called from Loop right after SetFLip (crossing correction)...
//	This routine assumes that the 1D is properly binned for the assumed projection! 
//
//	examples:
//			GenericProject(hR2dydphi,"Avg","Y",h1)		--> h1 is R2(dphi) from (dy,dphi), crossing correctable
//			GenericProject(hR2dydphi,"Avg","X",h1)		--> h1 is R2(dy) from (dy,dphi), crossing correctable
//			GenericProject(hR2y1y2,"Avg","D",h1)		--> h1 is R2(dy) from (y1,y2), NOT crossing correctable
//			GenericProject(hR2y1y2,"Avg","A",h1)		--> h1 is R2(<y>) from (y1,y2), NOT crossing correctable
//
//	now includes protection on h1 (should be empty at call!) !!
//
void GenericProject(TH2D* h2,TString steerav,TString steerax,TH1D* h1){
  //
  if (!h1){ cout<<"reader::GenericProject -- no h1 pointer, cannot project!"<<endl; exit(0); }
  if (!h2){ cout<<"reader::GenericProject -- no h2 pointer, cannot project!"<<endl; exit(0); }
  bool doadd=false, doavg=false;
  doadd	= steerav.Contains("Add",TString::kIgnoreCase);	
  doavg	= steerav.Contains("Avg",TString::kIgnoreCase);	
  if (!doadd && !doavg){ cout<<"reader::GenericProject steering issue."<<endl; exit(0); }
  if ( doadd &&  doavg){ cout<<"reader::GenericProject steering issue."<<endl; exit(0); }
  int kProjAxis	= -1;
  if (steerax.Contains("X",TString::kIgnoreCase)) kProjAxis = 0;	// project onto x axis
  if (steerax.Contains("Y",TString::kIgnoreCase)) kProjAxis = 1;	// project onto y axis
  if (steerax.Contains("A",TString::kIgnoreCase)) kProjAxis = 2;	// project onto axis-avg axis
  if (steerax.Contains("D",TString::kIgnoreCase)) kProjAxis = 3;	// project onto axis-diff axis
  int kR2			=  0;
  TString histtit	= (TString)h2->GetTitle();
  if (histtit.Contains("y_{1},y_{2}",TString::kIgnoreCase)) kR2 = 0;	// h2 is (y1,y2);
  if (histtit.Contains("dy,d#phi"   ,TString::kIgnoreCase)) kR2 = 1;	// h2 is (dy,dphi);
  if (histtit.Contains("dy,dq"      ,TString::kIgnoreCase)) kR2 = 2;	// h2 is (dy,dq);
  //kR2	= 1;	// hardwired in this B2 macro
  //
  h1			->Reset();
  h1->SetMaximum(1.0);
  h1->SetMinimum(-1.0);
  TH1D *h1e	= (TH1D*)h1->Clone("h1e");
  TH1D *h1N	= (TH1D*)h1->Clone("h1N");
  //
  if (kProjAxis==0){		// project onto X-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
      for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
	double x,y,a,d,val,vale;
	x		= h2->GetXaxis()->GetBinCenter(ibx);
	val		= h2->GetBinContent(ibx,iby);
	vale	= h2->GetBinError(ibx,iby);
	h1		->Fill(x, val);
	h1e		->Fill(x, vale*vale);
	h1N		->Fill(x, 1.0);
      }	// iby
    }	// ibx
    //
  } else if (kProjAxis==1){		// project onto Y-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
      for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
	double x,y,a,d,val,vale;
	y		= h2->GetYaxis()->GetBinCenter(iby);
	val		= h2->GetBinContent(ibx,iby);
	vale	= h2->GetBinError(ibx,iby);
	h1		->Fill(y, val);
	h1e		->Fill(y, vale*vale);
	h1N		->Fill(y, 1.0);
      }	// iby
    }	// ibx
  } else if (kProjAxis==2){		// project onto avg-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
      for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
	double x,y,a,d,val,vale;
	a		= (h2->GetXaxis()->GetBinCenter(ibx) + h2->GetYaxis()->GetBinCenter(iby))/2.;
	val		= h2->GetBinContent(ibx,iby);
	vale	= h2->GetBinError(ibx,iby);
	h1		->Fill(a, val);
	h1e		->Fill(a, vale*vale);
	h1N		->Fill(a, 1.0);
      }	// iby
    }	// ibx
  } else if (kProjAxis==3){		// project onto diff-axis of h2
    for (int ibx=1;ibx<=h2->GetXaxis()->GetNbins();ibx++){
      for (int iby=1;iby<=h2->GetYaxis()->GetNbins();iby++){
	double x,y,a,d,val,vale;
	d		= h2->GetXaxis()->GetBinCenter(ibx) - h2->GetYaxis()->GetBinCenter(iby);
	val		= h2->GetBinContent(ibx,iby);
	vale	= h2->GetBinError(ibx,iby);
	h1		->Fill(d, val);
	h1e		->Fill(d, vale*vale);
	h1N		->Fill(d, 1.0);
      }	// iby
    }	// ibx
  }
  //
  //---- now finally fill the projection plot h1
  for (int ibi=1;ibi<=h1->GetNbinsX();ibi++){
    double v,e,n;
    v	=       h1	->GetBinContent(ibi);
    e	= sqrt( h1e	->GetBinContent(ibi));
    n	=       h1N	->GetBinContent(ibi);
    if (n>0){
      if (doavg){
	h1->SetBinContent(ibi,v/n);					
	h1->SetBinError(ibi,e/n);
      } else if (doadd){
	h1->SetBinContent(ibi,v);					
	h1->SetBinError(ibi,e);
      }
    }
  }
  //
  delete h1e;
  delete h1N;
  //
}

double binomial_pdf(unsigned int k, double p, unsigned int n) {
  if (k > n) {
    return 0.0;
  } else {
    double coeff = ROOT::Math::lgamma(n+1) - ROOT::Math::lgamma(k+1) - ROOT::Math::lgamma(n-k+1);
    return std::exp(coeff + k * std::log(p) + (n - k) * ROOT::Math::log1p(-p));
  }
}
double negative_binomial_pdf(unsigned int k, double p, double n) {
  // impelment in term of gamma function
  if (n < 0)  return 0.0;
  if (p < 0 || p > 1.0) return 0.0;
  double coeff = ROOT::Math::lgamma(k+n) - ROOT::Math::lgamma(k+1.0) - ROOT::Math::lgamma(n);
  return std::exp(coeff + n * std::log(p) + double(k) * ROOT::Math::log1p(-p));
}
double poisson_pdf(unsigned int n, double mu) {
  if (n >  0){ 
    return std::exp (n*std::log(mu) - ROOT::Math::lgamma(n+1) - mu);
  } else  {
    //  when  n = 0 and mu = 0,  1 is returned
    if (mu >= 0) return std::exp(-mu);
    // return a nan for mu < 0 since it does not make sense
    return std::log(mu);
  } 
}
Double_t fBD(Double_t *x, Double_t *par) {
  int k	= x[0];
  double vpdf	= binomial_pdf(k,par[1],par[2]);
  return par[0]*vpdf;
}
Double_t fPoi(Double_t *x, Double_t *par) {
  int k	= x[0];
  double vpdf = poisson_pdf(k,par[1]);
  return par[0]*vpdf;
}
Double_t fNBD(Double_t *x, Double_t *par) {
  int k	= x[0];
  double vpdf	= negative_binomial_pdf(k,par[1],par[2]);
  return par[0]*vpdf;
}
Double_t fC2fit(Double_t *x, Double_t *par) {
  //
  double dy	= x[0];
  //double dytri= 1.0 - fabs(dy)/(YU-YL);
  double dphi	= x[1];
  double base	= par[0];	//0 global constant offset
  double a_h 	= par[1];	//A scale
  double a_mx	= 0.0;		//A posn of src peak in dy
  double a_sx	= par[2];	//A sigma-x (dy)
  double a_my	= 0.0;		//A posn of src peak in dphi
  double a_sy	= par[3];	//A sigma-y (dphi)
  double b_h	= par[4];	//B scale
  //
  double corA	= a_h*exp(-0.5*(((dy-a_mx)/a_sx)*((dy-a_mx)/a_sx)+((dphi-a_my)/a_sy)*((dphi-a_my)/a_sy)));
  double corB	= b_h*cos(2.*dphi/57.2958);
  //double val	= dytri*(base + corA + corB);
  double val	= base + corA + corB;
  //
  return val;
  //
}

