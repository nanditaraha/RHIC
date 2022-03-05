#include "MultCum.h"

//------------------------------------------------------------
MultCum::MultCum(TH1D* h){

	mNCumulants = 20;
	for (int i=0;i<mNCumulants;i++){
		mCumulants[i]	= 0.;
	}

	if (h==0){ cout<<"no histogram"<<endl; }

	mH	= h;
	setCumulants();
}

//------------------------------------------------------------
MultCum::~MultCum(){
}

//------------------------------------------------------------
//double MultCum::GetCumulant(int iord){
//	return mCumulants[iord];
//}

//------------------------------------------------------------
//------------------------------------------------------------
void MultCum::setCumulants(){
	//
	if (!mH){ cout<<"MultCum::setCumulants hist pointer problem... "<<endl; exit(0); }
	if (mH->Integral()==0){ return; }
	TH1D *mHsc	= (TH1D*)mH->Clone("mHsc");
		mHsc	->Scale(1./mH->Integral());
	//
	double adev1	= 0.0;
	double adev2	= 0.0;
	double adev3	= 0.0;
	double adev4	= 0.0;
	double adev5	= 0.0;
	double adev6	= 0.0;
	double den		= 0.0;
	double mean		= mH->GetMean();
	double f1		= 0;
	double f2		= 0;
	double f3		= 0;
	//
	//---- get <(dN)^n>....
	int nbx			 = mH->GetNbinsX();
	for (int ibx=1;ibx<=nbx;ibx++){
		double x 	 = mH->GetBinCenter(ibx);
		double w 	 = mH->GetBinContent(ibx);
		adev1		+= w*x;
		adev2		+= w*(x-mean)*(x-mean);
		adev3		+= w*(x-mean)*(x-mean)*(x-mean);
		adev4		+= w*(x-mean)*(x-mean)*(x-mean)*(x-mean);
		adev5		+= w*(x-mean)*(x-mean)*(x-mean)*(x-mean)*(x-mean);
		//adev6		+= w*(x-mean)*(x-mean)*(x-mean)*(x-mean)*(x-mean)*(x-mean);
		den			+= w;
		//
		double pri	 = mHsc->GetBinContent(ibx);
		if (x>-0.5) f1	+= pri*x;					// <N>
		if (x> 0.5) f2	+= pri*x*(x-1.);			// <N(N-1)>
		if (x> 1.5) f3	+= pri*x*(x-1.)*(x-2.);		// <N(N-1)(N-2)>
	}
	if (den==0.0){ cout<<"MultCum::setCumulants denominator is zero"<<endl; }
	adev1	/= den;
	adev2	/= den;
	adev3	/= den;
	adev4	/= den;
	adev5	/= den;
	//adev6	/= den;
	//
	mCumulants[0]		= adev1;					// C1
	mCumulants[1]		= adev2;					// C2
	mCumulants[2]		= adev3;					// C3
	mCumulants[3]		= adev4 - 3.*adev2*adev2;	// C4
	//
	if (adev1!=0.0){
		mCumulants[5]	= adev3/adev1;				// C3/C1
	} else {
		mCumulants[5]	= 0.0;
	}
	if (adev2!=0.0){
		mCumulants[4]	= adev2/adev1;				// C2/C1
		mCumulants[6]	= adev3/adev2;				// C3/C2
		mCumulants[7]	= mCumulants[3]/adev2;		// C4/C2
	} else {
		mCumulants[4]	= 0.0;
		mCumulants[6]	= 0.0;
		mCumulants[7]	= 0.0;
	}
	if (adev3!=0.0){
		mCumulants[8]	= mCumulants[3]/adev3;		// C4/C3
	} else {
		mCumulants[8]	= 0.0;
	}
	mCumulants[9]	= f1;
	mCumulants[10]	= f2;
	mCumulants[11]	= f3;
	//
	//for (int i=0;i<mNCumulants;i++){
	//	if (TMath::IsNaN(mCumulants[i])){
	//		cout<<"NaN, cumulant "<<i<<"  "<<mH->GetName()<<endl;
	//	}
	//}
	
}
