#ifndef __MultCum_h__
#define __MultCum_h__

#include <iostream>
#include "TROOT.h"
#include "TH1D.h"
#include "TMath.h"
using namespace std;

class MultCum{
	//
	public:
		//
		MultCum(TH1D* h=0);		// ctor
		virtual ~MultCum();		// dtor
		//
		double GetCumulant(int iord=0);
		//
	//
	private:
		//
		void setCumulants();
		//
		int mNCumulants;
		TH1D *mH;
		double mCumulants[20];
		//
	//
	//ClassDef(MultCum, 0)
};

inline double MultCum::GetCumulant(int iord){ return mCumulants[iord-1]; }

#endif
