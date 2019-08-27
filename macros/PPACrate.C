#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include "TChain.h"
#include "TCutG.h"
#include <TFile.h>
#include "TCanvas.h"
#include "TH2D.h"
#include "TH2.h"
#include "cmath"
#include "Riostream.h"

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamRaw.h"
#include "BeamBeam.h"

using namespace std;

const char *cutfname = NULL;
bool filecalled = false;



////This macro counts the hits in the TDC spectra for the F72 double PPACs. This can be used to infer beam rate
void PPACrate(){
  //style options
  gROOT->ForceStyle();
  //gStyle->SetOptStat(0);
  gStyle->SetOptDate(0);


  double TDCmax=164000.;
    
  //first and last run
  int first=2174;
  int last=3184;


  int i = last-first+1;
  int nafter[i];
  int npeak[i];
  int nbefore[i];
  double fPeak[i];
  double fSigma[i];
  
  int ii=0;
  TH1D *histSpectra=new TH1D("histSpectra","",1000,0,TDCmax);
  int numEvents=0;
    
    for(int runNo=first; runNo<=last; runNo++){    
      nafter[ii]=0;
      npeak[ii]=0;
      nbefore[ii]=0;
      fPeak[ii]=0;
      fSigma[ii]=0;

      ofstream myfile;
      if(runNo==first){
	myfile.open (Form("output/ppacTDCcount.%i.%i.csv",first,last));
      }
      if(runNo>first){
	myfile.open (Form("output/ppacTDCcount.%i.%i.csv",first,last),std::ios_base::app);
      }
      
      ifstream checkfile(Form("ridf/SMDAQ%d.ridf", runNo));
      if(!checkfile){
	ii++;
	myfile << Form("%i,%i,%i,%i,%i,%f,%f",runNo,numEvents,nbefore[ii],npeak[ii],nafter[ii],fPeak[ii],fSigma[ii]) << endl;
	continue;
      }

      auto eventStore = new TArtEventStore();
      TString ridfFile = Form("ridf/SMDAQ%d.ridf", runNo);
      eventStore -> Open(ridfFile.Data());
      
      auto storeMan = TArtStoreManager::Instance();
      TArtRawEventObject *rawevent = eventStore->GetRawEventObject();
      numEvents=0;
      while (eventStore -> GetNextEvent() && numEvents<100000) {
	numEvents++;
	if(true){
	  for(int i=0;i<rawevent->GetNumSeg();i++){
	    TArtRawSegmentObject *seg = rawevent->GetSegment(i);
	    int detector = seg->GetDetector();
	    if(detector == 1) {//PPACT is 1 
	      for(int j=0;j<seg->GetNumData();j++){
		TArtRawDataObject *d = seg->GetData(j);
		int ch = d->GetCh();
		int val = d->GetVal();
		if(val>164000) continue;
		int ind=(ch-16)/4;
		
		if(ch>=56 && ch<64){
		  histSpectra->Fill(val);
		}
	      }//for data segment
	    } // end of PPAC
	  }	  
	}//optional condition
      }//while event loop
    ii++;
    double peak;
    peak=(double)(histSpectra->GetMaximumBin())*TDCmax/1000.;
    cout << peak << endl;
    TF1 *fitPeak=new TF1("fitPeak","gaus",peak-2.,peak+2.);
    fitPeak->SetParameters(histSpectra->GetMaximum(),peak,1);
    histSpectra->Fit(fitPeak,"Q","",peak-500.,peak+500.);
    peak=fitPeak->GetParameter(1);
    double sigma=fitPeak->GetParameter(2);
    //cout << peak << ',' << sigma << endl;
    //cout << histSpectra->Integral((peak-5.*sigma)*1000./TDCmax,(peak+5.*sigma)*1000./TDCmax) << endl;
    //cout << histSpectra->Integral((peak+15.*sigma)*1000./TDCmax,(peak+15.*sigma+40000)*1000./TDCmax) << endl;
    //cout << histSpectra->Integral(0,(peak-15.*sigma)*1000./TDCmax) << endl;
    nbefore[ii]=(int)histSpectra->Integral(0,(peak-15.*sigma)*1000./TDCmax);
    npeak[ii]=(int)histSpectra->Integral((peak-5.*sigma)*1000./TDCmax,(peak+5.*sigma)*1000./TDCmax);
    nafter[ii]=(int)histSpectra->Integral((peak+15.*sigma)*1000./TDCmax,(peak+15.*sigma+10000)*1000./TDCmax);

    fPeak[ii]=peak;
    fSigma[ii]=sigma;



    histSpectra->Reset();
    cout << Form("%i,%i,%i,%i,%f,%f",runNo,nbefore[ii],npeak[ii],nafter[ii],fPeak[ii],fSigma[ii]) << endl;
    
    myfile << Form("%i,%i,%i,%i,%i,%f,%f",runNo,numEvents,nbefore[ii],npeak[ii],nafter[ii],fPeak[ii],fSigma[ii]) << endl;



    }//for run

    
}
