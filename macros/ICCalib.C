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
#include "TF1.h"
#include "TSpectrum.h"


// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "BeamRaw.h"
#include "BeamBeam.h"
#include "CalibSpectrum.h"

using namespace std;

/*
// 132Sn
double AoQmin = 2.6;
double AoQmax = 2.72;
*/
// 108Sn
double AoQmin = 2;
double AoQmax = 2.5;



//124Sn
//double AoQmin = 2.33;
//double AoQmax = 2.65;

double Zetmin = 46;
double Zetmax = 54;

double Betamin = 0.6;
double Betamax = 0.7;

double ICMeVmin = 2.;
double ICMeVmax = 3.;


void ICCalib(){
    gROOT->ForceStyle();
    
    //_________________________________________________________________________________________________________
    // Reading the root files, and plotting:
    //            PID (Z vs A/Q)                    - to look at the PID before calibration
    //            ADC[i] vs ADC[0] for IC           - Check if pedestal is 0, then no need to change ch2mev_1
    //            Z vs Beta dependence              - For ch2mev_0, if Beta dependence seen, need to correct it
    //BeamBeam *beam = new BeamBeam();
    //BeamRaw *raw = new BeamRaw();

    
    //ROOT Files
    //for(int ii=2840; ii<=3039; ii++){  // for 132Sn
    //for(int ii=2300; ii<=2400; ii++){  // for 124Sn

      BeamBeam *beam = new BeamBeam();
      BeamRaw *raw = new BeamRaw();
      
      int ii = 2257; //132Sn test run
          beam->fChainBeam->AddFile(Form("data/run%i.ridf.root",ii),0,"beam");
        raw->fChain->AddFile(Form("data/run%i.ridf.root",ii),0,"raw");
	//}
    beam->Init();
    raw->Init();

    Long64_t nentries = raw->fChain->GetEntriesFast();
    cout << " Entries = " << nentries << endl;
    
    //Canvases
    TCanvas *cvsADC = new TCanvas("cvsADC", "", 0, 0, 1400, 700);
    TCanvas *cvsPID = new TCanvas("cvsPID", "", 0, 0, 1000, 600);
    TCanvas *cvsZet = new TCanvas("cvsZet", "", 0, 0, 800, 600);

    //Histograms
    auto *histIC_ADC1vsADC0 = new TH2D("histIC_ADC1vsADC0", "histIC_ADC1vsADC0", 2000, 0, 8000, 2000, 0, 8000);
    auto *histIC_ADC2vsADC0 = new TH2D("histIC_ADC2vsADC0", "histIC_ADC2vsADC0", 2000, 0, 8000, 2000, 0, 8000);
    auto *histIC_ADC3vsADC0 = new TH2D("histIC_ADC3vsADC0", "histIC_ADC3vsADC0", 2000, 0, 8000, 2000, 0, 8000);
    auto *histIC_ADC4vsADC0 = new TH2D("histIC_ADC4vsADC0", "histIC_ADC4vsADC0", 2000, 0, 8000, 2000, 0, 8000);
    auto *histIC_ADC5vsADC0 = new TH2D("histIC_ADC5vsADC0", "histIC_ADC5vsADC0", 2000, 0, 8000, 2000, 0, 8000);
    
    auto *histBeamPID = new TH2D("histBeamPID", "histBeamPID", 1500, AoQmin, AoQmax, 1500, Zetmin, Zetmax);

    auto *histZetBeta = new TH2D("histZetBeta", "histZetBeta", 500, Betamin, Betamax, 500, Zetmin, Zetmax);

    auto *histZet = new TH1D("histZet", "histZet", 500, 45, 55);
    auto *histICMeV = new TH1D("histICMeV", "histICMeV", 200, ICMeVmin, ICMeVmax);

    auto *histIC_ADC_sq = new TH1D("histIC_ADC_sq", "histIC_ADC_sq", 2000,1000,8000);
    
    
    //////////////////////////// Entries loop ////////////////////////////
    for(int ientry=0; ientry < nentries; ientry++){
      //if(ientry%1000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
        beam->fChainBeam->GetEvent(ientry);
        raw->fChain->GetEvent(ientry);

        double AoQ = beam->aoq;
        double Zet = beam->z;
        double beta = beam->beta;
        double ICMeVSqSum = raw->BigRIPSIC_fCalMeVSqSum[2];
        
        double ionpair = raw->BigRIPSIC_ionpair[2];
        double de_v = TMath::Log(ionpair*beta*beta) - TMath::Log((1-beta*beta)) - beta*beta;
        
        double ICADC[6];
        for(int ii=0; ii<6; ii++){
            ICADC[ii] = raw->BigRIPSIC_fADC[2][ii];
        }
       
        histIC_ADC1vsADC0->Fill( ICADC[0],ICADC[1] );
        histIC_ADC2vsADC0->Fill( ICADC[0],ICADC[2] );
        histIC_ADC3vsADC0->Fill( ICADC[0],ICADC[3] );
        histIC_ADC4vsADC0->Fill( ICADC[0],ICADC[4] );
        histIC_ADC5vsADC0->Fill( ICADC[0],ICADC[5] );
        
        histBeamPID->Fill(AoQ,Zet);
        
        histZetBeta->Fill(beta,Zet);
        
        histZet->Fill(Zet);
        histICMeV->Fill( TMath::Sqrt(ICMeVSqSum/de_v)*beta );
    }

    
    //Plotting the histograms
    cvsADC->cd();
    cvsADC->Divide(3,2);
    cvsADC->cd(1);
    histIC_ADC1vsADC0->Draw("colz");
    histIC_ADC1vsADC0->GetXaxis()->SetTitle("F7IC ADC[0]");
    histIC_ADC1vsADC0->GetXaxis()->SetTitleSize(0.05);
    histIC_ADC1vsADC0->GetXaxis()->SetTitleOffset(0.75);
    histIC_ADC1vsADC0->GetXaxis()->SetLabelSize(0.04);
    histIC_ADC1vsADC0->GetXaxis()->SetTitleFont(22);
    histIC_ADC1vsADC0->GetXaxis()->SetLabelFont(132);
    histIC_ADC1vsADC0->GetXaxis()->CenterTitle(true);
    histIC_ADC1vsADC0->GetYaxis()->SetTitle("F7IC ADC[1]");
    histIC_ADC1vsADC0->GetYaxis()->SetTitleSize(0.05);
    histIC_ADC1vsADC0->GetYaxis()->SetTitleOffset(0.75);
    histIC_ADC1vsADC0->GetYaxis()->SetLabelSize(0.04);
    histIC_ADC1vsADC0->GetYaxis()->SetTitleFont(22);
    histIC_ADC1vsADC0->GetYaxis()->SetLabelFont(132);
    histIC_ADC1vsADC0->GetYaxis()->CenterTitle(true);
    cvsADC->cd(2);
    histIC_ADC2vsADC0->Draw("colz");
    histIC_ADC2vsADC0->GetXaxis()->SetTitle("F7IC ADC[0]");
    histIC_ADC2vsADC0->GetXaxis()->SetTitleSize(0.05);
    histIC_ADC2vsADC0->GetXaxis()->SetTitleOffset(0.75);
    histIC_ADC2vsADC0->GetXaxis()->SetLabelSize(0.04);
    histIC_ADC2vsADC0->GetXaxis()->SetTitleFont(22);
    histIC_ADC2vsADC0->GetXaxis()->SetLabelFont(132);
    histIC_ADC2vsADC0->GetXaxis()->CenterTitle(true);
    histIC_ADC2vsADC0->GetYaxis()->SetTitle("F7IC ADC[2]");
    histIC_ADC2vsADC0->GetYaxis()->SetTitleSize(0.05);
    histIC_ADC2vsADC0->GetYaxis()->SetTitleOffset(0.75);
    histIC_ADC2vsADC0->GetYaxis()->SetLabelSize(0.04);
    histIC_ADC2vsADC0->GetYaxis()->SetTitleFont(22);
    histIC_ADC2vsADC0->GetYaxis()->SetLabelFont(132);
    histIC_ADC2vsADC0->GetYaxis()->CenterTitle(true);
    cvsADC->cd(3);
    histIC_ADC3vsADC0->Draw("colz");
    histIC_ADC3vsADC0->GetXaxis()->SetTitle("F7IC ADC[0]");
    histIC_ADC3vsADC0->GetXaxis()->SetTitleSize(0.05);
    histIC_ADC3vsADC0->GetXaxis()->SetTitleOffset(0.75);
    histIC_ADC3vsADC0->GetXaxis()->SetLabelSize(0.04);
    histIC_ADC3vsADC0->GetXaxis()->SetTitleFont(22);
    histIC_ADC3vsADC0->GetXaxis()->SetLabelFont(132);
    histIC_ADC3vsADC0->GetXaxis()->CenterTitle(true);
    histIC_ADC3vsADC0->GetYaxis()->SetTitle("F7IC ADC[3]");
    histIC_ADC3vsADC0->GetYaxis()->SetTitleSize(0.05);
    histIC_ADC3vsADC0->GetYaxis()->SetTitleOffset(0.75);
    histIC_ADC3vsADC0->GetYaxis()->SetLabelSize(0.04);
    histIC_ADC3vsADC0->GetYaxis()->SetTitleFont(22);
    histIC_ADC3vsADC0->GetYaxis()->SetLabelFont(132);
    histIC_ADC3vsADC0->GetYaxis()->CenterTitle(true);
    cvsADC->cd(4);
    histIC_ADC4vsADC0->Draw("colz");
    histIC_ADC4vsADC0->GetXaxis()->SetTitle("F7IC ADC[0]");
    histIC_ADC4vsADC0->GetXaxis()->SetTitleSize(0.05);
    histIC_ADC4vsADC0->GetXaxis()->SetTitleOffset(0.75);
    histIC_ADC4vsADC0->GetXaxis()->SetLabelSize(0.04);
    histIC_ADC4vsADC0->GetXaxis()->SetTitleFont(22);
    histIC_ADC4vsADC0->GetXaxis()->SetLabelFont(132);
    histIC_ADC4vsADC0->GetXaxis()->CenterTitle(true);
    histIC_ADC4vsADC0->GetYaxis()->SetTitle("F7IC ADC[4]");
    histIC_ADC4vsADC0->GetYaxis()->SetTitleSize(0.05);
    histIC_ADC4vsADC0->GetYaxis()->SetTitleOffset(0.75);
    histIC_ADC4vsADC0->GetYaxis()->SetLabelSize(0.04);
    histIC_ADC4vsADC0->GetYaxis()->SetTitleFont(22);
    histIC_ADC4vsADC0->GetYaxis()->SetLabelFont(132);
    histIC_ADC4vsADC0->GetYaxis()->CenterTitle(true);
    cvsADC->cd(5);
    histIC_ADC5vsADC0->Draw("colz");
    histIC_ADC5vsADC0->GetXaxis()->SetTitle("F7IC ADC[0]");
    histIC_ADC5vsADC0->GetXaxis()->SetTitleSize(0.05);
    histIC_ADC5vsADC0->GetXaxis()->SetTitleOffset(0.75);
    histIC_ADC5vsADC0->GetXaxis()->SetLabelSize(0.04);
    histIC_ADC5vsADC0->GetXaxis()->SetTitleFont(22);
    histIC_ADC5vsADC0->GetXaxis()->SetLabelFont(132);
    histIC_ADC5vsADC0->GetXaxis()->CenterTitle(true);
    histIC_ADC5vsADC0->GetYaxis()->SetTitle("F7IC ADC[5]");
    histIC_ADC5vsADC0->GetYaxis()->SetTitleSize(0.05);
    histIC_ADC5vsADC0->GetYaxis()->SetTitleOffset(0.75);
    histIC_ADC5vsADC0->GetYaxis()->SetLabelSize(0.04);
    histIC_ADC5vsADC0->GetYaxis()->SetTitleFont(22);
    histIC_ADC5vsADC0->GetYaxis()->SetLabelFont(132);
    histIC_ADC5vsADC0->GetYaxis()->CenterTitle(true);
    cvsADC->cd(6);
    histZetBeta->Draw("colz");
    histZetBeta->GetXaxis()->SetTitle("Beta");
    histZetBeta->GetXaxis()->SetTitleSize(0.05);
    histZetBeta->GetXaxis()->SetTitleOffset(0.75);
    histZetBeta->GetXaxis()->SetLabelSize(0.04);
    histZetBeta->GetXaxis()->SetTitleFont(22);
    histZetBeta->GetXaxis()->SetLabelFont(132);
    histZetBeta->GetXaxis()->CenterTitle(true);
    histZetBeta->GetYaxis()->SetTitle("Z");
    histZetBeta->GetYaxis()->SetTitleSize(0.05);
    histZetBeta->GetYaxis()->SetTitleOffset(0.75);
    histZetBeta->GetYaxis()->SetLabelSize(0.04);
    histZetBeta->GetYaxis()->SetTitleFont(22);
    histZetBeta->GetYaxis()->SetLabelFont(132);
    histZetBeta->GetYaxis()->CenterTitle(true);
    cvsADC->Update();

    cvsPID->cd();
    histBeamPID->Draw("colz");
    histBeamPID->GetXaxis()->SetTitle("A/Q");
    histBeamPID->GetXaxis()->SetTitleSize(0.05);
    histBeamPID->GetXaxis()->SetTitleOffset(0.75);
    histBeamPID->GetXaxis()->SetLabelSize(0.04);
    histBeamPID->GetXaxis()->SetTitleFont(22);
    histBeamPID->GetXaxis()->SetLabelFont(132);
    histBeamPID->GetXaxis()->CenterTitle(true);
    histBeamPID->GetYaxis()->SetTitle("Z");
    histBeamPID->GetYaxis()->SetTitleSize(0.05);
    histBeamPID->GetYaxis()->SetTitleOffset(0.75);
    histBeamPID->GetYaxis()->SetLabelSize(0.04);
    histBeamPID->GetYaxis()->SetTitleFont(22);
    histBeamPID->GetYaxis()->SetLabelFont(132);
    histBeamPID->GetYaxis()->CenterTitle(true);
    cvsPID->Update();


    cvsZet->cd();
    cvsZet->Divide(2,1);
    cvsZet->cd(1);
    histZet->Draw();
    cvsZet->cd(2);
    histICMeV->Draw();
    cvsZet->Update();
    cvsZet->SaveAs(Form("./figIC/ICMeV%i.png",ii));    
    histICMeV->Clear();
    histZet->Clear();


    //_________________________________________________________________________________________________________
    //Peak finder in IC Energy (uses the CalibSpectrum class in the folder, compiled with BeamBeam & BeamRaw)
    CalibSpectrum mycal(histICMeV,ICMeVmin,ICMeVmax); // needs the histogram input + the range of the histogram (defined above)
    int nPeaks = 0;
    nPeaks = mycal.FindPeaks(); //finds the number of "Z" peaks seen with TSpectrum
    
    cerr << Form("Found %i candidates for the IC Energy peaks",nPeaks) << endl;
    mycal.FIT(0.04); //half range of the fits (in energy) as input, fits all peaks with gaussians
    
    // in a first step, hinder this part to cross-check which Z the first peak found should be
    //mycal.Calibration(49); //for 132Sn, output on terminal for zcoeff_0 & zcoeff_1
    mycal.Calibration(49); //for 124Sn, output on terminal for zcoeff_0 & zcoeff_1
    
    //    cvsCalib->SaveAs(Form("./figIC/ICcalib%i.png",ii));    


}

