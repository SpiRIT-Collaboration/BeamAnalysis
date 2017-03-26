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


using namespace std;




void standaloneICCalib(){
    gROOT->ForceStyle();
    int runNo=2843;
    double goalAoQ;
    if(runNo>=2174 && runNo<=2509) goalAoQ=108./50.;
    if(runNo>=2520 && runNo<=2653) goalAoQ=112./50.;
    if(runNo>=3044 && runNo<=3184) goalAoQ=124./50.;
    if(runNo>=2819 && runNo<=3039) goalAoQ=132./50.;

    double AoQmin = goalAoQ-0.01;
    double AoQmax = goalAoQ+0.01;

    double Zetmin = 49;
    double Zetmax = 51;

    double Betamin = 0.6;
    double Betamax = 0.7;

    double ICMeVmin = 5.4;
    double ICMeVmax = 5.8;

    //_________________________________________________________________________________________________________
    // Reading the root files, and plotting:
    //            PID (Z vs A/Q)                    - to look at the PID before calibration
    //            ADC[i] vs ADC[0] for IC           - Check if pedestal is 0, then no need to change ch2mev_1
    //            Z vs Beta dependence              - For ch2mev_0, if Beta dependence seen, need to correct it
    BeamBeam *beam = new BeamBeam();
    BeamRaw *raw = new BeamRaw();


    //ROOT Files
    //for(int ii=2840; ii<=3039; ii++){  // for 132Sn
    //for(int ii=3058; ii<=3061; ii++){  // for 124Sn
        //int ii = 2840; //132Sn test run
        int ii = runNo; //124Sn test run
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

    auto *histBeamPID = new TH2D("histBeamPID", "histBeamPID", 1500, AoQmin, AoQmax, 1500, Zetmin-1, Zetmax+1);

    auto *histZetBeta = new TH2D("histZetBeta", "histZetBeta", 500, Betamin, Betamax, 500, Zetmin, Zetmax);

    auto *histZet = new TH1D("histZet", "histZet", 500, Zetmin-1, Zetmax+1);
    auto *histICMeV = new TH1D("histICMeV", "histICMeV", 500, ICMeVmin, ICMeVmax);


    //////////////////////////// Entries loop ////////////////////////////
    for(int ientry=0; ientry < nentries; ientry++){
        if(ientry%100000 == 0) cout << "File read: " << ientry*100./nentries << "%" << endl;
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
        if(true && (AoQ>AoQmin && AoQ<AoQmax) && (not(Zet>50.5 && AoQ<goalAoQ+0.5))){
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


    ////////////////in line application of class CalibSpectrum//////////////////
    TCanvas *               cvsCalib;           // Canvas to put the gaussian fits & calib
    TH1D *                  hICMeV;             // Store the MeV histogram there
    double                  hICMeVmin=ICMeVmin;
    double                  hICMeVmax=ICMeVmax;
    int                     npeaks;             // # Peaks found with TSpectrum
    std::vector<double>     xpeaks;             // Peaks positions in MeV with TSpectrum
    std::vector<TF1*>       fitFunction;        // Fit ft with npeaks gaussians
    std::vector<double>     xpeaksfit;          // Peaks found in MeV AFTER fit
    std::vector<double>     Zvalues;

    cvsCalib = new TCanvas("cvsCalib", "cvsCalib", 0, 0, 1000, 600);
    cvsCalib->Divide(2,1);
    cvsCalib->cd(1);
    hICMeV = new TH1D("hICMeV","hICMeV", 500, ICMeVmin, ICMeVmax);
    hICMeV = histICMeV;
    hICMeV->Draw();

    TSpectrum *s = new TSpectrum(Zetmax-Zetmin+1);
    npeaks = s->Search(hICMeV,2,"",0.001);
    double *xp = s->GetPositionX();
    for(int i=0; i<npeaks; i++) xpeaks.push_back(xp[i]);
    sort(xpeaks.begin(), xpeaks.end());

    double fitRange=0.04;
    for(int i = 0; i<npeaks; i++){
        fitFunction.push_back( new TF1(Form("fitFunction_%i",i), "gaus", hICMeVmin, hICMeVmax) );

        fitFunction.at(i)->SetParLimits(0 + 3*i, 0, 2000);
        fitFunction.at(i)->SetParLimits(1 + 3*i, xpeaks[i]-fitRange, xpeaks[i]+fitRange);
        fitFunction.at(i)->SetParLimits(2 + 3*i, 0, 2);
        fitFunction.at(i)->SetParameter(1 + 3*i, xpeaks[i]);

        fitFunction.at(i)->SetLineColor(2);
        fitFunction.at(i)->SetLineStyle(i+1);

        hICMeV->Fit(Form("fitFunction_%i",i),"", "", xpeaks[i]-fitRange, xpeaks[i]+fitRange);

        xpeaksfit.push_back(fitFunction.at(i)->GetParameter(1));

        cout << "Fit results : " << fitFunction.at(i)->GetParameter(1) << endl;
        cout << endl;
    }


    for(int i = 0; i<npeaks; i++){
        cvsCalib->cd(1);
        fitFunction.at(i)->Draw("same");
    }


    int Zmin;
    Zmin=(int)Zetmin;
    for (int i=0; i<npeaks; i++) {
            Zvalues.push_back(Zmin+i);
            cerr << Form("%i -th peak fit: Z = %i",i+1,Zmin+i) << endl;
        }
        TGraph *gcalib = new TGraph(npeaks, &(xpeaksfit[0]),&(Zvalues[0]));
        gcalib->SetMarkerSize(1);
        gcalib->SetMarkerStyle(21);
        TF1 *fpol = new TF1("fpol", "pol1", hICMeVmin, hICMeVmax);
        gcalib->Fit("fpol");

        cvsCalib->cd(2);
        gcalib->Draw("AP");

        cout << "Calibration Fit results :::" << endl;
        cout << "    zcoeff_0 = " << fpol->GetParameter(1) << endl;
        cout << "    zcoeff_1 = "  << fpol->GetParameter(0) << endl;
/*
    //_________________________________________________________________________________________________________
    //Peak finder in IC Energy (uses the CalibSpectrum class in the folder, compiled with BeamBeam & BeamRaw)
    CalibSpectrum mycal(histICMeV,ICMeVmin,ICMeVmax); // needs the histogram input + the range of the histogram (defined above)
    int nPeaks = 0;
    nPeaks = mycal.FindPeaks(); //finds the number of "Z" peaks seen with TSpectrum
    cerr << Form("Found %i candidates for the IC Energy peaks",nPeaks) << endl;
    mycal.FIT(0.04); //half range of the fits (in energy) as input, fits all peaks with gaussians

    // in a first step, hinder this part to cross-check which Z the first peak found should be
    //mycal.Calibration(49); //for 132Sn, output on terminal for zcoeff_0 & zcoeff_1
    mycal.Calibration(47); //for 124Sn, output on terminal for zcoeff_0 & zcoeff_1


*/

}
