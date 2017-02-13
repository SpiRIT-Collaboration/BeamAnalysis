#include "CalibSpectrum.h"
#include "TH1D.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TF1.h"
#include <vector>
#include <iostream>

using namespace std;


ClassImp(CalibSpectrum);

//__________________________________________________________
CalibSpectrum::CalibSpectrum(TH1D *hMeV, double hmin, double hmax){
    cvsCalib = new TCanvas("cvsCalib", "cvsCalib", 0, 0, 1000, 600);
    cvsCalib->Divide(2,1);
    cvsCalib->cd(1);
    hICMeV = new TH1D("hICMeV","hICMeV", 500, hmin, hmax);
    
    hICMeV = hMeV;
    hICMeVmin = hmin;
    hICMeVmax = hmax;
    
    hICMeV->Draw();

}
CalibSpectrum::~CalibSpectrum(void){}

//__________________________________________________________
double CalibSpectrum::fitf(double *x, double *par){
    double fitval = 0.;
    int ij = 0;
    for(int i=0; i<npeaks; i++){
        ij = 3*i;
        fitval+= par[ij+0] * exp( -0.5*((x[0]-par[ij+1])/par[ij+2])*((x[0]-par[ij+1])/par[ij+2]) );
    }
    
    return fitval;
}



//__________________________________________________________
Int_t CalibSpectrum::FindPeaks(void){

    TSpectrum *s = new TSpectrum(10);
    npeaks = s->Search(hICMeV,2,"",0.05);
    double *xp = s->GetPositionX();
    for(int i=0; i<npeaks; i++) xpeaks.push_back(xp[i]);
    sort(xpeaks.begin(), xpeaks.end());     // need to sort the peaks founds by increasing energy
    return npeaks;
}

//__________________________________________________________
void CalibSpectrum::FIT(double fitRange){
    /*
    fitFunction = new TF1("fitFunction", CalibSpectrum::fitf, hICMeVmin, hICMeVmax, 3*npeaks);
    for (int i=0; i<npeaks; i++) {
        fitFunction->SetParLimits(0 + 3*i, 0, 2000);
        fitFunction->SetParLimits(1 + 3*i, xpeaks[i]-0.05, xpeaks[i]+0.05);
        fitFunction->SetParLimits(2 + 3*i, 0, 2);
        
        fitFunction->SetParameter(1 + 3*i, xpeaks[i]);
    }
     */
    
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
    
    
}

//__________________________________________________________
void CalibSpectrum::Calibration(int Zmin){
    
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

}

