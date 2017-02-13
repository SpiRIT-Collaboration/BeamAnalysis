#ifndef CalibSpectrum_h
#define CalibSpectrum_h

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "TF1.h"
#include "TGraph.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1D.h"
#include "TMath.h"
#include "TText.h"
#include "TSpectrum.h"


class CalibSpectrum {

public:
    CalibSpectrum(TH1D *hMeV, double hmin, double hmax);
    virtual ~CalibSpectrum();
    
    Int_t FindPeaks();                          // Find peaks with TSpectrum function
    double fitf(double *x, double *par);        // Definition of the fit function
    void FIT(double fitRange);                  // Create, initialize function, then fit
                                                // pull out the fitted xpeaksfit
    void Calibration(int Zmin);                 // Calibration of Z vs MeV with final peaks
                                                // pull out the zcoeff_0 & zcoeff_1 from it
    
    
protected:
    TCanvas *               cvsCalib;           // Canvas to put the gaussian fits & calib
    TH1D *                  hICMeV;             // Store the MeV histogram there
    double                  hICMeVmin;
    double                  hICMeVmax;
    int                     npeaks;             // # Peaks found with TSpectrum
    std::vector<double>     xpeaks;             // Peaks positions in MeV with TSpectrum
    std::vector<TF1*>       fitFunction;        // Fit ft with npeaks gaussians
    std::vector<double>     xpeaksfit;          // Peaks found in MeV AFTER fit
    std::vector<double>     Zvalues;            // Z values
    
    ClassDef(CalibSpectrum,1);
};

#endif
