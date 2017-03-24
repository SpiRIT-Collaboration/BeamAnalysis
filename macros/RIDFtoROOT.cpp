#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "TSystem.h"
#include "TNtuple.h"
#include "TArtEventStore.hh"
#include "TArtStoreManager.hh"
#include "TArtBigRIPSParameters.hh"
#include "TArtSAMURAIParameters.hh"
#include "TArtCalibBDC1Hit.hh"
#include "TArtCalibBDC1Track.hh"
#include "TArtCalibBDC2Hit.hh"
#include "TArtCalibBDC2Track.hh"
#include "TArtDCTrack.hh"
#include "TArtCalibIC.hh"
#include "TArtCalibPPAC.hh"
#include "TArtCalibPlastic.hh"
#include "TArtCalibFocalPlane.hh"
#include "TArtFocalPlane.hh"
#include "TArtEventInfo.hh"
#include "TArtPlastic.hh"
#include "TArtPPAC.hh"
#include "TArtCalibPID.hh"
#include "TArtRecoPID.hh"
#include "TArtRecoRIPS.hh"
#include "TArtRecoTOF.hh"
#include "TArtRecoBeam.hh"
#include "TArtPPACPara.hh"
#include "TArtBeam.hh"
#include "TArtTOF.hh"
#include "TArtIC.hh"
#include "TArtRIPS.hh"
#include "TTree.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TGraph.h"
#include "TClonesArray.h"

using namespace std;

// function to exit loop at keyboard interrupt.
bool stoploop = false;

void stop_interrupt(){
    printf("keyboard interrupt\n");
    stoploop = true;
}

int main(int argc, char *argv[]) {

    /*
     while (1) {
     cout << endl;
     cout << " Is it okay to use 'dc_tpf_" << dctpfRun << ".root' file? (y/n): ";
     char ans;
     cin >> ans;
     if (ans == 'y')
     break;
     else if (ans == 'n') {
     cout << " Stopped!" << endl;
     exit(0);
     }
     }
     */

    float customTOF;
    if(argc < 2){
        cerr << "Missing RIDF file argument" << endl;
    }
    if(argc < 3){
        cerr << "Missing TOF argument, using 280 ns" << endl;
        customTOF=280.;
    }
    else{
        char* charTOF = argv[2];
        customTOF = atof(charTOF);
    }
    char* charRun = argv[1];
    int runNo;
    runNo = atoi(charRun);
    Int_t dctpfRun = runNo;

    cerr << "Run Number =" << runNo << endl;

    // Create StoreManager both for calibration "TArtCalib..." and treatment "TArtReco..."
    //------------------------------------------------------------------------------------
    TArtStoreManager * storeMan = TArtStoreManager::Instance();

    // Create EventStore to control the loop and get the EventInfo
    //------------------------------------------------------------
    TArtEventStore * eventStore = new TArtEventStore();
    eventStore->SetInterrupt(&stoploop);
    TString ridfFile = Form("ridf/SMDAQ%i.ridf", runNo);
    eventStore -> Open(ridfFile.Data());
    //  eventStore -> Open();// for online analysis

    // Create BigRIPSParameters to get Plastics, PPACs, ICs and FocalPlanes parameters from ".xml" files
    //--------------------------------------------------------------------------------------------------
    TArtBigRIPSParameters *bigripsParameters = TArtBigRIPSParameters::Instance();
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSPPAC.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSPlastic.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/BigRIPSIC.xml");
    bigripsParameters -> LoadParameter((Char_t *) "db/FocalPlane.xml");

    // Create SAMURAIParameters
    //-------------------------

    // Create CalibPID to get and calibrate raw data ( CalibPID ->
    //[CalibPPAC , CalibIC, CalibPlastic , CalibFocalPlane] )
    TArtCalibPID *calibPID  = new TArtCalibPID();
    TArtCalibPPAC *calibPPAC = calibPID -> GetCalibPPAC();
    TArtCalibPlastic *calibPLA  = calibPID -> GetCalibPlastic();
    TArtCalibFocalPlane *calibFP   = calibPID -> GetCalibFocalPlane();
    TArtCalibIC *calibIC   = calibPID -> GetCalibIC();

    // Create RecoPID to get calibrated data and to reconstruct TOF, AoQ, Z, ... (RecoPID ->
    //[ RecoTOF , RecoRIPS , RecoBeam] )
    TArtRecoPID *recoPID = new TArtRecoPID();


    // Construct BigRIPS segment
    cout << "Defining bigrips parameters" << endl;
    // Construct transport matrix from the file
    // Set the BigRIPS parameters: up & downstream focal plane, center brho
    // Load the focal plane parameters
    //  TArtRIPS *rips35 = recoPID -> DefineNewRIPS(3, 5, (Char_t *) "matrix/mat1.mat", 7.2929); // F3 - F5
    //  TArtRIPS *rips57 = recoPID -> DefineNewRIPS(5, 7, (Char_t *) "matrix/mat2.mat", 7.3000); // F5 - F7
    TArtRIPS *rips35 = recoPID -> DefineNewRIPS(3, 5, (Char_t *) "matrix/mat1.mat", (Char_t *) "D3"); // F3 - F5
    TArtRIPS *rips57 = recoPID -> DefineNewRIPS(5, 7, (Char_t *) "matrix/mat2.mat", (Char_t *) "D5"); // F5 - F7

    // Construct ToF segment
    // Calculate the flight length from upstream to downstream focal plane.
    // Calculate the flight length of up to middle and middle to downstream focal plane (last argument is middle focal plane)
    // Load the focal plane parameters
    //TArtTOF *tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 242.3 + 90 - 25.8 + 9, 5); // online
    //TArtTOF *tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 280.06, 5); // first offline
    TArtTOF *tof37 = new TArtTOF();
    tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", customTOF, 5);

    //uncomment this section to set TOF for set run intervals
    /*
    if(runNo>=2805 && runNo<=2907) tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 276.9685, 5); // runs 2840-2907
    else if(runNo>=2908 && runNo<=2956) tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 277.0925, 5); // runs 2908-2956
    else if(runNo>=2957 && runNo<=2991) tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 277.204, 5); // runs 2957-2991
    else if(runNo>=2992 && runNo<=3184) tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 277.2805, 5); // runs 2992-3184
    else if(runNo>=3185 && runNo<=3211) tof37  = recoPID -> DefineNewTOF((Char_t *) "F3pl", (Char_t *) "F7pl", 315, 5); // runs 3185-3211(Cocktail)
    else cout << "Problem in run number !!! TOF offset not set !!!" << endl;
    */
    TArtBeam *beam_br37 = recoPID -> DefineNewBeam(rips35, rips57, tof37, (Char_t *) "F7IC");
    TArtBeam *beam_br57 = recoPID -> DefineNewBeam(rips57, tof37, (Char_t *) "F7IC");

    //Output File
    TFile *outfile = new TFile(Form("data/run%d.ridf.root", runNo), "recreate");

    // define data nodes which are supposed to be dumped to tree
    TTree *beam = new TTree("beam", "beam");
    TTree *raw = new TTree("raw", "raw");


    //BigRIPS
    TClonesArray *ppac_array = (TClonesArray *)storeMan->FindDataContainer("BigRIPSPPAC");
    cout<<ppac_array->GetName()<<endl;
    raw->Branch(ppac_array->GetName(),&ppac_array);

    TClonesArray *pla_array = (TClonesArray *)storeMan->FindDataContainer("BigRIPSPlastic");
    cout<<pla_array->GetName()<<endl;
    raw->Branch(pla_array->GetName(),&pla_array);

    TClonesArray * ic_array =
    (TClonesArray *)storeMan->FindDataContainer("BigRIPSIC");
    raw->Branch(ic_array->GetName(),&ic_array);

    TClonesArray * fpl_array =
    (TClonesArray *)storeMan->FindDataContainer("BigRIPSFocalPlane");
    raw->Branch(fpl_array->GetName(),&fpl_array);


    //PID reconstructed data:
    TClonesArray *rips_array =
    (TClonesArray *)storeMan->FindDataContainer("BigRIPSRIPS");
    cout<<rips_array->GetName()<<endl;
    beam->Branch(rips_array->GetName(),&rips_array);

    TClonesArray *tof_array  =
    (TClonesArray *)storeMan->FindDataContainer("BigRIPSTOF");
    cout<<tof_array->GetName()<<endl;
    beam->Branch(tof_array->GetName(),&tof_array);

    TClonesArray *beam_array =
    (TClonesArray *)storeMan->FindDataContainer("BigRIPSBeam");
    cout<<beam_array->GetName()<<endl;
    beam->Branch(beam_array->GetName(),&beam_array);

    //The Focalplane:
    Double_t F3X=-9999; Double_t F3A=-9999; Double_t F3Y=-9999; Double_t F3B=-9999;
    Double_t F5X=-9999; Double_t F5A=-9999; Double_t F5Y=-9999; Double_t F5B=-9999;
    Double_t F7X=-9999; Double_t F7A=-9999; Double_t F7Y=-9999; Double_t F7B=-9999;

    Double_t aoq, zet, tof, beta;


    beam->Branch("F3X",&F3X,"F3X/D");
    beam->Branch("F3A",&F3A,"F3A/D");
    beam->Branch("F3Y",&F3Y,"F3Y/D");
    beam->Branch("F3B",&F3B,"F3B/D");

    beam->Branch("F5X",&F5X,"F5X/D");
    beam->Branch("F5A",&F5A,"F5A/D");
    beam->Branch("F5Y",&F5Y,"F5Y/D");
    beam->Branch("F5B",&F5B,"F5B/D");

    beam->Branch("F7X",&F7X,"F7X/D");
    beam->Branch("F7A",&F7A,"F7A/D");
    beam->Branch("F7Y",&F7Y,"F7Y/D");
    beam->Branch("F7B",&F7B,"F7B/D");


    beam->Branch("aoq",&aoq,"aoq/D");
    beam->Branch("z",&zet,"z/D");
    beam->Branch("tof",&tof,"tof/D");
    beam->Branch("beta",&beta,"beta/D");


    Int_t numEvents = 0;
    while (eventStore -> GetNextEvent()) {
        if (numEvents%1000 == 0) {
            std::cout << "Event number = " << numEvents << std::endl;
        }

        calibPID -> ClearData();
        recoPID -> ClearData();
        calibPID -> ReconstructData();
        recoPID -> ReconstructData();

        aoq = beam_br37 -> GetAoQ();
        zet = beam_br37 -> GetZet();
        tof = tof37 -> GetTOF();
        beta = tof37 -> GetBeta();



        F3X = -9999; F3A = -9999; F3Y = -9999; F3B = -9999;
        F5X = -9999; F5A = -9999; F5Y = -9999; F5B = -9999;
        F7X = -9999; F7A = -9999; F7Y = -9999; F7B = -9999;




        TArtFocalPlane *fp = nullptr;
        TVectorD *vec;
        fp = calibFP -> FindFocalPlane(3);
        if (fp) {
            vec = fp->GetOptVector();
            F3X = (*vec)(0); F3A = (*vec)(1);
            F3Y = (*vec)(2); F3B = (*vec)(3);
        }
        fp = calibFP -> FindFocalPlane(5);
        if (fp) {
            vec = fp -> GetOptVector();
            F5X = (*vec)(0); F5A = (*vec)(1);
            F5Y = (*vec)(2); F5B = (*vec)(3);
        }
        fp = calibFP -> FindFocalPlane(7);
        if (fp) {
            vec = fp -> GetOptVector();
            F7X = (*vec)(0); F7A = (*vec)(1);
            F7Y = (*vec)(2); F7B = (*vec)(3);
        }

        beam -> Fill();
        raw->Fill();

        numEvents++;

    }

    outfile -> Write();

    //gSystem -> Exec(Form("ln -sf ../data/run%d.ridf.root plaData_experiment132Sn/run%d.ridf.root", runNo, runNo));
}
