 void VertexCompare(int runNo = 2843){
   auto TGTmag = new TChain("TGTmag");
  auto cbmsim = new TChain("cbmsim");

  for (Int_t i = 0; i < 42; i++) {
  cbmsim -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
  }
  TGTmag -> AddFile(Form("./output/BDC/BDCout.%i.ridf.root",runNo));

  cbmsim -> AddFriend(TGTmag);

  TClonesArray *vertexArray = nullptr;
  cbmsim -> SetBranchAddress("STVertex", &vertexArray);

  Double_t bdc_x,bdc_y;
  cbmsim->SetBranchAddress("TGT_x_0_5T",&bdc_x);
  cbmsim->SetBranchAddress("TGT_y_0_5T",&bdc_y);
  Double_t bdc_px,bdc_py;
  cbmsim->SetBranchAddress("TGT_px_0_5T",&bdc_px);
  cbmsim->SetBranchAddress("TGT_py_0_5T",&bdc_py);
  
  
  auto *histXdif= new TH1D("histXdif","TPC vertex X - BDC vertex X",200,-10,10);
  auto *histYdif= new TH1D("histYdif","TPC vertex Y - BDC vertex Y",200,-10,10);
  auto *histXYprofile = new TH2D("XYprofile","beam profile reconstructed from TPC",200,-50,50,200,-50,50);
  auto *histpx= new TH1D("histpx","momentum px (MeV/c)",200,-1000,1000);
  auto *histpy= new TH1D("histpy","momentum py (MeV/c)",200,-1000,1000);
  
  
  Int_t nEvents = cbmsim -> GetEntries();
  cout << "Number of events: " << nEvents << endl;

  for(int iEvent=0; iEvent<nEvents; iEvent++){
    cbmsim->GetEntry(iEvent);
    Int_t nVertices = vertexArray -> GetEntries();
    if (nVertices != 1)
      continue;

    STVertex *vertex = (STVertex *) vertexArray -> At(0);
    TVector3 posVertex = vertex -> GetPos();
    //cout << posVertex.Y()/10. +23<<"," << bdc_y<< endl;
    if( posVertex.X()-bdc_x>4 /*|| posVertex.X()-bdc_x>4  || posVertex.Y()-bdc_y+227<-2 || posVertex.Y()-bdc_y+227>2*/ ){
      histXdif->Fill(posVertex.X()-bdc_x);
      histYdif->Fill(posVertex.Y()+227-bdc_y);
      histXYprofile->Fill(bdc_x,bdc_y);
      histpx->Fill(bdc_px);
      histpy->Fill(bdc_py);
    }
  }

  
  auto *cvs1 = new TCanvas("cvs", "", 800, 800);
  cvs1->Divide(2,2);
  cvs1->cd(1);
  histXYprofile->Draw("colz");
  //cbmsim -> Draw("bdc_y:bdc_x >> hist(300, -100, 100, 300, -100, 100)", "", "colz");
  cvs1->cd(2);
  cbmsim -> Draw("STVertex.fPos.Y():STVertex.fPos.X() >> hist(300, -100, 100, 300, -330, -130)", "STVertex.fPos.Z()>-18&&STVertex.fPos.Z()<-5", "colz");
  cvs1->cd(3);
  histXdif->Draw();
  cvs1->cd(4);
  histYdif->Draw();
  

 }
