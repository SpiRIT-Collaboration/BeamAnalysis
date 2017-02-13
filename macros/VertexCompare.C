 void VertexCompare(int runNo = 2843){
   auto TGTmag = new TChain("TGTmag");
  auto cbmsim = new TChain("cbmsim");

  for (Int_t i = 0; i < 42; i++) {
  cbmsim -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
  }
  TGTmag -> AddFile(Form("./output/BDC/BDCout.%i.ridf.root",runNo));

  cbmsim -> AddFriend(TGTMag);

  TClonesArray *eventArray = nullptr;
  cbmsim -> SetBranchAddress("STVertex", &eventArray);

  auto *cvs1 = new TCanvas("cvs", "", 800, 800);
  cvs1->Divide(2,2);
  cvs1->cd(1);
  cbmsim -> Draw("TGT_y_0_5T:TGT_x_0_5T >> hist(300, -100, 100, 300, -330, -130)", "", "colz");
  cvs1->cd(2);
  cbmsim -> Draw("STVertex.fPos.Y():STVertex.fPos.X() >> hist(300, -100, 100, 300, -330, -130)", "STVertex.fPos.Z()>-18&&STVertex.fPos.Z()<-5", "colz");
/*

   auto c = new TChain("cbmsim");
   for (Int_t i = 0; i < 42; i++) {
   c -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
   }
   //c -> Draw("STVertex.fPos.Y():STVertex.fPos.X() >> hist(300, -100, 100, 300, -330, -130)", "STVertex.fPos.Z()>-18&&STVertex.fPos.Z()<-5", "colz");
*/
 }
