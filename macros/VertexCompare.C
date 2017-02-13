 void VertexCompare(int runNo = 2843){
   auto beam = new TChain("beam");
  auto cbmsim = new TChain("cbmsim");

  for (Int_t i = 0; i < 42; i++) {
  cbmsim -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
  }
  beam -> AddFile(Form("./output/BDC/BDCout.%i.ridf.root",runNo));

  cbmsim -> AddFriend(beam);

  TClonesArray *eventArray = nullptr;
  cbmsim -> SetBranchAddress("STEvent", &eventArray);

  auto *cvs1 = new TCanvas("cvs", "", 800, 800);

  cbmsim -> Draw("z:aoq >> histBeamPID", "", "colz");

/*

   auto c = new TChain("cbmsim");
   for (Int_t i = 0; i < 42; i++) {
   c -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
   }
   //c -> Draw("STVertex.fPos.Y():STVertex.fPos.X() >> hist(300, -100, 100, 300, -330, -130)", "STVertex.fPos.Z()>-18&&STVertex.fPos.Z()<-5", "colz");
*/
 }
