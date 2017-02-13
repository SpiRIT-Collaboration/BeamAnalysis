 void VertexCompare(int runNo = 2843){

   auto c = new TChain("cbmsim");
   for (Int_t i = 0; i < 42; i++) {
   c -> AddFile(Form("/mnt/spirit/analysis/changj/SpiRITROOT/macros/data/run2843_s%d.reco.v1.04.root", i));
   }
   c -> Draw("STVertex.fPos.Y():STVertex.fPos.X() >> hist(300, -100, 100, 300, -330, -130)", "STVertex.fPos.Z()>-18&&STVertex.fPos.Z()<-5", "colz");
 }
