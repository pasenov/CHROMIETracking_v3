{

 TFile *f = new TFile("PixelTelescope_BeamData_DQM.root");
 f->cd("DQMData/run100000"); 

 TCanvas *c = new TCanvas("c","X Coordinate Layer 0",1200,1000);
 c->Divide(4,2);

 c->cd(1);
 cluster3DTree->Draw("x","detId==353375236");
 
 c->cd(2);
 cluster3DTree->Draw("x","detId==353113092");

 c->cd(3);
 cluster3DTree->Draw("x","detId==352850948");

 c->cd(4);
 cluster3DTree->Draw("x","detId==352588804");

 c->cd(5);
 cluster3DTree->Draw("x","detId==344201220");

 c->cd(6);
 cluster3DTree->Draw("x","detId==344463364");

 c->cd(7);
 cluster3DTree->Draw("x","detId==344725508");

 c->cd(8);
 cluster3DTree->Draw("x","detId==344987652");

// 

 TCanvas *c1 = new TCanvas("c1","X Coordinate Layer 1",1200,1000);
 c1->Divide(4,2);

 c1->cd(1);
 cluster3DTree->Draw("x","detId==353376260");
 
 c1->cd(2);
 cluster3DTree->Draw("x","detId==353114116");

 c1->cd(3);
 cluster3DTree->Draw("x","detId==352851972");

 c1->cd(4);
 cluster3DTree->Draw("x","detId==352589828");

 c1->cd(5);
 cluster3DTree->Draw("x","detId==344200196");

 c1->cd(6);
 cluster3DTree->Draw("x","detId==344462340");

 c1->cd(7);
 cluster3DTree->Draw("x","detId==344724484");

 c1->cd(8);
 cluster3DTree->Draw("x","detId==344986628");
}
