{

  TH2D *hxz = new TH2D("hxz","hxz; X(mm); Z(mm)",1000,-300,300,1000,-300,300);
  TH2D *hyz = new TH2D("hyz","hyz; Y(mm); Z(mm)",1000,-300,300,1000,-300,300);
  TBox *b1 = new TBox(-150,-150,150,150);
  b1->SetLineColor(1);
  b1->SetFillStyle(0);
  TCanvas *c1 = new TCanvas("c1","c1",1400,600);
  c1->Divide(2,1);
  trk->SetMarkerStyle(20);
  trk->SetMarkerSize(0.2);
  for (int iev=0; iev<100; iev++){
    c1->cd(1);
    hxz->Draw();
    b1->Draw();
    trk->SetMarkerColor(1);
    trk->Draw("z:-x",Form("e>0.005&&mat==1&&event==%i",iev),"same");
    trk->SetMarkerColor(2);
    trk->Draw("z:-x",Form("e>0.005&&mat==2&&event==%i",iev),"same");
    c1->cd(2);
    hyz->Draw();
    b1->Draw();
    trk->SetMarkerColor(1);
    trk->Draw("z:y",Form("e>0.005&&mat==1&&event==%i",iev),"same");
    trk->SetMarkerColor(2);
    trk->Draw("z:y",Form("e>0.005&&mat==2&&event==%i",iev),"same");
    c1->Update();
    c1->Print(Form("../output/plot/plot%i.png",iev));
  }
  
}
