void wavetest(float rsource=-2){
  utils::type system=utils::pn;

  float radius = 3.2;
  float V0 = 17.4;
  float spinFact = 3./8;
  float bindingE=-2.22;

  waveUtils a;
  a.setParams(V0, radius, 1.44, rsource,spinFact);
  a.init();

  float mred = 938/2.;

  bool coal = true;

  if(system==utils::pn || system==utils::pp || system==utils::nn) {
    if(system==utils::pn || system==utils::nn){
      a.setCharges(1,0);
      if(system==utils::nn){
        coal = false;
      }
    } else {
      a.setCharges(1,1);
      coal = false;
    }
  }

  if(system==utils::Dn || system==utils::Dp){
    mred = 2*938/3.;
    bindingE = -8.48 + 2.22;
    if(system==utils::Dn){
      a.setCharges(2,0);
    } else {
      a.setCharges(2,1);
    }
  }
  if(system==utils::Tn || system==utils::Tp || system==utils::Hen || system==utils::Hep){
    mred = 3*938/4.;
    bindingE = -28.3 + 8.48;
    if(system==utils::Tn || system==utils::Hen){
      a.setCharges(3,0);
      if(system==utils::Tn){
        coal = false;
      }
    } else if(system==utils::Tp) {
      a.setCharges(3,1);
    } else {
      a.setCharges(3,2);
      coal = false;
    }
  }
  if(system==utils::DD){
    mred = 938;
    bindingE = -28.3 + 2.22*2;
    a.setCharges(4,2);
  }

  float kstar = 50;

  a.setKstar(kstar,1,system);
  a.getUSource()->Draw();
  a.getUDeuteron()->Draw("SAME");

  a.getUDeuteron()->SetLineColor(4);
  printf("kin Source = %f vs %f\n",a.kineticSource(),kstar*kstar/938);
  printf("kin deuteron = %f - V deuteron = %f - E = %f\n",a.kineticDeuteron(),a.potentialDeuteron(),a.kineticDeuteron()+a.potentialDeuteron());
  double err;
  printf("Norm = %f\n",a.getDeuteron2int()->IntegralOneDim(0,20,1E-8,1E-8,err));

  TH1F *hSE = new TH1F("hSE",";k* (Mev/c);",40,0,400);
  TH1F *hME = new TH1F("hME",";k* (Mev/c);",40,0,400);
  hME->SetLineColor(2);

  hSE->Sumw2();
  hME->Sumw2();

  int np=0;
  float x[1000],y[1000],z[1000],t[1000],s[1000],k[1000];
  float coalInt = 0;
  for(kstar = 1; kstar < 1000; kstar+=1){
    x[np] = kstar;
    a.setKstar(kstar,1,system);
    s[np] = a.calcProb() * coal;
    coalInt += s[np]*kstar*kstar;
    y[np] = (a.kineticSource()/mred*(938/2.) + a.potentialSource() - bindingE*s[np])/(1-s[np]);
    k[np] = a.kineticSource()/mred*(938/2.);
    z[np] = kstar*kstar/2/mred;
    t[np] = a.getKstarFinal(s[np], mred, bindingE);
    for(int jj=0; jj < 100; jj++){
      hSE->Fill(t[np],kstar*kstar);
      hME->Fill(x[np],kstar*kstar);
    }

    if(kstar > 199){
//      a.getCoalRe()->Draw();
//      return;
    }
    np++;
  }
  printf("coalInt(Rsource=%.1f)/coalInt(0) = %f\n",std::abs(rsource),coalInt/6052541);

  new TCanvas;
  TGraph *g = new TGraph(np,x,y);
  g->Draw("AP");
  g->SetMarkerStyle(20);
//  TGraph *g2 = new TGraph(np,z,y);
  TGraph *g2 = new TGraph(np,z,y);
  g2->Draw("AP");
  g2->SetMarkerStyle(20);
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  new TCanvas;
  TGraph *g3 = new TGraph(np,x,x);
  g3->Draw("AP");
  g3->SetMarkerStyle(20);
  TGraph *g4 = new TGraph(np,x,t);
  g4->Draw("P");
  g4->SetMarkerStyle(20);
  g4->SetMarkerColor(2);
  g4->SetLineColor(2);
  new TCanvas;
  TGraph *g5 = new TGraph(np,x,s);
  g5->Draw("AP");
  g5->SetMarkerStyle(20);
  new TCanvas;
  hSE->Divide(hME);
  hSE->Draw();

  TFile *fout = new TFile(Form("wave_%.1f.root",rsource),"RECREATE");
  g->Write();
  g2->Write();
  g3->Write();
  g4->Write();
  g5->Write();
  hSE->Write();
  fout->Close();
}
