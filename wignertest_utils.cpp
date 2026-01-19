void wignertest_utils()
{
    wignerUtils::setParams();
    wignerUtils::init();
    wignerUtils::setIntegrationRanges(0., 20, 0., 0.6);
    wignerUtils::setSourceRadius(1.0); // Equivalent to setR0(1.0)

    TString s = "data.root";
    TFile *file = new TFile(s, "RECREATE");
    std::cout << "Creating " << s << "\n";

    TTree *tree = new TTree("tree", "W x W");

    double r0, WW, norm, wH, wK, wV, k;
    tree->Branch("r0", &r0, "r0/D");
    tree->Branch("WxW", &WW, "WW/D");
    tree->Branch("norm", &norm, "norm/D");
    tree->Branch("wH", &wH, "wH/D");
    tree->Branch("wK", &wK, "wK/D");
    tree->Branch("wV", &wV, "wV/D");
    tree->Branch("k", &k, "k/D");

    std::cout << "Deuteron integral: " << wignerUtils::getDeuteronInt() << "\n";

    float xk[1000], xP[1000], xK[1000], xE[1000], xV[1000], xkf[1000];
    int np = 0;
    double mRed = utils::getMass(2212) / 2;

    for (double i = 1E-3 / 2; i < 0.5; i += 2E-3)
    {
        xk[np] = i;
        wignerUtils::setRadiusK(i);
        r0 = wignerUtils::getRadius();
        k = i;
        norm = wignerUtils::getNorm();
        WW = wignerUtils::checkWxW();
        wK = wignerUtils::kineticSource();
        wV = wignerUtils::potentialSource();
        wH = wignerUtils::getWH();

        xK[np] = i * i * 0.5 / mRed;
        xV[np] = wV;
        xP[np] = wignerUtils::getcoal();

        if (xK[np] > xV[np])
        {
            xE[np] = xK[np] + xV[np];
            xkf[np] = wignerUtils::getKstarFinal(xP[np], mRed); //sqrt(xE[np] * 2 * mRed);
        }
        else
        {
            xE[np] = 0;
            xkf[np] = 0;
        }

        if (WW < 0.00001)
            WW = 1;

        std::cout << "i : " << i << " coal: " << xP[np] << " r0:  " << r0 << " k*: " << k
                  << " Norm: " << norm << " Check: " << WW << " K: " << wK
                  << " V: " << wV << " H: " << wH << "\n";

        tree->Fill();
        ++np;
    }

    tree->Write();
    file->Close();

    int clean_np = 0;
    double xk_clean[1000], xP_clean[1000], xV_clean[1000];

    for (int i = 0; i < np; ++i)
    {
        if (xkf[i] == 0 || xP[i] == 0)
            continue;
        xk_clean[clean_np] = xk[i];
        xP_clean[clean_np] = xP[i];
        xV_clean[clean_np] = xV[i];
        clean_np++;
    }

    TCanvas *c = new TCanvas("wignertest_utils", "wignertest_utils");
    c->Divide(2, 2);
    c->cd(1);
    TGraph *g = new TGraph(clean_np, xk_clean, xP_clean);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle("Coal prob R_{deuteron} = 3.2 fm;k*_{0} (GeV/c);P(k*_{0})");
    g->SetMinimum(0);

    c->cd(2);
    g = new TGraph(np, xk, xkf);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle("Effect of interaction;k*_{0} (GeV/c);k*_{f} (GeV/c)");
    TLine *l = new TLine(0, 0, 0.5, 0.5);
    l->SetLineColor(2);
    l->Draw("SAME");

    c->cd(3);
    g = new TGraph(np, xK, xE);
    g->SetMarkerStyle(20);
    g->SetTitle("Effect of interaction;K.E._{0} (GeV);E_{f} (GeV)");
    g->Draw("AP");
    l = new TLine(0, 0, 0.25, 0.25);
    l->SetLineColor(2);
    l->Draw("SAME");

    c->cd(4);
    g = new TGraph(clean_np, xk_clean, xV_clean);
    g->SetMarkerStyle(20);
    g->Draw("AP");
    g->SetTitle("V interaction;k*_{0} (GeV/c);V (GeV)");
}
