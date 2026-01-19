int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron
int pdgAntiPr = -2212; // antiproton
int pdgAntiNe = -2112; // antineutron
int pdgAntiDe = -4324; // antideuteron

bool isWigner = false;

#define MAXD(...) std::max(std::initializer_list<double>{__VA_ARGS__})

void runAll(TString listname, TString outfolder, TString outname, double R_D = 3.2, double R_0 = 1.125)
{
  // set file reader
  vreader *reader = new readerFE();
  //reader->setEtaRange(-100.,100);
  reader->openFile(listname, true); // true to read a collection
  int nev = reader->getNevents();
  printf("N events found = %d\n", nev);

  // set analisys
  int nmixedEvents = 10;
  vrun *runPrPr = new vrun(nmixedEvents);
  runPrPr->setIDName("PrPr");
  runPrPr->selectPDG(pdgPr, pdgPr); // pdg of the two particles to be correlated
  runPrPr->init();                  // init analysis

  vrun *runPrNe = new vrun(nmixedEvents);
  runPrNe->setIDName("PrNe");
  runPrNe->selectPDG(pdgPr, pdgNe); // pdg of the two particles to be correlated
  runPrNe->init();                  // init analysis

  vrun *runPrDe = new vrun(nmixedEvents);
  runPrDe->setIDName("PrDe");
  runPrDe->selectPDG(pdgPr, pdgDe); // pdg of the two particles to be correlated
  runPrDe->init();                  // init analysis

  vrun *runAntiPrPr = new vrun(nmixedEvents);
  runAntiPrPr->setIDName("AntiPrPr");
  runAntiPrPr->selectPDG(pdgAntiPr, pdgAntiPr); // pdg of the two particles to be correlated
  runAntiPrPr->init();                  // init analysis

  vrun *runAntiPrNe = new vrun(nmixedEvents);
  runAntiPrNe->setIDName("AntiPrNe");
  runAntiPrNe->selectPDG(pdgAntiPr, pdgAntiNe); // pdg of the two particles to be correlated
  runAntiPrNe->init();                  // init analysis

  vrun *runAntiPrDe = new vrun(nmixedEvents);
  runAntiPrDe->setIDName("AntiPrDe");
  runAntiPrDe->selectPDG(pdgAntiPr, pdgAntiDe); // pdg of the two particles to be correlated
  runAntiPrDe->init();                  // init analysis

  vrun *runPrAntiPr = new vrun(nmixedEvents);
  runPrAntiPr->setIDName("PrAntiPr");
  runPrAntiPr->selectPDG(pdgPr, pdgAntiPr); // pdg of the two particles to be correlated
  runPrAntiPr->init();                  // init analysis

  vrun *runPrAntiDe = new vrun(nmixedEvents);
  runPrAntiDe->setIDName("PrAntiDe");
  runPrAntiDe->selectPDG(pdgPr, pdgAntiDe); // pdg of the two particles to be correlated
  runPrAntiDe->init();                  // init analysis

  vrun *runDeAntiPr = new vrun(nmixedEvents);
  runDeAntiPr->setIDName("DeAntiPr");
  runDeAntiPr->selectPDG(pdgDe, pdgAntiPr); // pdg of the two particles to be correlated
  runDeAntiPr->init();                  // init analysis

  runSpectrum *spectra = new runSpectrum(); //spectrum analysis
  //spectra->setRapidityRange(-100.,100);
  spectra->init();

  //waveUtils::setParams(10, R);

  //waveUtils::setParams(waveUtils::calculateV0(), R);

  //-------- interactor ------------------------
  double vStrong = 0;
  vfempto *interactor;
  auto R_0_eff = R_0*2*(-1); 
  if(! isWigner){
   interactor = new femptoSource;      // Lenard-Jones  strong potential
   interactor->setParams(17.4, R_D, 1.44, R_0_eff, 3./8);
   vStrong = waveUtils::calculateV0();
   interactor->setParams(vStrong, R_D, 1.44, R_0_eff, 3./8);
  } else {
   interactor = new femptoSource<wignerUtils>;      // Lenard-Jones trong potential
   interactor->setParams(17.4E-3, 3.2, 1.44E-3, 1.5, 3./8);
  }
  interactor->setThreshold(0.4);

  //-------- interactor ------------------------
  TFile *weightFilePr = new TFile("weight/weightPr.root", "READ");
  auto hWeightTotalPr = (TH1D *)weightFilePr->Get("ratioPr");

  std::cout<<"Strong Radius: "<< R_D << "\t Strong Potential: "<< vStrong <<"\n";

  for (int i = 0; i < nev; i++)
  {
    // read event
    reader->NextEvent();
    std::vector<particleMC> &event = reader->getParticles();
    auto PrPt = returnPt(event, -0.5, 0.5, pdgPr);
    auto NePt = returnPt(event, -0.5, 0.5, pdgNe);
    auto AntiPrPt = returnPt(event, -0.5, 0.5, pdgAntiPr);
    auto AntiNePt = returnPt(event, -0.5, 0.5, pdgAntiNe);
    auto ptValue = MAXD(PrPt, NePt, AntiPrPt, AntiNePt);
    double weight;
    if (ptValue < 0) weight = 1;
    else weight = computeWeight(ptValue, hWeightTotalPr);

    if (i % 100000 == 0)
      std::cout << "Event: " << i << "/" << nev << " - " << ((double)i / (double)nev) * 100 << "\n";

    // for MC add it the effect of afterburner
    interactor->doInteractAll(event, false, true);
    //---------------------------------------

    // then analyze
    runPrPr->setEvent(event); // pass the event to the analyzer
    runPrPr->doAnalysis();       // process the event

    runPrNe->setEvent(event); // pass the event to the analyzer
    runPrNe->doAnalysis();       // process the event

    runPrDe->setEvent(event); // pass the event to the analyzer
    runPrDe->doAnalysis();       // process the event

    runAntiPrPr->setEvent(event); // pass the event to the analyzer
    runAntiPrPr->doAnalysis();       // process the event

    runAntiPrNe->setEvent(event); // pass the event to the analyzer
    runAntiPrNe->doAnalysis();       // process the event

    runAntiPrDe->setEvent(event); // pass the event to the analyzer
    runAntiPrDe->doAnalysis();       // process the event

    runPrAntiPr->setEvent(event); // pass the event to the analyzer
    runPrAntiPr->doAnalysis();       // process the event

    runPrAntiDe->setEvent(event); // pass the event to the analyzer
    runPrAntiDe->doAnalysis();       // process the event

    runDeAntiPr->setEvent(event); // pass the event to the analyzer
    runDeAntiPr->doAnalysis();       // process the event

    spectra->setEvent(event);
    spectra->doAnalysis(weight);
  }

  // finalization
  runPrPr->finalize(); // finalize the analysis
  TFile *foutPrPr = new TFile(outfolder + "pp_" + outname, "RECREATE");
  runPrPr->write(); // write outputs of the analysis
  //hCorrelationPrPr->Write();
  foutPrPr->Close();

  // finalization
  runPrNe->finalize(); // finalize the analysis
  TFile *foutPrNe = new TFile(outfolder + "pn_" + outname, "RECREATE");
  runPrNe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrNe->Close();

  runPrDe->finalize(); // finalize the analysis
  TFile *foutPrDe = new TFile(outfolder + "pd_" + outname, "RECREATE");
  runPrDe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrDe->Close();

  runAntiPrPr->finalize(); // finalize the analysis
  TFile *foutAntiPrPr = new TFile(outfolder + "Antipp_" + outname, "RECREATE");
  runAntiPrPr->write(); // write outputs of the analysis
  //hCorrelationPrPr->Write();
  foutAntiPrPr->Close();

  // finalization
  runAntiPrNe->finalize(); // finalize the analysis
  TFile *foutAntiPrNe = new TFile(outfolder + "Antipn_" + outname, "RECREATE");
  runAntiPrNe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutAntiPrNe->Close();

  runAntiPrDe->finalize(); // finalize the analysis
  TFile *foutAntiPrDe = new TFile(outfolder + "Antipd_" + outname, "RECREATE");
  runAntiPrDe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutAntiPrDe->Close();

  runPrAntiPr->finalize(); // finalize the analysis
  TFile *foutPrAntiPr = new TFile(outfolder + "pAntip_" + outname, "RECREATE");
  runPrAntiPr->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrAntiPr->Close();

  runPrAntiDe->finalize(); // finalize the analysis
  TFile *foutPrAntiDe = new TFile(outfolder + "pAntid_" + outname, "RECREATE");
  runPrAntiDe->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutPrAntiDe->Close();

  runDeAntiPr->finalize(); // finalize the analysis
  TFile *foutDeAntiPr = new TFile(outfolder + "dAntip_" + outname, "RECREATE");
  runDeAntiPr->write(); // write outputs of the analysis
  //hCorrelationPrNe->Write();
  foutDeAntiPr->Close();

  TFile *foutSpec = new TFile(outfolder + "Spectrum" + outname , "RECREATE");
  spectra->write();
  foutSpec->Close();

}
