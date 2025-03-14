int pdgPr = 2212; // proton
int pdgNe = 2112; // neutron
int pdgDe = 4324; // deuteron

void run(){
  // set file reader
  readerMC *reader = new readerMC();
  reader->openFile("listaMC",true); // true to read a collection
  int nev = reader->getNevents();
  printf("N events found = %d\n",nev);

  // set analisys
  int nmixedEvents = 10;
  vrun *run = new vrun(nmixedEvents);
  run->selectPDG(pdgPr,pdgPr);       // pdg of the two particles to be correlated
  run->init(); // init analysys

  for(int i=0; i < nev; i++){
    // read event
    reader->NextEvent();
    std::vector<particleMC>& event = reader->getParticles();

    // for MC add it the effect of afterburner

    //---------------------------------------

    // then analyze
    run->setEvent(event); // pass the event to the analyzer
    run->process();       // process the event
  }

  // finalization
  run->finalize();  // finalize the analysis
  TFile *fout = new TFile("runOutput.root","RECREATE");
  run->write();     // write outputs of the analysis
  fout->Close();
}
