#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class particle;
#pragma link C++ class particleMC;
#pragma link C++ class particleCand;
#pragma link C++ class std::vector<particleCand>;
#pragma link C++ class std::vector<particle>;
#pragma link C++ class std::vector<particleMC>;
#pragma link C++ class utils++;
#pragma link C++ class waveUtils++;
#pragma link C++ class vfempto++;
#pragma link C++ class fempto++;
#pragma link C++ class femptoSource<waveUtils>++;
#pragma link C++ class femptoSource<wignerUtils>++;
#pragma link C++ class vreader++;
#pragma link C++ class readerMC++;
#pragma link C++ class readerFE++;

#pragma link C++ class vrun++;
#pragma link C++ class runSpectrum++;
#pragma link C++ class wignerUtils++;
#pragma link C++ class wignerSource++;

#endif
