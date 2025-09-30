#ifndef CRUNSPECTRUM
#define CRUNSPECTRUM


#include "Cvrun.h"
#include "TH1D.h"

class runSpectrum : public vrun
{
    public:
      runSpectrum(int nmix=10, TString dirName = "runSpectrum") : vrun(nmix, dirName){}
      void setRapidityRange(float minRapidity, float maxRapidity) {mMinRapidity = minRapidity, mMaxRapidity = maxRapidity;}
    private:
      TH1D *mHProtonSpectrum = nullptr;
      TH1D *mHNeutronSpectrum = nullptr;
      TH1D *mHDeuteronSpectrum = nullptr;
      TH1D *mHTritiumSpectrum = nullptr;
      TH1D *mHHelium3Spectrum = nullptr;
      TH1D *mHHelium4Spectrum = nullptr;

      TH1D *mHProtonSpectrumA = nullptr;
      TH1D *mHNeutronSpectrumA = nullptr;
      TH1D *mHDeuteronSpectrumA = nullptr;
      TH1D *mHTritiumSpectrumA = nullptr;
      TH1D *mHHelium3SpectrumA = nullptr;
      TH1D *mHHelium4SpectrumA = nullptr;

      TH1D *mHProtonSpectrumM = nullptr;
      TH1D *mHNeutronSpectrumM = nullptr;
      TH1D *mHDeuteronSpectrumM = nullptr;
      TH1D *mHTritiumSpectrumM = nullptr;
      TH1D *mHHelium3SpectrumM = nullptr;
      TH1D *mHHelium4SpectrumM = nullptr;

      TH2D *mHProtonSpectrumARapidity = nullptr;
      TH2D *mHNeutronSpectrumARapidity = nullptr;
      TH2D *mHDeuteronSpectrumARapidity = nullptr;
      TH2D *mHTritiumSpectrumARapidity = nullptr;
      TH2D *mHHelium3SpectrumARapidity = nullptr;
      TH2D *mHHelium4SpectrumARapidity = nullptr;

      TH2D *mHProtonSpectrumMRapidity = nullptr;
      TH2D *mHNeutronSpectrumMRapidity = nullptr;
      TH2D *mHDeuteronSpectrumMRapidity = nullptr;
      TH2D *mHTritiumSpectrumMRapidity = nullptr;
      TH2D *mHHelium3SpectrumMRapidity = nullptr;
      TH2D *mHHelium4SpectrumMRapidity = nullptr;

      TH2D *mHProtonSpectrumAPEta = nullptr;
      TH2D *mHNeutronSpectrumAPEta = nullptr;
      TH2D *mHDeuteronSpectrumAPEta = nullptr;
      TH2D *mHTritiumSpectrumAPEta = nullptr;
      TH2D *mHHelium3SpectrumAPEta = nullptr;
      TH2D *mHHelium4SpectrumAPEta = nullptr;

      TH2D *mHProtonSpectrumMPEta = nullptr;
      TH2D *mHNeutronSpectrumMPEta = nullptr;
      TH2D *mHDeuteronSpectrumMPEta = nullptr;
      TH2D *mHTritiumSpectrumMPEta = nullptr;
      TH2D *mHHelium3SpectrumMPEta = nullptr;
      TH2D *mHHelium4SpectrumMPEta = nullptr;

      TH2D *mHProtonSpectrumAPtEta = nullptr;
      TH2D *mHNeutronSpectrumAPtEta = nullptr;
      TH2D *mHDeuteronSpectrumAPtEta = nullptr;
      TH2D *mHTritiumSpectrumAPtEta = nullptr;
      TH2D *mHHelium3SpectrumAPtEta = nullptr;
      TH2D *mHHelium4SpectrumAPtEta = nullptr;

      TH2D *mHProtonSpectrumMPtEta = nullptr;
      TH2D *mHNeutronSpectrumMPtEta = nullptr;
      TH2D *mHDeuteronSpectrumMPtEta = nullptr;
      TH2D *mHTritiumSpectrumMPtEta = nullptr;
      TH2D *mHHelium3SpectrumMPtEta = nullptr;
      TH2D *mHHelium4SpectrumMPtEta = nullptr;

      TH1D *mHAntiProtonSpectrum = nullptr;
      TH1D *mHAntiNeutronSpectrum = nullptr;
      TH1D *mHAntiDeuteronSpectrum = nullptr;
      TH1D *mHAntiTritiumSpectrum = nullptr;
      TH1D *mHAntiHelium3Spectrum = nullptr;
      TH1D *mHAntiHelium4Spectrum = nullptr;
      
      TH1D *mHAntiProtonSpectrumA = nullptr;
      TH1D *mHAntiNeutronSpectrumA = nullptr;
      TH1D *mHAntiDeuteronSpectrumA = nullptr;
      TH1D *mHAntiTritiumSpectrumA = nullptr;
      TH1D *mHAntiHelium3SpectrumA = nullptr;
      TH1D *mHAntiHelium4SpectrumA = nullptr;

      TH1D *mHAntiProtonSpectrumM = nullptr;
      TH1D *mHAntiNeutronSpectrumM = nullptr;
      TH1D *mHAntiDeuteronSpectrumM = nullptr;
      TH1D *mHAntiTritiumSpectrumM = nullptr;
      TH1D *mHAntiHelium3SpectrumM = nullptr;
      TH1D *mHAntiHelium4SpectrumM = nullptr;

      TH2D *mHAntiProtonSpectrumARapidity = nullptr;
      TH2D *mHAntiNeutronSpectrumARapidity = nullptr;
      TH2D *mHAntiDeuteronSpectrumARapidity = nullptr;
      TH2D *mHAntiTritiumSpectrumARapidity = nullptr;
      TH2D *mHAntiHelium3SpectrumARapidity = nullptr;
      TH2D *mHAntiHelium4SpectrumARapidity = nullptr;

      TH2D *mHAntiProtonSpectrumMRapidity = nullptr;
      TH2D *mHAntiNeutronSpectrumMRapidity = nullptr;
      TH2D *mHAntiDeuteronSpectrumMRapidity = nullptr;
      TH2D *mHAntiTritiumSpectrumMRapidity = nullptr;
      TH2D *mHAntiHelium3SpectrumMRapidity = nullptr;
      TH2D *mHAntiHelium4SpectrumMRapidity = nullptr;

      TH2D *mHAntiProtonSpectrumAPEta = nullptr;
      TH2D *mHAntiNeutronSpectrumAPEta = nullptr;
      TH2D *mHAntiDeuteronSpectrumAPEta = nullptr;
      TH2D *mHAntiTritiumSpectrumAPEta = nullptr;
      TH2D *mHAntiHelium3SpectrumAPEta = nullptr;
      TH2D *mHAntiHelium4SpectrumAPEta = nullptr;

      TH2D *mHAntiProtonSpectrumMPEta = nullptr;
      TH2D *mHAntiNeutronSpectrumMPEta = nullptr;
      TH2D *mHAntiDeuteronSpectrumMPEta = nullptr;
      TH2D *mHAntiTritiumSpectrumMPEta = nullptr;
      TH2D *mHAntiHelium3SpectrumMPEta = nullptr;
      TH2D *mHAntiHelium4SpectrumMPEta = nullptr;

      TH2D *mHAntiProtonSpectrumAPtEta = nullptr;
      TH2D *mHAntiNeutronSpectrumAPtEta = nullptr;
      TH2D *mHAntiDeuteronSpectrumAPtEta = nullptr;
      TH2D *mHAntiTritiumSpectrumAPtEta = nullptr;
      TH2D *mHAntiHelium3SpectrumAPtEta = nullptr;
      TH2D *mHAntiHelium4SpectrumAPtEta = nullptr;

      TH2D *mHAntiProtonSpectrumMPtEta = nullptr;
      TH2D *mHAntiNeutronSpectrumMPtEta = nullptr;
      TH2D *mHAntiDeuteronSpectrumMPtEta = nullptr;
      TH2D *mHAntiTritiumSpectrumMPtEta = nullptr;
      TH2D *mHAntiHelium3SpectrumMPtEta = nullptr;
      TH2D *mHAntiHelium4SpectrumMPtEta = nullptr;

      TH2D *mHProtonPEta = nullptr;
      TH2D *mHNeutronPEta = nullptr;
      TH2D *mHDeuteronPEta = nullptr;
      TH2D *mHTritiumPEta = nullptr;
      TH2D *mHHelium3PEta = nullptr;
      TH2D *mHHelium4PEta = nullptr;

      TH2D *mHAntiProtonPEta = nullptr;
      TH2D *mHAntiNeutronPEta = nullptr;
      TH2D *mHAntiDeuteronPEta = nullptr;
      TH2D *mHAntiTritiumPEta = nullptr;
      TH2D *mHAntiHelium3PEta = nullptr;
      TH2D *mHAntiHelium4PEta = nullptr;

      TH2D *mHProtonPtEta = nullptr;
      TH2D *mHNeutronPtEta = nullptr;
      TH2D *mHDeuteronPtEta = nullptr;
      TH2D *mHTritiumPtEta = nullptr;
      TH2D *mHHelium3PtEta = nullptr;
      TH2D *mHHelium4PtEta = nullptr;

      TH2D *mHAntiProtonPtEta = nullptr;
      TH2D *mHAntiNeutronPtEta = nullptr;
      TH2D *mHAntiDeuteronPtEta = nullptr;
      TH2D *mHAntiTritiumPtEta = nullptr;
      TH2D *mHAntiHelium3PtEta = nullptr;
      TH2D *mHAntiHelium4PtEta = nullptr;

      void initHistos() override;
      void initEventsHisto() override;
      void process() override;
      void process(double weight) override;
      void writeHistos() override;

      int selectP(const particleCand& p);

      int mPDGPr = 2212;
      int mPDGNe = 2112;
      int mPDGDe = 4324;
      int mPDGT = 2212 + 2*2112;
      int mPDGHe3 = 2*2212 + 2112;
      int mPDGHe4 = 2*2212 + 2*2112;

      int mPDGAntiPr = -2212;
      int mPDGAntiNe = -2112;
      int mPDGAntiDe = -4324;
      int mPDGAntiT = -2212 - 2*2112;
      int mPDGAntiHe3 = -2*2212 - 2112;
      int mPDGAntiHe4 = -2*2212 - 2*2112;

      float protonMass = 0.938272088;
      float neutronMass = 0.939565420;
      float deuteronMass = 1.875612942;
      float helium3Mass = 2.808391607;
      float tritiumMass = 2.808921132;
      float helium4Mass = 3.727379378;

      float mMinRapidity = -0.5;
      float mMaxRapidity = 0.5;

};

#endif