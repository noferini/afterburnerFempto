#include "CrunSpectrum.h"
#include "TFile.h"

void runSpectrum::initHistos()
{
    mHProtonSpectrum = new TH1D("mHProtonSpectrum" + mName, "Spectrum;p_{T} (GeV/c)",200, 0., 10.); //smaller width
    mHNeutronSpectrum = new TH1D("mHNeutronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHDeuteronSpectrum = new TH1D("mHDeuteronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHTritiumSpectrum = new TH1D("mHTritiumSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHHelium3Spectrum = new TH1D("mHHelium3Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHHelium4Spectrum = new TH1D("mHHelium4Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);

    mHProtonSpectrumA = new TH1D("mHProtonSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHNeutronSpectrumA = new TH1D("mHNeutronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHDeuteronSpectrumA = new TH1D("mHDeuteronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHTritiumSpectrumA = new TH1D("mHTritiumSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHHelium3SpectrumA = new TH1D("mHHelium3SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHHelium4SpectrumA = new TH1D("mHHelium4SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);

    mHProtonSpectrumM = new TH1D("mHProtonSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHNeutronSpectrumM = new TH1D("mHNeutronSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHDeuteronSpectrumM = new TH1D("mHDeuteronSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHTritiumSpectrumM = new TH1D("mHTritiumSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHHelium3SpectrumM = new TH1D("mHHelium3SpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHHelium4SpectrumM = new TH1D("mHHelium4SpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);

    mHProtonSpectrumARapidity = new TH2D("mHProtonSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHNeutronSpectrumARapidity = new TH2D("mHNeutronSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHDeuteronSpectrumARapidity = new TH2D("mHDeuteronSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHTritiumSpectrumARapidity = new TH2D("mHTritiumSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHHelium3SpectrumARapidity = new TH2D("mHHelium3SpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHHelium4SpectrumARapidity = new TH2D("mHHelium4SpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);

    mHProtonSpectrumMRapidity = new TH2D("mHProtonSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHNeutronSpectrumMRapidity = new TH2D("mHNeutronSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHDeuteronSpectrumMRapidity = new TH2D("mHDeuteronSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHTritiumSpectrumMRapidity = new TH2D("mHTritiumSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHHelium3SpectrumMRapidity = new TH2D("mHHelium3SpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHHelium4SpectrumMRapidity = new TH2D("mHHelium4SpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);

    mHProtonSpectrumAPEta = new TH2D("mHProtonSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHNeutronSpectrumAPEta = new TH2D("mHNeutronSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHDeuteronSpectrumAPEta = new TH2D("mHDeuteronSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHTritiumSpectrumAPEta = new TH2D("mHTritiumSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium3SpectrumAPEta = new TH2D("mHHelium3SpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium4SpectrumAPEta = new TH2D("mHHelium4SpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);

    mHProtonSpectrumMPEta = new TH2D("mHProtonSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHNeutronSpectrumMPEta = new TH2D("mHNeutronSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHDeuteronSpectrumMPEta = new TH2D("mHDeuteronSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHTritiumSpectrumMPEta = new TH2D("mHTritiumSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium3SpectrumMPEta = new TH2D("mHHelium3SpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium4SpectrumMPEta = new TH2D("mHHelium4SpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);

    mHProtonSpectrumAPtEta = new TH2D("mHProtonSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHNeutronSpectrumAPtEta = new TH2D("mHNeutronSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHDeuteronSpectrumAPtEta = new TH2D("mHDeuteronSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHTritiumSpectrumAPtEta = new TH2D("mHTritiumSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium3SpectrumAPtEta = new TH2D("mHHelium3SpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium4SpectrumAPtEta = new TH2D("mHHelium4SpectrumAEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);

    mHProtonSpectrumMPtEta = new TH2D("mHProtonSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHNeutronSpectrumMPtEta = new TH2D("mHNeutronSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHDeuteronSpectrumMPtEta = new TH2D("mHDeuteronSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHTritiumSpectrumMPtEta = new TH2D("mHTritiumSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium3SpectrumMPtEta = new TH2D("mHHelium3SpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHHelium4SpectrumMPtEta = new TH2D("mHHelium4SpectrumMEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);

    mHAntiProtonSpectrum = new TH1D("mHAntiProtonSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHAntiNeutronSpectrum = new TH1D("mHAntiNeutronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHAntiDeuteronSpectrum = new TH1D("mHAntiDeuteronSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHAntiTritiumSpectrum = new TH1D("mHAntiTritiumSpectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHAntiHelium3Spectrum = new TH1D("mHAntiHelium3Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);
    mHAntiHelium4Spectrum = new TH1D("mHAntiHelium4Spectrum" + mName, "Spectrum;p_{T} (GeV/c)", 200, 0., 10.);

    mHAntiProtonSpectrumA = new TH1D("mHAntiProtonSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHAntiNeutronSpectrumA = new TH1D("mHAntiNeutronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHAntiDeuteronSpectrumA = new TH1D("mHAntiDeuteronSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHAntiTritiumSpectrumA = new TH1D("mHAntiTritiumSpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHAntiHelium3SpectrumA = new TH1D("mHAntiHelium3SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);
    mHAntiHelium4SpectrumA = new TH1D("mHAntiHelium4SpectrumA" + mName, "Spectrum;p_{T}/A (GeV/c)", 200, 0., 10.);

    mHAntiProtonSpectrumM = new TH1D("mHAntiProtonSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHAntiNeutronSpectrumM = new TH1D("mHAntiNeutronSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHAntiDeuteronSpectrumM = new TH1D("mHAntiDeuteronSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHAntiTritiumSpectrumM = new TH1D("mHAntiTritiumSpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHAntiHelium3SpectrumM = new TH1D("mHAntiHelium3SpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);
    mHAntiHelium4SpectrumM = new TH1D("mHAntiHelium4SpectrumM" + mName, "Spectrum;p_{T}/M", 200, 0., 10.);

    mHAntiProtonSpectrumARapidity = new TH2D("mHAntiProtonSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiNeutronSpectrumARapidity = new TH2D("mHAntiNeutronSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiDeuteronSpectrumARapidity = new TH2D("mHAntiDeuteronSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiTritiumSpectrumARapidity = new TH2D("mHAntiTritiumSpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiHelium3SpectrumARapidity = new TH2D("mHAntiHelium3SpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiHelium4SpectrumARapidity = new TH2D("mHAntiHelium4SpectrumARapidity" + mName, "Spectrum;p_{T}/A (GeV/c);Rapidity", 100, 0., 10., 400, -10., 10.);

    mHAntiProtonSpectrumMRapidity = new TH2D("mHAntiProtonSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiNeutronSpectrumMRapidity = new TH2D("mHAntiNeutronSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiDeuteronSpectrumMRapidity = new TH2D("mHAntiDeuteronSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiTritiumSpectrumMRapidity = new TH2D("mHAntiTritiumSpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiHelium3SpectrumMRapidity = new TH2D("mHAntiHelium3SpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);
    mHAntiHelium4SpectrumMRapidity = new TH2D("mHAntiHelium4SpectrumMRapidity" + mName, "Spectrum;p_{T}/M;Rapidity", 100, 0., 10., 400, -10., 10.);

    mHAntiProtonSpectrumAPEta = new TH2D("mHAntiProtonSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiNeutronSpectrumAPEta = new TH2D("mHAntiNeutronSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiDeuteronSpectrumAPEta = new TH2D("mHAntiDeuteronSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiTritiumSpectrumAPEta = new TH2D("mHAntiTritiumSpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium3SpectrumAPEta = new TH2D("mHAntiHelium3SpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium4SpectrumAPEta = new TH2D("mHAntiHelium4SpectrumAPEta" + mName, "Spectrum;p/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);

    mHAntiProtonSpectrumMPEta = new TH2D("mHAntiProtonSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiNeutronSpectrumMPEta = new TH2D("mHAntiNeutronSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiDeuteronSpectrumMPEta = new TH2D("mHAntiDeuteronSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiTritiumSpectrumMPEta = new TH2D("mHAntiTritiumSpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium3SpectrumMPEta = new TH2D("mHAntiHelium3SpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium4SpectrumMPEta = new TH2D("mHAntiHelium4SpectrumMPEta" + mName, "Spectrum;p/M;Eta", 100, 0., 10., 1000, -50., 50);

    mHAntiProtonSpectrumAPtEta = new TH2D("mHAntiProtonSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiNeutronSpectrumAPtEta = new TH2D("mHAntiNeutronSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiDeuteronSpectrumAPtEta = new TH2D("mHAntiDeuteronSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiTritiumSpectrumAPtEta = new TH2D("mHAntiTritiumSpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium3SpectrumAPtEta = new TH2D("mHAntiHelium3SpectrumAPtEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium4SpectrumAPtEta = new TH2D("mHAntiHelium4SpectrumAEta" + mName, "Spectrum;p_{T}/A (GeV/c);Eta", 100, 0., 10., 1000, -50., 50);

    mHAntiProtonSpectrumMPtEta = new TH2D("mHAntiProtonSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiNeutronSpectrumMPtEta = new TH2D("mHAntiNeutronSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiDeuteronSpectrumMPtEta = new TH2D("mHAntiDeuteronSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiTritiumSpectrumMPtEta = new TH2D("mHAntiTritiumSpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium3SpectrumMPtEta = new TH2D("mHAntiHelium3SpectrumMPtEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);
    mHAntiHelium4SpectrumMPtEta = new TH2D("mHAntiHelium4SpectrumMEta" + mName, "Spectrum;p_{T}/M;Eta", 100, 0., 10., 1000, -50., 50);

    mHProtonPEta = new TH2D("mHProtonPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHNeutronPEta = new TH2D("mHNeutronPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHDeuteronPEta = new TH2D("mHDeuteronPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHTritiumPEta = new TH2D("mHTritiumPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHHelium3PEta = new TH2D("mHHelium3PEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHHelium4PEta = new TH2D("mHHelium4PEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    
    mHAntiProtonPEta = new TH2D("mHAntiProtonPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHAntiNeutronPEta = new TH2D("mHAntiNeutronPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHAntiDeuteronPEta = new TH2D("mHAntiDeuteronPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHAntiTritiumPEta = new TH2D("mHAntiTritiumPEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHAntiHelium3PEta = new TH2D("mHAntiHelium3PEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);
    mHAntiHelium4PEta = new TH2D("mHAntiHelium4PEta" + mName, "p vs Eta distribution;p (GeV/c);Eta", 1000 , 0., 100., 1000, -50., 50);

    mHProtonPtEta = new TH2D("mHProtonPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHNeutronPtEta = new TH2D("mHNeutronPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHDeuteronPtEta = new TH2D("mHDeuteronPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHTritiumPtEta = new TH2D("mHTritiumPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHHelium3PtEta = new TH2D("mHHelium3PtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHHelium4PtEta = new TH2D("mHHelium4PtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    
    mHAntiProtonPtEta = new TH2D("mHAntiProtonPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHAntiNeutronPtEta = new TH2D("mHAntiNeutronPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHAntiDeuteronPtEta = new TH2D("mHAntiDeuteronPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHAntiTritiumPtEta = new TH2D("mHAntiTritiumPtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHAntiHelium3PtEta = new TH2D("mHAntiHelium3PtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
    mHAntiHelium4PtEta = new TH2D("mHAntiHelium4PtEta" + mName, "p_{T} vs Eta distribution;p_{T} (GeV/c);Eta", 100 , 0., 10., 1000, -50., 50);
}
//_________________________________________________________________________
void runSpectrum::initEventsHisto()
{
    mEvents = new TH1D("mEvents" + mName, "Number of events", 16, 0, 16);
    mEvents->Fill("Number of unweighted Events", 0);
    mEvents->Fill("Number of Events", 0);
    mEvents->Fill("Number of accepted Events", 0);
    mEvents->Fill("Selected Particles", 0);
    mEvents->Fill("Selected Protons", 0);
    mEvents->Fill("Selected Neutrons", 0);
    mEvents->Fill("Selected Deuterons", 0);
    mEvents->Fill("Selected Tritons", 0);
    mEvents->Fill("Selected Helium3", 0);
    mEvents->Fill("Selected Helium4", 0);
    mEvents->Fill("Selected Antiprotons", 0);
    mEvents->Fill("Selected Antineutrons", 0);
    mEvents->Fill("Selected Antideuterons", 0);
    mEvents->Fill("Selected Antitritons", 0);
    mEvents->Fill("Selected Antihelium3", 0);
    mEvents->Fill("Selected Antihelium4", 0);
}
//_________________________________________________________________________
void runSpectrum::process()
{
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdg = selectP(p);
        if (ipdg == -1)
                continue;
        float y = p.q[ipdg].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity)
        {
            mEvents->Fill("Selected Particles", 1);
            if (p.pdgOptions[ipdg] == mPDGPr)
            {
                mEvents->Fill("Selected Protons", 1);
                mHProtonSpectrum->Fill(p.q[ipdg].Pt());
                mHProtonSpectrumA->Fill(p.q[ipdg].Pt());
                mHProtonSpectrumM->Fill(p.q[ipdg].Pt()/protonMass);
                mHProtonSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity());
                mHProtonSpectrumMRapidity->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Rapidity());
                mHProtonPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHProtonPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHProtonSpectrumAPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHProtonSpectrumMPEta->Fill(p.q[ipdg].P()/protonMass, p.q[ipdg].Eta());
                mHProtonSpectrumAPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHProtonSpectrumMPtEta->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Eta());    
            }else if(p.pdgOptions[ipdg] == mPDGNe){
                mEvents->Fill("Selected Neutrons", 1);
                mHNeutronSpectrum->Fill(p.q[ipdg].Pt());
                mHNeutronSpectrumA->Fill(p.q[ipdg].Pt());
                mHNeutronSpectrumM->Fill(p.q[ipdg].Pt()/neutronMass);
                mHNeutronSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity());
                mHNeutronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Rapidity());
                mHNeutronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHNeutronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHNeutronSpectrumAPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHNeutronSpectrumMPEta->Fill(p.q[ipdg].P()/neutronMass, p.q[ipdg].Eta());
                mHNeutronSpectrumAPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHNeutronSpectrumMPtEta->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGDe){
                mEvents->Fill("Selected Deuterons", 1);
                mHDeuteronSpectrum->Fill(p.q[ipdg].Pt());
                mHDeuteronSpectrumA->Fill(p.q[ipdg].Pt()/2);
                mHDeuteronSpectrumM->Fill(p.q[ipdg].Pt()/deuteronMass);
                mHDeuteronSpectrumARapidity->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Rapidity());
                mHDeuteronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Rapidity());
                mHDeuteronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHDeuteronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHDeuteronSpectrumAPEta->Fill(p.q[ipdg].P()/2, p.q[ipdg].Eta());
                mHDeuteronSpectrumMPEta->Fill(p.q[ipdg].P()/deuteronMass, p.q[ipdg].Eta());
                mHDeuteronSpectrumAPtEta->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Eta());
                mHDeuteronSpectrumMPtEta->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGT){
                mEvents->Fill("Selected Tritons", 1);
                mHTritiumSpectrum->Fill(p.q[ipdg].Pt());
                mHTritiumSpectrumA->Fill(p.q[ipdg].Pt()/3);
                mHTritiumSpectrumM->Fill(p.q[ipdg].Pt()/tritiumMass);
                mHTritiumSpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity());
                mHTritiumSpectrumMRapidity->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Rapidity());
                mHTritiumPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHTritiumPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHTritiumSpectrumAPEta->Fill(p.q[ipdg].P()/3, p.q[ipdg].Eta());
                mHTritiumSpectrumMPEta->Fill(p.q[ipdg].P()/tritiumMass, p.q[ipdg].Eta());
                mHTritiumSpectrumAPtEta->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Eta());
                mHTritiumSpectrumMPtEta->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGHe3){
                mEvents->Fill("Selected Helium3", 1);
                mHHelium3Spectrum->Fill(p.q[ipdg].Pt());
                mHHelium3SpectrumA->Fill(p.q[ipdg].Pt()/3);
                mHHelium3SpectrumM->Fill(p.q[ipdg].Pt()/helium3Mass);
                mHHelium3SpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity());
                mHHelium3SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Rapidity());
                mHHelium3PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHHelium3PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHHelium3SpectrumAPEta->Fill(p.q[ipdg].P()/3, p.q[ipdg].Eta());
                mHHelium3SpectrumMPEta->Fill(p.q[ipdg].P()/helium3Mass, p.q[ipdg].Eta());
                mHHelium3SpectrumAPtEta->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Eta());
                mHHelium3SpectrumMPtEta->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGHe4){
                mEvents->Fill("Selected Helium4", 1);
                mHHelium4Spectrum->Fill(p.q[ipdg].Pt());
                mHHelium4SpectrumA->Fill(p.q[ipdg].Pt()/4);
                mHHelium4SpectrumM->Fill(p.q[ipdg].Pt()/helium4Mass);
                mHHelium4SpectrumARapidity->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Rapidity());
                mHHelium4SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Rapidity());
                mHHelium4PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHHelium4PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHHelium4SpectrumAPEta->Fill(p.q[ipdg].P()/4, p.q[ipdg].Eta());
                mHHelium4SpectrumMPEta->Fill(p.q[ipdg].P()/helium4Mass, p.q[ipdg].Eta());
                mHHelium4SpectrumAPtEta->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Eta());
                mHHelium4SpectrumMPtEta->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Eta());
            }else if (p.pdgOptions[ipdg] == mPDGAntiPr){
                mEvents->Fill("Selected Antiprotons", 1);
                mHAntiProtonSpectrum->Fill(p.q[ipdg].Pt());
                mHAntiProtonSpectrumA->Fill(p.q[ipdg].Pt());
                mHAntiProtonSpectrumM->Fill(p.q[ipdg].Pt()/protonMass);
                mHAntiProtonSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity());
                mHAntiProtonSpectrumMRapidity->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Rapidity());
                mHAntiProtonPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiProtonPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiProtonSpectrumAPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiProtonSpectrumMPEta->Fill(p.q[ipdg].P()/protonMass, p.q[ipdg].Eta());
                mHAntiProtonSpectrumAPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiProtonSpectrumMPtEta->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Eta());  
            }else if(p.pdgOptions[ipdg] == mPDGAntiNe){
                mEvents->Fill("Selected Antineutrons", 1);
                mHAntiNeutronSpectrum->Fill(p.q[ipdg].Pt());
                mHAntiNeutronSpectrumA->Fill(p.q[ipdg].Pt());
                mHAntiNeutronSpectrumM->Fill(p.q[ipdg].Pt()/neutronMass);
                mHAntiNeutronSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity());
                mHAntiNeutronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Rapidity());
                mHAntiNeutronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiNeutronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiNeutronSpectrumAPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiNeutronSpectrumMPEta->Fill(p.q[ipdg].P()/neutronMass, p.q[ipdg].Eta());
                mHAntiNeutronSpectrumAPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiNeutronSpectrumMPtEta->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGAntiDe){
                mEvents->Fill("Selected Antideuterons", 1);
                mHAntiDeuteronSpectrum->Fill(p.q[ipdg].Pt());
                mHAntiDeuteronSpectrumA->Fill(p.q[ipdg].Pt()/2);
                mHAntiDeuteronSpectrumM->Fill(p.q[ipdg].Pt()/deuteronMass);
                mHAntiDeuteronSpectrumARapidity->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Rapidity());
                mHAntiDeuteronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Rapidity());
                mHAntiDeuteronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiDeuteronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiDeuteronSpectrumAPEta->Fill(p.q[ipdg].P()/2, p.q[ipdg].Eta());
                mHAntiDeuteronSpectrumMPEta->Fill(p.q[ipdg].P()/deuteronMass, p.q[ipdg].Eta());
                mHAntiDeuteronSpectrumAPtEta->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Eta());
                mHAntiDeuteronSpectrumMPtEta->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGAntiT){
                mEvents->Fill("Selected Antitritons", 1);
                mHAntiTritiumSpectrum->Fill(p.q[ipdg].Pt());
                mHAntiTritiumSpectrumA->Fill(p.q[ipdg].Pt()/3);
                mHAntiTritiumSpectrumM->Fill(p.q[ipdg].Pt()/tritiumMass);
                mHAntiTritiumSpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity());
                mHAntiTritiumSpectrumMRapidity->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Rapidity());
                mHAntiTritiumPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiTritiumPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiTritiumSpectrumAPEta->Fill(p.q[ipdg].P()/3, p.q[ipdg].Eta());
                mHAntiTritiumSpectrumMPEta->Fill(p.q[ipdg].P()/tritiumMass, p.q[ipdg].Eta());
                mHAntiTritiumSpectrumAPtEta->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Eta());
                mHAntiTritiumSpectrumMPtEta->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGAntiHe3){
                mEvents->Fill("Selected Antihelium3", 1);
                mHAntiHelium3Spectrum->Fill(p.q[ipdg].Pt());
                mHAntiHelium3SpectrumA->Fill(p.q[ipdg].Pt()/3);
                mHAntiHelium3SpectrumM->Fill(p.q[ipdg].Pt()/helium3Mass);
                mHAntiHelium3SpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity());
                mHAntiHelium3SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Rapidity());
                mHAntiHelium3PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiHelium3PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiHelium3SpectrumAPEta->Fill(p.q[ipdg].P()/3, p.q[ipdg].Eta());
                mHAntiHelium3SpectrumMPEta->Fill(p.q[ipdg].P()/helium3Mass, p.q[ipdg].Eta());
                mHAntiHelium3SpectrumAPtEta->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Eta());
                mHAntiHelium3SpectrumMPtEta->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Eta());
            }else if(p.pdgOptions[ipdg] == mPDGAntiHe4){
                mEvents->Fill("Selected Antihelium4", 1);
                mHAntiHelium4Spectrum->Fill(p.q[ipdg].Pt());
                mHAntiHelium4SpectrumA->Fill(p.q[ipdg].Pt()/4);
                mHAntiHelium4SpectrumM->Fill(p.q[ipdg].Pt()/helium4Mass);
                mHAntiHelium4SpectrumARapidity->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Rapidity());
                mHAntiHelium4SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Rapidity());
                mHAntiHelium4PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta());
                mHAntiHelium4PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta());
                mHAntiHelium4SpectrumAPEta->Fill(p.q[ipdg].P()/4, p.q[ipdg].Eta());
                mHAntiHelium4SpectrumMPEta->Fill(p.q[ipdg].P()/helium4Mass, p.q[ipdg].Eta());
                mHAntiHelium4SpectrumAPtEta->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Eta());
                mHAntiHelium4SpectrumMPtEta->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Eta());
            }
        }
    }
}
//_________________________________________________________________________
void runSpectrum::process(double weight)
{
    for (int i = 0; i < mVect.size(); ++i)
    {
        const particleCand &p = mVect[i];
        int ipdg = selectP(p);
        if (ipdg == -1)
                continue;
        float y = p.q[ipdg].Rapidity();
        if (y >= mMinRapidity && y <= mMaxRapidity)
        {
            mEvents->Fill("Selected Particles", weight);
            if (p.pdgOptions[ipdg] == mPDGPr)
            {
                mEvents->Fill("Selected Protons", weight);
                mHProtonSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHProtonSpectrumA->Fill(p.q[ipdg].Pt(), weight);
                mHProtonSpectrumM->Fill(p.q[ipdg].Pt()/protonMass, weight);
                mHProtonSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity(), weight);
                mHProtonSpectrumMRapidity->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Rapidity(), weight);
                mHProtonPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHProtonPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGNe){
                mEvents->Fill("Selected Neutrons", weight);
                mHNeutronSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHNeutronSpectrumA->Fill(p.q[ipdg].Pt(), weight);
                mHNeutronSpectrumM->Fill(p.q[ipdg].Pt()/neutronMass, weight);
                mHNeutronSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity(), weight);
                mHNeutronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Rapidity(), weight);
                mHNeutronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHNeutronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGDe){
                mEvents->Fill("Selected Deuterons", weight);
                mHDeuteronSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHDeuteronSpectrumA->Fill(p.q[ipdg].Pt()/2, weight);
                mHDeuteronSpectrumM->Fill(p.q[ipdg].Pt()/deuteronMass, weight);
                mHDeuteronSpectrumARapidity->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Rapidity(), weight);
                mHDeuteronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Rapidity(), weight);
                mHDeuteronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHDeuteronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGT){
                mEvents->Fill("Selected Tritons", weight);
                mHTritiumSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHTritiumSpectrumA->Fill(p.q[ipdg].Pt()/3, weight);
                mHTritiumSpectrumM->Fill(p.q[ipdg].Pt()/tritiumMass, weight);
                mHTritiumSpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity(), weight);
                mHTritiumSpectrumMRapidity->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Rapidity(), weight);
                mHTritiumPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHTritiumPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGHe3){
                mEvents->Fill("Selected Helium3", weight);
                mHHelium3Spectrum->Fill(p.q[ipdg].Pt(), weight);
                mHHelium3SpectrumA->Fill(p.q[ipdg].Pt()/3, weight);
                mHHelium3SpectrumM->Fill(p.q[ipdg].Pt()/helium3Mass, weight);
                mHHelium3SpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity(), weight);
                mHHelium3SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Rapidity(), weight);
                mHHelium3PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHHelium3PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGHe4){
                mEvents->Fill("Selected Helium4", weight);
                mHHelium4Spectrum->Fill(p.q[ipdg].Pt(), weight);
                mHHelium4SpectrumA->Fill(p.q[ipdg].Pt()/4, weight);
                mHHelium4SpectrumM->Fill(p.q[ipdg].Pt()/helium4Mass, weight);
                mHHelium4SpectrumARapidity->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Rapidity(), weight);
                mHHelium4SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Rapidity(), weight);
                mHHelium4PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHHelium4PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if (p.pdgOptions[ipdg] == mPDGAntiPr){
                mEvents->Fill("Selected Antiprotons", weight);
                mHAntiProtonSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiProtonSpectrumA->Fill(p.q[ipdg].Pt(), weight);
                mHAntiProtonSpectrumM->Fill(p.q[ipdg].Pt()/protonMass, weight);
                mHAntiProtonSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity(), weight);
                mHAntiProtonSpectrumMRapidity->Fill(p.q[ipdg].Pt()/protonMass, p.q[ipdg].Rapidity(), weight);
                mHAntiProtonPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiProtonPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGAntiNe){
                mEvents->Fill("Selected Antineutrons", weight);
                mHAntiNeutronSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiNeutronSpectrumA->Fill(p.q[ipdg].Pt(), weight);
                mHAntiNeutronSpectrumM->Fill(p.q[ipdg].Pt()/neutronMass, weight);
                mHAntiNeutronSpectrumARapidity->Fill(p.q[ipdg].Pt(), p.q[ipdg].Rapidity(), weight);
                mHAntiNeutronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/neutronMass, p.q[ipdg].Rapidity(), weight);
                mHAntiNeutronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiNeutronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGAntiDe){
                mEvents->Fill("Selected Antideuterons", weight);
                mHAntiDeuteronSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiDeuteronSpectrumA->Fill(p.q[ipdg].Pt()/2, weight);
                mHAntiDeuteronSpectrumM->Fill(p.q[ipdg].Pt()/deuteronMass, weight);
                mHAntiDeuteronSpectrumARapidity->Fill(p.q[ipdg].Pt()/2, p.q[ipdg].Rapidity(), weight);
                mHAntiDeuteronSpectrumMRapidity->Fill(p.q[ipdg].Pt()/deuteronMass, p.q[ipdg].Rapidity(), weight);
                mHAntiDeuteronPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiDeuteronPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGAntiT){
                mEvents->Fill("Selected Antitritons", weight);
                mHAntiTritiumSpectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiTritiumSpectrumA->Fill(p.q[ipdg].Pt()/3, weight);
                mHAntiTritiumSpectrumM->Fill(p.q[ipdg].Pt()/tritiumMass, weight);
                mHAntiTritiumSpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity(), weight);
                mHAntiTritiumSpectrumMRapidity->Fill(p.q[ipdg].Pt()/tritiumMass, p.q[ipdg].Rapidity(), weight);
                mHAntiTritiumPEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiTritiumPtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGAntiHe3){
                mEvents->Fill("Selected Antihelium3", weight);
                mHAntiHelium3Spectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiHelium3SpectrumA->Fill(p.q[ipdg].Pt()/3, weight);
                mHAntiHelium3SpectrumM->Fill(p.q[ipdg].Pt()/helium3Mass, weight);
                mHAntiHelium3SpectrumARapidity->Fill(p.q[ipdg].Pt()/3, p.q[ipdg].Rapidity(), weight);
                mHAntiHelium3SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium3Mass, p.q[ipdg].Rapidity(), weight);
                mHAntiHelium3PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiHelium3PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }else if(p.pdgOptions[ipdg] == mPDGAntiHe4){
                mEvents->Fill("Selected Antihelium4", weight);
                mHAntiHelium4Spectrum->Fill(p.q[ipdg].Pt(), weight);
                mHAntiHelium4SpectrumA->Fill(p.q[ipdg].Pt()/4, weight);
                mHAntiHelium4SpectrumM->Fill(p.q[ipdg].Pt()/helium4Mass, weight);
                mHAntiHelium4SpectrumARapidity->Fill(p.q[ipdg].Pt()/4, p.q[ipdg].Rapidity(), weight);
                mHAntiHelium4SpectrumMRapidity->Fill(p.q[ipdg].Pt()/helium4Mass, p.q[ipdg].Rapidity(), weight);
                mHAntiHelium4PEta->Fill(p.q[ipdg].P(), p.q[ipdg].Eta(), weight);
                mHAntiHelium4PtEta->Fill(p.q[ipdg].Pt(), p.q[ipdg].Eta(), weight);
            }
        }
    }
}
//_________________________________________________________________________
void runSpectrum::writeHistos()
{
    mHProtonSpectrum->Write();
    mHNeutronSpectrum->Write();
    mHDeuteronSpectrum->Write();
    mHTritiumSpectrum->Write();
    mHHelium3Spectrum->Write();
    mHHelium4Spectrum->Write();
    mHProtonSpectrumA->Write();
    mHNeutronSpectrumA->Write();
    mHDeuteronSpectrumA->Write();
    mHTritiumSpectrumA->Write();
    mHHelium3SpectrumA->Write();
    mHHelium4SpectrumA->Write();
    mHProtonSpectrumM->Write();
    mHNeutronSpectrumM->Write();
    mHDeuteronSpectrumM->Write();
    mHTritiumSpectrumM->Write();
    mHHelium3SpectrumM->Write();
    mHHelium4SpectrumM->Write();
    mHProtonSpectrumARapidity->Write();
    mHNeutronSpectrumARapidity->Write();
    mHDeuteronSpectrumARapidity->Write();
    mHTritiumSpectrumARapidity->Write();
    mHHelium3SpectrumARapidity->Write();
    mHHelium4SpectrumARapidity->Write();
    mHProtonSpectrumMRapidity->Write();
    mHNeutronSpectrumMRapidity->Write();
    mHDeuteronSpectrumMRapidity->Write();
    mHTritiumSpectrumMRapidity->Write();
    mHHelium3SpectrumMRapidity->Write();
    mHHelium4SpectrumMRapidity->Write();
    mHProtonSpectrumAPEta->Write();
    mHNeutronSpectrumAPEta->Write();
    mHDeuteronSpectrumAPEta->Write();
    mHTritiumSpectrumAPEta->Write();
    mHHelium3SpectrumAPEta->Write();
    mHHelium4SpectrumAPEta->Write();
    mHProtonSpectrumMPEta->Write();
    mHNeutronSpectrumMPEta->Write();
    mHDeuteronSpectrumMPEta->Write();
    mHTritiumSpectrumMPEta->Write();
    mHHelium3SpectrumMPEta->Write();
    mHHelium4SpectrumMPEta->Write();
    mHProtonSpectrumAPtEta->Write();
    mHNeutronSpectrumAPtEta->Write();
    mHDeuteronSpectrumAPtEta->Write();
    mHTritiumSpectrumAPtEta->Write();
    mHHelium3SpectrumAPtEta->Write();
    mHHelium4SpectrumAPtEta->Write();
    mHProtonSpectrumMPtEta->Write();
    mHNeutronSpectrumMPtEta->Write();
    mHDeuteronSpectrumMPtEta->Write();
    mHTritiumSpectrumMPtEta->Write();
    mHHelium3SpectrumMPtEta->Write();
    mHHelium4SpectrumMPtEta->Write();
    mHAntiProtonSpectrum->Write();
    mHAntiNeutronSpectrum->Write();
    mHAntiDeuteronSpectrum->Write();
    mHAntiTritiumSpectrum->Write();
    mHAntiHelium3Spectrum->Write();
    mHAntiHelium4Spectrum->Write();
    mHAntiProtonSpectrumA->Write();
    mHAntiNeutronSpectrumA->Write();
    mHAntiDeuteronSpectrumA->Write();
    mHAntiTritiumSpectrumA->Write();
    mHAntiHelium3SpectrumA->Write();
    mHAntiHelium4SpectrumA->Write();
    mHAntiProtonSpectrumM->Write();
    mHAntiNeutronSpectrumM->Write();
    mHAntiDeuteronSpectrumM->Write();
    mHAntiTritiumSpectrumM->Write();
    mHAntiHelium3SpectrumM->Write();
    mHAntiHelium4SpectrumM->Write();
    mHAntiProtonSpectrumARapidity->Write();
    mHAntiNeutronSpectrumARapidity->Write();
    mHAntiDeuteronSpectrumARapidity->Write();
    mHAntiTritiumSpectrumARapidity->Write();
    mHAntiHelium3SpectrumARapidity->Write();
    mHAntiHelium4SpectrumARapidity->Write();
    mHAntiProtonSpectrumMRapidity->Write();
    mHAntiNeutronSpectrumMRapidity->Write();
    mHAntiDeuteronSpectrumMRapidity->Write();
    mHAntiTritiumSpectrumMRapidity->Write();
    mHAntiHelium3SpectrumMRapidity->Write();
    mHAntiHelium4SpectrumMRapidity->Write();
    mHAntiProtonSpectrumAPEta->Write();
    mHAntiNeutronSpectrumAPEta->Write();
    mHAntiDeuteronSpectrumAPEta->Write();
    mHAntiTritiumSpectrumAPEta->Write();
    mHAntiHelium3SpectrumAPEta->Write();
    mHAntiHelium4SpectrumAPEta->Write();
    mHAntiProtonSpectrumMPEta->Write();
    mHAntiNeutronSpectrumMPEta->Write();
    mHAntiDeuteronSpectrumMPEta->Write();
    mHAntiTritiumSpectrumMPEta->Write();
    mHAntiHelium3SpectrumMPEta->Write();
    mHAntiHelium4SpectrumMPEta->Write();
    mHAntiProtonSpectrumAPtEta->Write();
    mHAntiNeutronSpectrumAPtEta->Write();
    mHAntiDeuteronSpectrumAPtEta->Write();
    mHAntiTritiumSpectrumAPtEta->Write();
    mHAntiHelium3SpectrumAPtEta->Write();
    mHAntiHelium4SpectrumAPtEta->Write();
    mHAntiProtonSpectrumMPtEta->Write();
    mHAntiNeutronSpectrumMPtEta->Write();
    mHAntiDeuteronSpectrumMPtEta->Write();
    mHAntiTritiumSpectrumMPtEta->Write();
    mHAntiHelium3SpectrumMPtEta->Write();
    mHAntiHelium4SpectrumMPtEta->Write();
    mHProtonPEta->Write();
    mHNeutronPEta->Write();
    mHDeuteronPEta->Write();
    mHTritiumPEta->Write();
    mHHelium3PEta->Write();
    mHHelium4PEta->Write();
    mHAntiProtonPEta->Write();
    mHAntiNeutronPEta->Write();
    mHAntiDeuteronPEta->Write();
    mHAntiTritiumPEta->Write();
    mHAntiHelium3PEta->Write();
    mHAntiHelium4PEta->Write();
    mHProtonPtEta->Write();
    mHNeutronPtEta->Write();
    mHDeuteronPtEta->Write();
    mHTritiumPtEta->Write();
    mHHelium3PtEta->Write();
    mHHelium4PtEta->Write();
    mHAntiProtonPtEta->Write();
    mHAntiNeutronPtEta->Write();
    mHAntiDeuteronPtEta->Write();
    mHAntiTritiumPtEta->Write();
    mHAntiHelium3PtEta->Write();
    mHAntiHelium4PtEta->Write();
}
//_________________________________________________________________________
int runSpectrum::selectP(const particleCand &p)
{
    for (int i = 0; i < p.pdgOptions.size(); ++i)
    {
        bool particleCondition = p.pdgOptions[i] == mPDGPr || p.pdgOptions[i] == mPDGNe || p.pdgOptions[i] == mPDGDe || p.pdgOptions[i] == mPDGT || p.pdgOptions[i] == mPDGHe3 || p.pdgOptions[i] == mPDGHe4;
        bool antiparticleCondition = p.pdgOptions[i] == mPDGAntiPr || p.pdgOptions[i] == mPDGAntiNe || p.pdgOptions[i] == mPDGAntiDe || p.pdgOptions[i] == mPDGAntiT || p.pdgOptions[i] == mPDGAntiHe3 || p.pdgOptions[i] == mPDGAntiHe4;
        if (particleCondition || antiparticleCondition)
        {
            return i;
        }
    }
    return -1;
}
