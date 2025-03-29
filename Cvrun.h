#ifndef CVRUN
#define CVRUN

#define MAXMIXEDEVENTS 1000

#include "Cutils.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"

class vrun
{
  public:
    vrun(int nmix=10) : mNmix(TMath::Min(nmix,MAXMIXEDEVENTS)) {}

    void setNmixedEvents(int nmix) { mNmix = TMath::Min(nmix,MAXMIXEDEVENTS); }
    void init();
//    void setEvent(your custom input);
    void setEvent(const std::vector<particleMC>& vect);     // to be defined for your special input but then it should be mapped on particleCand
    void setIDName(const TString& name);
    virtual void process();
    virtual void finalize();
    virtual void write();

    void selectPDG(int pdg1, int pdg2) { mPDG1 = pdg1, mPDG2 = pdg2, mIsSamePart = (pdg1 == pdg2); }

    virtual TH2D* getKstarSE();
    virtual TH2D* getKstarME();
    virtual TH2D* getDPhiDEtaSE();
    virtual TH2D* getDPhiDEtaME();

  protected:
    std::vector<particleCand> mVect;
    virtual void initHistos();
    virtual void finalizeHistos();
    virtual int selectP1(const particleCand& p);
    virtual int selectP2(const particleCand& p);
  private:
    TString mName = "";
    int mNmix;       // number of mixed events, max is MAXMIXEDEVENTS
    std::vector<particleCand> mEvPrev[MAXMIXEDEVENTS];
    TH2D *mHkstarSE;
    TH2D *mHkstarME;
    TH2D *mDPhiDEtaSE;
    TH2D *mDPhiDEtaME;

    int mPDG1 = 2212;
    int mPDG2 = 2212;
    bool mIsSamePart = true;
};

#endif
