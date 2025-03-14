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
    virtual void process();
    virtual void finalize();
    virtual void write();

    void selectPDG(int pdg1, int pdg2) { mPDG1 = pdg1, mPDG2 = pdg2, mIsSamePart = (pdg1 == pdg2); }

  protected:
    std::vector<particleCand> mVect;
    virtual void initHistos();
    virtual void finalizeHistos();
    virtual int selectP1(const particleCand& p);
    virtual int selectP2(const particleCand& p);
  private:
    int mNmix;       // number of mixed events, max is MAXMIXEDEVENTS
    std::vector<particleCand> mEvPrev[MAXMIXEDEVENTS];
    TH2D *mHkstarSE;
    TH2D *mHkstarME;

    int mPDG1 = 2212;
    int mPDG2 = 2212;
    bool mIsSamePart = true;
};

#endif
