#ifndef CWAVEUTILS
#define CWAVEUTILS

#include "TF1.h"
#include "Cutils.h"

class waveUtils
{
  public:
    static void init();
    static float getCoalProb(const particleMC& p1, const particleMC& p2);
    static float calcProb();
    static void setParams(float strong=17.4, float strongR=3.2, float coloumb=1.44, float sourceRadius=0, float spinFact=3./8) { mStrong = strong, mStrongR = strongR, mCoulomb = coloumb, mSourceRadius = sourceRadius, mSpinCoalFactor = spinFact; }

    static void setSourceRadius(float radius);
    static void setKstar(float kstar, float kt=1.0, utils::type system=utils::pn);
    static void setCharges(float cS, float cC);

    static void setNuclearRadius(float radius);

    static void setIsRadiusMtDependent(bool val=true) { mIsRadiusMtDependent=val; }
    static bool isRadiusMtDependent() { return mIsRadiusMtDependent; }

    static TF1 *getCoalRe() { return mCoalescenceRe; }
    static TF1 *getCoalIm() { return mCoalescenceIm; }

    static TF1 *getUDeuteron() { return mUDeuteron; }
    static TF1 *getDeuteron() { return mDeuteron; }
    static TF1 *getDeuteron2int() { return mDeuteron2int; }
    static TF1 *getDeuteronKin() { return mDeuteronKin; }
    static TF1 *getDeuteronV() { return mDeuteronV; }

    static TF1 *getSource () { return mSource; }
    static TF1 *getUSource () { return mUSource; }
    static TF1 *getSource2int () { return mSource; }
    static TF1 *getSourceKin () { return mSourceKin; }
    static TF1 *getSourceV () { return mSourceV; }
    static float kineticSource() { double err; return mSourceKin->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    static float potentialSource() { double err; return mSourceV->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    static float getKstarFinal(float coalProb = 0., float massRed=MRED_NN, float boundE=EBOUND_D);
    static float kineticDeuteron() { double err; return mDeuteronKin->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }
    static float potentialDeuteron() { double err; return mDeuteronV->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err); }

    static bool mIsInitialized;

    static float calculateRadius(float EbindingEff=-2.22, float mred=938/2., float chargeStrong=1, float chargeCoulomb=0, const char* title="");
    static float calculateV0(float EbindingEff=-2.22, float mred=938/2., float chargeStrong=1, float chargeCoulomb=0, const char* title="");

    static constexpr double HCUT = 197.3;       // [fm*Mev]
    static constexpr double EBOUND_D = -2.22;   // [MeV] deuteron binding energy
    static constexpr double EBOUND_T = -8.48;   // [MeV] deuteron binding energy
    static constexpr double EBOUND_3HE = -7.72; // [MeV] deuteron binding energy
    static constexpr double EBOUND_HE = -28.3;  // [MeV] deuteron binding energy

    static constexpr float MRED_NN = 938/2.;    // [MeV] reducted mass for Nucleon-Nucleon
    static constexpr float MRED_DN = 2*938/3.;  // [MeV] reducted mass for Nucleon-Nucleon
    static constexpr float MRED_DD = 938;       // [MeV] reducted mass for Nucleon-Nucleon
    static constexpr float MRED_TN = 3*938/4.;  // [MeV] reducted mass for Nucleon-Nucleon
    static constexpr float MRED_HEN = 3*938/4.; // [MeV] reducted mass for Nucleon-Nucleon

    static float mStrong;                       // attractive potential
    static float mStrongR;                      // radius of box potential
    static float mCoulomb;
    static float mSourceRadius;
    static float mSpinCoalFactor;
    static float mStrongDn;
    static float mStrongDD;
    static float mStrongTn;
    static float mStrongHen;
    static float mStrongRDn;                      // radius of box potential
    static float mStrongRDD;                      // radius of box potential
    static float mStrongRTn;                      // radius of box potential
    static float mStrongRHen;                     // radius of box potential

    // deuteron functions
    static double uDeuteron(double *x,double *pm);
    static double waveDeuteron(double *x,double *pm);
    static double intDeuteron(double *x,double *pm);
    static double deuteronKin(double *x,double *pm);
    static double deuteronV(double *x,double *pm);
    static TF1 *mDeuteron;
    static TF1 *mUDeuteron;
    static TF1 *mDeuteron2int;
    static TF1 *mDeuteronKin;
    static TF1 *mDeuteronV;

    // source functions
    static float mMaxIntRange;
    static double source(double *x,double *pm);
    static double usource(double *x,double *pm);
    static double source2int(double *x,double *pm);
    static double sourceKin(double *x,double *pm);
    static double sourceV(double *x,double *pm);
    static double sourceCos(double *x,double *pm);
    static double sourceSin(double *x,double *pm);
    static TF1 *mSource;
    static TF1 *mUSource;
    static TF1 *mSource2int;
    static TF1 *mSourceKin;
    static TF1 *mSourceV;

    // coalescence
    static double coalescenceRe(double *x,double *pm);
    static double coalescenceIm(double *x,double *pm);
    static TF1 *mCoalescenceRe;
    static TF1 *mCoalescenceIm;

  private:
    static double mK1;        // sqrt(2*MRED*(EBOUND+V0))/HCUT;
    static double mK2;        // sqrt(-2*MRED*EBOUND)/HCUT;
    static double mNormRight; // sin(k1*mStrongR) * TMath::Exp(k2*mStrongR);
    static bool mIsRadiusMtDependent;
};

#endif
