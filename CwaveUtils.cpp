#include "CwaveUtils.h"

bool waveUtils::mIsInitialized = false;
float waveUtils::mStrong = 17.4;        // attractive potential
float waveUtils::mStrongR = 3.2;        // radius of box potential
float waveUtils::mCoulomb = 1.44E-3;
float waveUtils::mSourceRadius = -1.5;
float waveUtils::mSpinCoalFactor = 3./8;
float waveUtils::mMaxIntRange = 20;
double waveUtils::mK1 = 0;              // sqrt(2*MRED*(EBOUND+V0))/HCUT;
double waveUtils::mK2 = 0;              // sqrt(-2*MRED*EBOUND)/HCUT;
double waveUtils::mNormRight = 0;       // sin(k1*mStrongR) * TMath::Exp(k2*mStrongR);

TF1 *waveUtils::mDeuteron = nullptr;
TF1 *waveUtils::mUDeuteron = nullptr;
TF1 *waveUtils::mDeuteron2int = nullptr;
TF1 *waveUtils::mDeuteronKin = nullptr;
TF1 *waveUtils::mDeuteronV = nullptr;
TF1 *waveUtils::mSource = nullptr;
TF1 *waveUtils::mUSource = nullptr;
TF1 *waveUtils::mSource2int = nullptr;
TF1 *waveUtils::mSourceKin = nullptr;
TF1 *waveUtils::mSourceV = nullptr;
TF1 *waveUtils::mCoalescenceRe = nullptr;
TF1 *waveUtils::mCoalescenceIm = nullptr;

float waveUtils::getKstarFinal(float coalProb, float massRed, float boundE) {
  float Efin = TMath::Max(float(0.), (kineticSource() + potentialSource() - boundE*coalProb)/(1 - coalProb));

//  printf("K=%f - V=%f E=%f -> k*=%f\n",kineticSource(),potentialSource(),Efin,sqrt(Efin*2*massRed));

  return sqrt(Efin*2*massRed);
}
//_________________________________________________________________________
void waveUtils::setCharges(float cS, float cC){
  mSourceV->SetParameter(2,cS);
  mSourceV->SetParameter(3,cC);

  mDeuteronV->SetParameter(5,cS);
  mDeuteronV->SetParameter(6,cC);
}
//_________________________________________________________________________
void waveUtils::setKstar(float kstar, float kt) {
  if(! mIsInitialized){
    init();
  }
  if(mSourceRadius < 0){
    float radiusSource = std::abs(mSourceRadius/kt);
    static const float factor = sqrt(3*0.5) * HCUT;
    float radiusWave = factor/kstar;
    float radius = sqrt(radiusWave*radiusWave + radiusSource*radiusSource);
    float kstarWave = factor/radius;
    float kstarEff = sqrt(kstar*kstar - kstarWave*kstarWave);
    setSourceRadius(-radius);
    mSourceKin->SetParameter(2, kstarEff);
    mCoalescenceRe->SetParameter(2, kstarEff);
    mCoalescenceIm->SetParameter(2, kstarEff);
  } else {
    mSourceKin->SetParameter(2, kstar);
    mCoalescenceRe->SetParameter(2, kstar);
    mCoalescenceIm->SetParameter(2, kstar);
    setSourceRadius(mSourceRadius/kt); // to trigger new normalization
  }

}
//_________________________________________________________________________
void waveUtils::setSourceRadius(float radius) {
  if(! mIsInitialized){
    init();
  }
  if(radius > 0){
    mSourceRadius = radius;
  } else {
//    mSourceRadius = radius;
    radius *= -1;
  }

  mSource2int->SetParameter(1,radius);
  mSource->SetParameter(1,radius);
  mUSource->SetParameter(1,radius);
  mSourceKin->SetParameter(1,radius);
  mSourceV->SetParameter(1,radius);
  mCoalescenceRe->SetParameter(1,radius);
  mCoalescenceIm->SetParameter(1,radius);

  mSource2int->SetParameter(0,1);

  mMaxIntRange = TMath::Max(20., double(radius*3));
  if(mMaxIntRange > 100){
    mMaxIntRange = 100;
  }

  double err;
  double norm = 1./sqrt(mSource2int->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err));

  mSource->SetParameter(0,norm);
  mUSource->SetParameter(0,norm);
  mSource2int->SetParameter(0,norm);
  mSourceKin->SetParameter(0,norm);
  mSourceV->SetParameter(0,norm);
  mCoalescenceRe->SetParameter(0,norm);
  mCoalescenceIm->SetParameter(0,norm);
}
//_________________________________________________________________________
double waveUtils::coalescenceRe(double *x,double *pm){
  double r = x[0];
  double pDeu[5] = {pm[3], pm[4], pm[5], pm[6], pm[7]};

  return waveDeuteron(x,pDeu) * sourceCos(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double waveUtils::coalescenceIm(double *x,double *pm){
  double r = x[0];
  double pDeu[5] = {pm[3], pm[4], pm[5], pm[6], pm[7]};

  return waveDeuteron(x,pDeu) * sourceSin(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double waveUtils::uDeuteron(double *x,double *pm){
    double r = x[0];

    double res = pm[0];
    double radius = pm[1];
    double k1 = pm[2];
    double k2 = pm[3];
    double normRight = pm[4];

    if(r < radius){
      res *= sin(k1*r);
    } else {
      res *= TMath::Exp(-k2*r) * normRight;
    }

    return res;
}
//_________________________________________________________________________
double waveUtils::waveDeuteron(double *x,double *pm){
    double r = x[0];

    return uDeuteron(x,pm) / r;
}
//_________________________________________________________________________
double waveUtils::intDeuteron(double *x,double *pm){
    double r = x[0];

    return waveDeuteron(x,pm) * waveDeuteron(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double waveUtils::deuteronKin(double *x,double *pm){
    double r = x[0];

    double res = 0;
    double radius = pm[1];
    double k1 = pm[2];
    double k2 = pm[3];

    if(r < radius){
      res = k1*k1*HCUT*HCUT*0.5/MRED_NN;
    } else {
      res = -k2*k2*HCUT*HCUT*0.5/MRED_NN;
    }

  return res * intDeuteron(x,pm);
}
//_________________________________________________________________________
double waveUtils::deuteronV(double *x,double *pm){
  double r = x[0];
  float chargeStrong = pm[5];
  float chargeColoumb = pm[6];
  double val = 0;
  if(r < mStrongR){
    val = - mStrong * chargeStrong;
  }
  val += chargeColoumb * mCoulomb / r;

  return val * intDeuteron(x,pm);
}
//_________________________________________________________________________
double waveUtils::source(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    return norm * TMath::Exp(-r*r*0.5/(radius*radius));
}
//_________________________________________________________________________
double waveUtils::sourceCos(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    double kstar = pm[2];

//    return norm * TMath::Exp(-r*r*0.5/(radius*radius)) * cos(kstar*r/HCUT);
    return norm * TMath::Exp(-r*r*0.5/(radius*radius)) * sin(kstar*r/HCUT) * HCUT/(kstar*r);
}
//_________________________________________________________________________
double waveUtils::sourceSin(double *x,double *pm){
    double r = x[0];
    double norm = pm[0];
    double radius = pm[1];
    double kstar = pm[2];

//    return norm * TMath::Exp(-r*r*0.5/(radius*radius)) * sin(kstar*r/HCUT);
    return 0;
}
//_________________________________________________________________________
double waveUtils::usource(double *x,double *pm){
    double r = x[0];

    return source(x,pm) * r;
}
//_________________________________________________________________________
double waveUtils::source2int(double *x,double *pm){
    double r = x[0];

    return source(x,pm) * source(x,pm) * 4 * TMath::Pi() * r * r;
}
//_________________________________________________________________________
double waveUtils::sourceKin(double *x,double *pm){
  double r = x[0];
  double radius = pm[1];
  double kstar = pm[2];
  double R02inv = 1./(radius*radius);

  double val = 3 * R02inv * (HCUT*HCUT)/(2*MRED_NN);
  val -= r*r*R02inv*R02inv * (HCUT*HCUT)/(2*MRED_NN);
  val += kstar*kstar / (2*MRED_NN);

  return val * source2int(x,pm);
}
//_________________________________________________________________________
double waveUtils::sourceV(double *x,double *pm){
  double r = x[0];
  float chargeStrong = pm[2];
  float chargeColoumb = pm[3];
  double val = - (r < mStrongR) * mStrong * chargeStrong;
  val += chargeColoumb * mCoulomb / r;

  return val * source2int(x,pm);
}
//_________________________________________________________________________
void waveUtils::init(){
  mIsInitialized = true;
  mK1 = sqrt(2*MRED_NN*(EBOUND_D+mStrong))/HCUT;
  mK2 = sqrt(-2*MRED_NN*EBOUND_D)/HCUT;
  printf("Deuteron k1=%f - k2=%f - Radius=%f - V = %f (E=%f)\n",mK1,mK2,mStrongR,mStrong,EBOUND_D);
  mNormRight = sin(mK1*mStrongR) * TMath::Exp(mK2*mStrongR);

  mSource = new TF1("fSource",waveUtils::source,0,20,2);
  mUSource = new TF1("fUSource",waveUtils::usource,0,20,2);
  mSource2int = new TF1("fSource2int",waveUtils::source2int,0,20,2);
  mSourceKin = new TF1("fSourceKin",waveUtils::sourceKin,0,20,3);
  mSourceV = new TF1("fSourceV",waveUtils::sourceV,0,20,4);

  mDeuteron2int = new TF1("fDeuteron2int",waveUtils::intDeuteron,0,20,5);
  mDeuteron2int->SetParameter(0,1);
  mDeuteron2int->SetParameter(1,mStrongR);
  mDeuteron2int->SetParameter(2,mK1);
  mDeuteron2int->SetParameter(3,mK2);
  mDeuteron2int->SetParameter(4,mNormRight);
  double err;
  double norm = 1./sqrt(mDeuteron2int->IntegralOneDim(0,20,1E-8,1E-8,err));
  mDeuteron2int->SetParameter(0,norm);

  mDeuteron = new TF1("fDeuteron",waveUtils::waveDeuteron,0,20,5);
  mDeuteron->SetParameter(0,norm);
  mDeuteron->SetParameter(1,mStrongR);
  mDeuteron->SetParameter(2,mK1);
  mDeuteron->SetParameter(3,mK2);
  mDeuteron->SetParameter(4,mNormRight);

  mUDeuteron = new TF1("fUDeuteron",waveUtils::uDeuteron,0,20,5);
  mUDeuteron->SetParameter(0,norm);
  mUDeuteron->SetParameter(1,mStrongR);
  mUDeuteron->SetParameter(2,mK1);
  mUDeuteron->SetParameter(3,mK2);
  mUDeuteron->SetParameter(4,mNormRight);

  mDeuteronKin = new TF1("fDeuteronKin",waveUtils::deuteronKin,0,20,5);
  mDeuteronKin->SetParameter(0,norm);
  mDeuteronKin->SetParameter(1,mStrongR);
  mDeuteronKin->SetParameter(2,mK1);
  mDeuteronKin->SetParameter(3,mK2);
  mDeuteronKin->SetParameter(4,mNormRight);

  mDeuteronV = new TF1("fDeuteronV",waveUtils::deuteronV,0,20,7);
  mDeuteronV->SetParameter(0,norm);
  mDeuteronV->SetParameter(1,mStrongR);
  mDeuteronV->SetParameter(2,mK1);
  mDeuteronV->SetParameter(3,mK2);
  mDeuteronV->SetParameter(4,mNormRight);

  mCoalescenceRe = new TF1("fCoalescenceRe",waveUtils::coalescenceRe,0,20,8);
  mCoalescenceRe->SetParameter(3, norm);
  mCoalescenceRe->SetParameter(4, mStrongR);
  mCoalescenceRe->SetParameter(5, mK1);
  mCoalescenceRe->SetParameter(6, mK2);
  mCoalescenceRe->SetParameter(7, mNormRight);
  mCoalescenceIm = new TF1("fCoalescenceIm",waveUtils::coalescenceIm,0,20,8);
  mCoalescenceIm->SetParameter(3, norm);
  mCoalescenceIm->SetParameter(4, mStrongR);
  mCoalescenceIm->SetParameter(5, mK1);
  mCoalescenceIm->SetParameter(6, mK2);
  mCoalescenceIm->SetParameter(7, mNormRight);

  if(mSourceRadius >= 0){
    setSourceRadius(mSourceRadius);
  } else {
    setSourceRadius(0);
  }
}
//_________________________________________________________________________
float waveUtils::calcProb(){
  double err;
  float probRe = mCoalescenceRe->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err);
  probRe *= probRe;
  float probIm = mCoalescenceIm->IntegralOneDim(0,mMaxIntRange,1E-8,1E-8,err);
  probIm *= probIm;

  return (probRe + probIm) * mSpinCoalFactor;
}
//_________________________________________________________________________
float waveUtils::getCoalProb(const particleMC& p1, const particleMC& p2) {
  if(p1.pdg * p2.pdg < 0){ // baryon-antibaryon cannot do coalescence
    return 0;
  }
  if(p1.pdg == p2.pdg){ // identical particle cannot do coalescence
    return 0;
  }
  double kstar = utils::getKstar(p1,p2) * 1E3; // to MeV
  double kt = utils::getKt(p1,p2); // in GeV

  setKstar(kstar,kt);

  return calcProb();
}
