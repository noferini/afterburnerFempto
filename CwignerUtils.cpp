#include "CwignerUtils.h"
#include "TMath.h"
#include "TF2.h"

double wignerUtils::mHCut = 0.1973; // GeV fm
double wignerUtils::mRMin = 0.;
double wignerUtils::mRMax = 20.;
double wignerUtils::mPMin = 0.;
double wignerUtils::mPMax = 0.6;
double wignerUtils::mDx = 0.01;
double wignerUtils::mDp = 0.001;
double wignerUtils::mFactor = sqrt(3. / 8) * mHCut;
float wignerUtils::mSpinCoalFactor = 1;

TF2 *wignerUtils::mW = nullptr;
TF2 *wignerUtils::mWxW = nullptr;
TF2 *wignerUtils::mWxJ = nullptr;
TF2 *wignerUtils::mWxJforItself = nullptr;
TF2 *wignerUtils::mK = nullptr;
TF2 *wignerUtils::mV = nullptr;
TF2 *wignerUtils::mH = nullptr;
TF2 *wignerUtils::mWK = nullptr;
TF2 *wignerUtils::mWV = nullptr;
TF2 *wignerUtils::mWH = nullptr;
TF2 *wignerUtils::mD = nullptr;
TF2 *wignerUtils::mDInt = nullptr;
TF2 *wignerUtils::mC = nullptr;

double wignerUtils::mR0 = 1.;
double wignerUtils::mRadius = 1.;
double wignerUtils::mKStar = 0.050;
double wignerUtils::mKin = 0.050;
double wignerUtils::mNorm = 1.;
double wignerUtils::mMu = utils::getMass(2212) / 2;
double wignerUtils::mRWidth = 3.2;  // 2.1
double wignerUtils::mV0 = 17.4E-3; //-0.0337
double wignerUtils::mBoundE = 2.2;  // MeV

TFile *wignerUtils::mFileDeuteron = new TFile("$cpath/deuteronFunction/wigner2.root", "READ");
TH2D *wignerUtils::mDeuteronH = (TH2D *)mFileDeuteron->Get("h");
//_________________________________________________________________________
void wignerUtils::init()
{
    mW = new TF2("w", wignerSource, mRMin, mRMax, mPMin, mPMax, 3);
    mWxJ = new TF2("wxj", jacobianFun, mRMin, mRMax, mPMin, mPMax, 3);
    mWxJforItself = new TF2("mWxJforItself", jacobianW2, mRMin, mRMax, mPMin, mPMax, 3);
    mWxW = new TF2("wxw", wignerSource2, mRMin, mRMax, mPMin, mPMax, 3);
    mK = new TF2("K", kineticEnergy, mRMin, mRMax, mPMin, mPMax, 4);
    mV = new TF2("V", potentialEnergy, mRMin, mRMax, mPMin, mPMax, 6);
    mH = new TF2("H", hamiltonian, mRMin, mRMax, mPMin, mPMax, 6);
    mWK = new TF2("WxK", wK, mRMin, mRMax, mPMin, mPMax, 4);
    mWV = new TF2("WxV", wV, mRMin, mRMax, mPMin, mPMax, 6);
    mWH = new TF2("WxH", wH, mRMin, mRMax, mPMin, mPMax, 6);
    mD = new TF2("WD", wignerDeuteron, mRMin, mRMax, mPMin, mPMax, 0);
    mDInt = new TF2("WDInt", wignerDeuteronIntegral, mRMin, mRMax, mPMin, mPMax, 0);
    mC = new TF2("CoalescenceProb", coalescenceProbability, mRMin, mRMax, mPMin, mPMax, 3);
    setFunctionsParameters();
    setIntegrationRanges(0., 20, 0., 0.6);
}
//_________________________________________________________________________
float wignerUtils::calcProb()
{
  return getcoal();
}
//_________________________________________________________________________
void wignerUtils::setCharges(float cS, float cC)
{
}
//_________________________________________________________________________
float wignerUtils::kineticDeuteron(){
  return 0;
}
//_________________________________________________________________________
float wignerUtils::potentialDeuteron(){
  return 0;
}//_________________________________________________________________________
float wignerUtils::getKstarFinal(float coalProb, float massRed, float boundE) {
  float kinetic = kineticSource();
  float potential = potentialSource();
  float Efin = TMath::Max(float(0.), (kinetic + potential - boundE * coalProb) / (1 - coalProb));
  return sqrt(Efin * 2 * massRed);
}
//_________________________________________________________________________
float wignerUtils::getCoalProb(const particleMC& p1, const particleMC& p2) {
  if (p1.pdg * p2.pdg < 0 || p1.pdg == p2.pdg){
    return 0;
  }

  if(p1.pdg == p2.pdg){ // identical particle cannot do coalescence
    return 0;
  }

  double kstar = utils::getKstar(p1, p2);
  double kt = utils::getKt(p1, p2);

  int pdgM = std::abs(p1.pdg);
  int pdgL = std::abs(p2.pdg);

  if(pdgM < pdgL){
    int dummy = pdgL;
    pdgL = pdgM;
    pdgM = dummy;
  }
  utils::type system=utils::none;

  if(pdgM == 2112){ // nn
    system = utils::nn;
  } else if (pdgM == 2212) { // pn or pp
    if(pdgL == 2112){
      system = utils::pn;
    } else {
      system = utils::pp;
    }
  } else if(pdgM == 4324) { // Dn or Dp or DD
    if(pdgL == 2112){
      system = utils::Dn;
    } else if(pdgL == 2212){
      system = utils::Dp;
    } else {
      system = utils::DD;
    }
  } else if(pdgM == 6436) { // Tn or Tp
    if(pdgL == 2112){
      system = utils::Tn;
    } else if(pdgL == 2212){
      system = utils::Tp;
    }
  } else if(pdgM == 6536) { // 3Hen o 3Hep
    if(pdgL == 2112){
      system = utils::Hen;
    } else if(pdgL == 2212){
      system = utils::Hep;
    }
  }

  if(system == utils::none){
    return 0;
  }

  setKstar(kstar, kt);
  return calcProb();
}
//_________________________________________________________________________
void wignerUtils::setFunctionsParameters()
{
    setThreeParam(mWxJ);
    normalization();
    setThreeParam(mW);
    setThreeParam(mWxJforItself);
    setThreeParam(mWxW);
    setThreeParam(mC);
    setFourParam(mK);
    setSixParam(mV);
    setSixParam(mH);
    setFourParam(mWK);
    setSixParam(mWV);
    setSixParam(mWH);
}
//_________________________________________________________________________
void wignerUtils::setRanges(double xmin, double ymin, double xmax, double ymax)
{
    if (xmin < wignerUtils::getMinX() || xmax > wignerUtils::getMaxX() || ymin < wignerUtils::getMinP() || ymax > wignerUtils::getMaxP())
    {
        std::cout << "Invalid ranges specified; previous ranges kept.\n";
    }
    else
    {
        mW->SetRange(xmin, ymin, xmax, ymax);
        mRMin = xmin;
        mRMax = xmax;
        mPMin = ymin;
        mPMax = ymax;
    }
}
//_________________________________________________________________________
void wignerUtils::setMu(double mu)
{
    mMu = mu;
    reSetMu();
}
//_________________________________________________________________________
void wignerUtils::setRWidth(double rWidth)
{
    mRWidth = rWidth;
    reSetRWidth();
}
//_________________________________________________________________________
void wignerUtils::setV0(double v0)
{
    mV0 = v0;
    reSetV0();
}
//_________________________________________________________________________
void wignerUtils::setThreeParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
}
//_________________________________________________________________________
void wignerUtils::setFourParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerUtils::setSixParam(TF2 *function)
{
    function->SetParameter(0, mNorm);
    function->SetParameter(1, mRadius);
    function->SetParameter(2, mKStar);
    function->SetParameter(3, mMu);
    function->SetParameter(4, mRWidth);
    function->SetParameter(5, mV0);
}
//_________________________________________________________________________
double wignerUtils::wignerSource(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double norm = 1. / (TMath::Pi() * mHCut);
    norm *= norm * norm;

    return pm[0] * norm * TMath::Exp(-r * r * 0.25 / (pm[1] * pm[1]) - 4 * (p * p + pm[2] * pm[2] - 2 * p * pm[2]) * (pm[1] * pm[1]) / (mHCut * mHCut));
}
//_________________________________________________________________________
double wignerUtils::wignerSource2(double *x, double *pm)
{
    return wignerSource(x, pm) * jacobianW2(x, pm);
}
//_________________________________________________________________________
double wignerUtils::jacobianFun(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 8 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);

    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;

    return wignerSource(x, pm) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::jacobianW2(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 16 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);
    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;
    return jacobian * wignerSource(x, pm);
}
//_________________________________________________________________________
double wignerUtils::kineticEnergy(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return (p * p) / (2 * pm[3]);
}
//_________________________________________________________________________
double wignerUtils::potentialEnergy(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double E = 0;
    if (r < pm[4])
    {
        E = pm[5];
    }
    return E;
}
//_________________________________________________________________________
double wignerUtils::hamiltonian(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return kineticEnergy(x, pm) + potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wK(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * kineticEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wV(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * potentialEnergy(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wH(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return jacobianFun(x, pm) * hamiltonian(x, pm);
}
//_________________________________________________________________________
double wignerUtils::wignerDeuteron(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    return mDeuteronH->Interpolate(r, p);
}
//_________________________________________________________________________
double wignerUtils::wignerDeuteronIntegral(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = 16 * TMath::Pi() * TMath::Pi() * r * r * p * p;
    return mDeuteronH->Interpolate(r, p) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::coalescenceProbability(double *x, double *pm)
{
    double r = x[0];
    double p = x[1];

    double jacobian = (r * p) * (r * p) * 16 * TMath::Pi() * TMath::Pi();

    float kstarP = pm[2] * p;
    if (kstarP < 1E-16)
    {
        kstarP = 1E-16;
    }
    double alpha = 8 * kstarP * pm[1] * pm[1] / (mHCut * mHCut);

    jacobian *= 0.5 * (1 - TMath::Exp(-2 * alpha)) / alpha;

    return wignerDeuteron(x, pm) * wignerSource(x, pm) * jacobian;
}
//_________________________________________________________________________
double wignerUtils::radius(double k, double r0)
{
    double radiusWave = mFactor / k;
    return sqrt(radiusWave * radiusWave + r0 * r0);
}
//_________________________________________________________________________
double wignerUtils::kStarEff(double k, double radius)
{
    float radiusWave = mFactor/k;
    mRadius = sqrt(radiusWave*radiusWave + radius*radius);
    float kstarWave = mFactor/mRadius;
    float kstarEff = sqrt(k*k - kstarWave*kstarWave);
    return kstarEff; //radius = sqrt(k * k - kstarWave * kstarWave);
}
//_________________________________________________________________________
double wignerUtils::getMinX()
{
    return mRMin;
}
//_________________________________________________________________________
double wignerUtils::getMaxX()
{
    return mRMax;
}
//_________________________________________________________________________
double wignerUtils::getMinP()
{
    return mPMin;
}
//_________________________________________________________________________
double wignerUtils::getHCut()
{
    return mHCut;
}
//_________________________________________________________________________
double wignerUtils::getMaxP()
{
    return mPMax;
}
//_________________________________________________________________________
void wignerUtils::setMinX(double minX)
{
    mRMin = minX;
}
//_________________________________________________________________________
void wignerUtils::setMaxX(double maxX)
{
    mRMax = maxX;
}
//_________________________________________________________________________
void wignerUtils::setMinP(double minP)
{
    mPMin = minP;
}
//_________________________________________________________________________
void wignerUtils::setMaxP(double maxP)
{
    mPMax = maxP;
}
//_________________________________________________________________________
void wignerUtils::setIntegrationRanges(double minX, double maxX, double minP, double maxP)
{
    setMinX(minX);
    setMaxX(maxX);
    setMinP(minP);
    setMaxP(maxP);
}
//_________________________________________________________________________
double wignerUtils::integral(TF2 *function, double minX, double maxX, double minP, double maxP)
{
    double res = 0;
    // res = function->Integral(minX, maxX, minP, maxP);
    for (float x = mDx / 2 + minX; x < maxX; x += mDx)
    {
        for (float p = mDp / 2 + minP; p < maxP; p += mDp)
        {
            res += function->Eval(x, p);
        }
    }
    res *= mDx * mDp;
    return res;
}
//_________________________________________________________________________
void wignerUtils::normalization()
{
    double norm = 1. / wignerUtils::integral(mWxJ, 0., TMath::Max(5. * mRadius, 20.), 0., 0.6);
    mNorm = norm;
}
//_________________________________________________________________________
void wignerUtils::setParams(float strong, float strongR, float coloumb, float sourceRadius, float spinFact)
{
    mV0 = strong;
    mRWidth = strongR;
    mR0 = std::abs(sourceRadius);
    mSpinCoalFactor = spinFact;
}
//_________________________________________________________________________
void wignerUtils::setSourceRadius(float radius)
{
//    mR0 = std::abs(radius);
}
//_________________________________________________________________________
void wignerUtils::setKstar(float kstar, float kt, utils::type system)
{
    setRadiusK(double(kstar));
/*
    mKStar = kStarEff(kstar, mR0);
    mWxJ->SetParameter(0, 1.);
    reSetKStar();
    normalization();
    reSetNorm();
*/
}
//_________________________________________________________________________
void wignerUtils::setRadiusK(double k)
{
    mKin = k;
//    mRadius = wignerUtils::radius(k, mR0);
    mKStar = wignerUtils::kStarEff(k, mR0);
    mWxJ->SetParameter(0, 1.);
    reSetRadius();
    reSetKStar();
    normalization();
    reSetNorm();
}
//_________________________________________________________________________
void wignerUtils::reSetKStar()
{
    mW->SetParameter(2, mKStar);
    mWxJ->SetParameter(2, mKStar);
    mWxJforItself->SetParameter(2, mKStar);
    mWxW->SetParameter(2, mKStar);
    mK->SetParameter(2, mKStar);
    mV->SetParameter(2, mKStar);
    mH->SetParameter(2, mKStar);
    mWK->SetParameter(2, mKStar);
    mWV->SetParameter(2, mKStar);
    mWH->SetParameter(2, mKStar);
    mC->SetParameter(2, mKStar);
}
//_________________________________________________________________________
void wignerUtils::reSetNorm()
{
    mW->SetParameter(0, mNorm);
    mWxJ->SetParameter(0, mNorm);
    mWxJforItself->SetParameter(0, mNorm);
    mWxW->SetParameter(0, mNorm);
    mK->SetParameter(0, mNorm);
    mV->SetParameter(0, mNorm);
    mH->SetParameter(0, mNorm);
    mWK->SetParameter(0, mNorm);
    mWV->SetParameter(0, mNorm);
    mWH->SetParameter(0, mNorm);
    mC->SetParameter(0, mNorm);
}
//_________________________________________________________________________
void wignerUtils::reSetMu()
{
    mK->SetParameter(3, mMu);
    mV->SetParameter(3, mMu);
    mH->SetParameter(3, mMu);
    mWK->SetParameter(3, mMu);
    mWV->SetParameter(3, mMu);
    mWH->SetParameter(3, mMu);
}
//_________________________________________________________________________
void wignerUtils::reSetRWidth()
{
    mK->SetParameter(4, mRWidth);
    mV->SetParameter(4, mRWidth);
    mH->SetParameter(4, mRWidth);
    mWK->SetParameter(4, mRWidth);
    mWV->SetParameter(4, mRWidth);
    mWH->SetParameter(4, mRWidth);
}
//_________________________________________________________________________
void wignerUtils::reSetV0()
{
    mK->SetParameter(5, mV0);
    mV->SetParameter(5, mV0);
    mH->SetParameter(5, mV0);
    mWK->SetParameter(5, mV0);
    mWV->SetParameter(5, mV0);
    mWH->SetParameter(5, mV0);
}
//_________________________________________________________________________
void wignerUtils::reSetRadius()
{
    mW->SetParameter(1, mRadius);
    mWxJ->SetParameter(1, mRadius);
    mWxJforItself->SetParameter(1, mRadius);
    mWxW->SetParameter(1, mRadius);
    mK->SetParameter(1, mRadius);
    mV->SetParameter(1, mRadius);
    mH->SetParameter(1, mRadius);
    mWK->SetParameter(1, mRadius);
    mWV->SetParameter(1, mRadius);
    mWH->SetParameter(1, mRadius);
    mC->SetParameter(1, mRadius);
}
//_________________________________________________________________________
float wignerUtils::kineticSource()
{
    return integral(mWK);
}
//_________________________________________________________________________
float wignerUtils::potentialSource()
{
    return -integral(mWV);
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerFunction()
{
    return mW;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerFunctionForItself()
{
    return mWxW;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerFunctionForJacobian()
{
    return mWxJ;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerFunction2ForJacobian()
{
    return mWxJforItself;
}
//_________________________________________________________________________
TF2 *wignerUtils::getKineticEnergyFunction()
{
    return mK;
}
//_________________________________________________________________________
TF2 *wignerUtils::getPotentialEnergyFunction()
{
    return mV;
}
//_________________________________________________________________________
TF2 *wignerUtils::getHamiltonianFunction()
{
    return mH;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerKinetic()
{
    return mWK;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerPotential()
{
    return mWV;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerHamiltonan()
{
    return mWH;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerDeuteron()
{
    return mD;
}
//_________________________________________________________________________
TF2 *wignerUtils::getWignerDeuteronIntegral()
{
    return mDInt;
}
//_________________________________________________________________________
TF2 *wignerUtils::getCoalescenceProbability()
{
    return mC;
}
//_________________________________________________________________________
double wignerUtils::getNorm()
{
    return mNorm;
}
//_________________________________________________________________________
double wignerUtils::getWH()
{
    return integral(mWH);
}
//_________________________________________________________________________
double wignerUtils::checkWxW()
{
    return integral(mWxW) * (getHCut() * 2 * TMath::Pi()) * (getHCut() * 2 * TMath::Pi()) * (getHCut() * 2 * TMath::Pi());
}
//_________________________________________________________________________
double wignerUtils::getRadius()
{
    return mRadius;
}
//_________________________________________________________________________
double wignerUtils::getKStar()
{
    return mKStar;
}
//_________________________________________________________________________
double wignerUtils::getRMin()
{
    return mRMin;
}
//_________________________________________________________________________
double wignerUtils::getRMax()
{
    return mRMax;
}
//_________________________________________________________________________
double wignerUtils::getPMin()
{
    return mPMin;
}
//_________________________________________________________________________
double wignerUtils::getPMax()
{
    return mPMax;
}
//_________________________________________________________________________
double wignerUtils::getMu()
{
    return mMu;
}
//_________________________________________________________________________
double wignerUtils::getRWidth()
{
    return mRWidth;
}
//_________________________________________________________________________
double wignerUtils::getV0()
{
    return mV0;
}
//_________________________________________________________________________
double wignerUtils::getcoal()
{
//    printf("R=%f - k*=%f\n",mKStar,mRadius);

    return integral(mC) * (getHCut() * 2 * TMath::Pi()) * (getHCut() * 2 * TMath::Pi()) * (getHCut() * 2 * TMath::Pi()) * mSpinCoalFactor;
}
//_________________________________________________________________________
double wignerUtils::getDeuteronInt()
{
    return integral(mDInt);
}
