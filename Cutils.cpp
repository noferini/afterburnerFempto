#include "Cutils.h"

double utils::getMass(int pdg){
  double mass = 0;
  const double bindingDe = 2.22;
  const double bindingTr = 8.482;
  const double bindingHe3 = 7.714;
  const double bindingHe4 = 28.3;

  uint pdgAbs = std::abs(pdg);
  if(pdgAbs < 4000){
    mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  } else { // nuclei (as mapped from us) are not known
    int A = pdgAbs / 2000;
    int Z = (pdgAbs - A * 2100) / 100;
    int N = A - Z;
    mass = TDatabasePDG::Instance()->GetParticle(2212)->Mass() * Z;
    mass += TDatabasePDG::Instance()->GetParticle(2212)->Mass() * N;

    if(A == 2 && Z == 1){
      mass -= bindingDe;
    } else if(A == 3 && Z == 1) {
      mass -= bindingTr;
    } else if(A == 3 && Z == 2) {
      mass -= bindingHe3;
    } else if(A == 4 && Z == 2) {
      mass -= bindingHe4;
    }
  }
  return mass;
}
//_________________________________________________________________________
double utils::getKstar(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2){
  TLorentzVector pSum = p1.q[iPdg1] + p2.q[iPdg2];
  double Minv = pSum.M();
  double mass1square = p1.q[iPdg1].M()*p1.q[iPdg2].M();
  double mass2square = p2.q[iPdg1].M()*p2.q[iPdg2].M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKt(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2){
  TLorentzVector pSum = 0.5*(p1.q[iPdg1] + p2.q[iPdg2]);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getKstarAsPr(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2){
  float m1scaling = 0.938/p1.q[iPdg1].M();
  float m2scaling = 0.938/p2.q[iPdg2].M();

  TLorentzVector pSum = p1.q[iPdg1]*m1scaling + p2.q[iPdg2]*m2scaling;
  double Minv = pSum.M();
  double mass1square = p1.q[iPdg1].M()*p1.q[iPdg1].M()*m1scaling*m1scaling;
  double mass2square = p2.q[iPdg2].M()*p2.q[iPdg2].M()*m2scaling*m2scaling;
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKtAsPr(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2){
  float m1scaling = 0.938/p1.q[iPdg1].M();
  float m2scaling = 0.938/p2.q[iPdg2].M();
  TLorentzVector pSum = 0.5*(p1.q[iPdg1]*m1scaling + p2.q[iPdg2]*m2scaling);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getKstar(const particleMC& p1, const particleMC& p2){
  TLorentzVector pSum = p1.q + p2.q;
  double Minv = pSum.M();
  double mass1square = p1.q.M()*p1.q.M();
  double mass2square = p2.q.M()*p2.q.M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKt(const particleMC& p1, const particleMC& p2){
  TLorentzVector pSum = 0.5*(p1.q + p2.q);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getKstarAsPr(const particleMC& p1, const particleMC& p2){
  float m1scaling = 0.938/p1.q.M();
  float m2scaling = 0.938/p2.q.M();

  TLorentzVector pSum = p1.q*m1scaling + p2.q*m2scaling;
  double Minv = pSum.M();
  double mass1square = p1.q.M()*p1.q.M();
  double mass2square = p2.q.M()*p2.q.M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKtAsPr(const particleMC& p1, const particleMC& p2){
  float m1scaling = 0.938/p1.q.M();
  float m2scaling = 0.938/p2.q.M();
  TLorentzVector pSum = 0.5*(p1.q*m1scaling + p2.q*m2scaling);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getDPhi(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2, double rangeMin, double rangeMax, double shift){
  double phi1 = p1.q[iPdg1].Phi();
  double phi2 = p2.q[iPdg2].Phi();

  double dPhi = phi1 - phi2;
  if(dPhi < rangeMin) dPhi = dPhi + shift;
  if(dPhi > rangeMax) dPhi = dPhi - shift;

  return dPhi;
}
double utils::getDEta(const particleCand& p1, const particleCand& p2, int iPdg1, int iPdg2){
  double eta1 = p1.q[iPdg1].Eta();
  double eta2 = p2.q[iPdg2].Eta();

  double dEta = eta1 - eta2;
  return dEta;
}
//_________________________________________________________________________
void particleMC::print() const {
  printf("pdg = %d - mother = %d - charges: Strong = %d, Coloumb = %d - Ndaugthers = %lu\n",pdg, mother, StrongC, ColoumbC, daughters.size());
  q.Print();
}
//_________________________________________________________________________
void particleCand::addOption(double px, double py, double pz, int pdg){
  double mass = utils::getMass(pdg);

  double energy = sqrt(px*px + py*py + pz*pz + mass*mass);
  pdgOptions.push_back(pdg);
  q.emplace_back(px,py,pz,energy);
}
