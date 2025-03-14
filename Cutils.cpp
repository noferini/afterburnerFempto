#include "Cutils.h"

double utils::getKstar(const particleCand& p1, const particleCand& p2, int iPdg){
  TLorentzVector pSum = p1.q[iPdg] + p2.q[iPdg];
  double Minv = pSum.M();
  double mass1square = p1.q[iPdg].M()*p1.q[iPdg].M();
  double mass2square = p2.q[iPdg].M()*p2.q[iPdg].M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKt(const particleCand& p1, const particleCand& p2, int iPdg){
  TLorentzVector pSum = 0.5*(p1.q[iPdg] + p2.q[iPdg]);
  return pSum.Pt();
}
//_________________________________________________________________________
double utils::getKstarAsPr(const particleCand& p1, const particleCand& p2, int iPdg){
  float m1scaling = 0.938/p1.q[iPdg].M();
  float m2scaling = 0.938/p2.q[iPdg].M();

  TLorentzVector pSum = p1.q[iPdg]*m1scaling + p2.q[iPdg]*m2scaling;
  double Minv = pSum.M();
  double mass1square = p1.q[iPdg].M()*p1.q[iPdg].M();
  double mass2square = p2.q[iPdg].M()*p2.q[iPdg].M();
  double A = Minv*Minv - mass1square - mass2square;

  return sqrt(A*A - 4 * mass1square * mass2square) * 0.5 / Minv;
}
//_________________________________________________________________________
double utils::getKtAsPr(const particleCand& p1, const particleCand& p2, int iPdg){
  float m1scaling = 0.938/p1.q[iPdg].M();
  float m2scaling = 0.938/p2.q[iPdg].M();
  TLorentzVector pSum = 0.5*(p1.q[iPdg]*m1scaling + p2.q[iPdg]*m2scaling);
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
void particleMC::print() const {
  printf("pdg = %d - mother = %d - charges: Strong = %d, Coloumb = %d - Ndaugthers = %lu\n",pdg, mother, StrongC, ColoumbC, daughters.size());
  q.Print();
}
//_________________________________________________________________________
void particleCand::addOption(double px, double py, double pz, int pdg){
  double mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  double energy = sqrt(px*px + py*py + pz*pz + mass*mass);
  pdgOptions.push_back(pdg);
  q.emplace_back(px,py,pz,energy);
}
