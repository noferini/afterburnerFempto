#include "Cvfempto.h"
#include "TRandom.h"

void vfempto::doInteractAll(std::vector<particleMC>& part, bool doScattering, bool doCoal){
  unsigned long npart = part.size();

  doInteractAllStep(0, part, doScattering, doCoal);

  for(int i=1; i < mNmaxSteps; i++){
    if(part.size() > npart){ // go ahead only if new particles were produced via coalescence
      npart = part.size();
      doInteractAllStep(i, part, false, doCoal);
    }
  }
}
//_________________________________________________________________________
void vfempto::doInteractAllStep(int step, std::vector<particleMC>& part, bool doScattering, bool doCoal){
  if(! mIsInitialized){
    init();
  }

  // build pairs only if both charged
  std::vector<int> groupBelongTo;
  std::vector<std::pair<int,int>> intPairs;
  for(int i1 = 0; i1 < part.size(); i1++){
    groupBelongTo.push_back(-1); // inizializing to none
    int sC = part[i1].StrongC;
    int cC = part[i1].ColoumbC;

    if(part[i1].mother > -1){ // only mother inteacts
      continue;
    }

    if(!(sC || cC)){ // particle can interact otherwise skip
      continue;
    }

    for(int i2 = i1+1; i2 < part.size(); i2++){
      if(part[i2].mother > -1){ // only mother inteacts
        continue;
      }

      if(!(part[i2].StrongC && sC || part[i2].ColoumbC && cC)){ // two particles can interact otherwise skip
        continue;
      }

      double kstar = utils::getKstar(part[i1],part[i2]);
      if(kstar > mThreshold){
        continue;
      }

      intPairs.push_back(std::make_pair(i1,i2));
    }
  }

  // sorting by kstar
  std::sort(intPairs.begin(), intPairs.end(), [part](const auto &a, const auto &b)
  {
     double ks1 = utils::getKstar(part[a.first],part[a.second]);
     double ks2 = utils::getKstar(part[b.first],part[b.second]);
     return ks1 < ks2;
  });

  std::vector<std::vector<int>> intPairsGrouped;
  std::vector<std::vector<int>> intPairsGroupedCoal;
  for(const auto& o : intPairs){ // do interactions starting from lower kstar
    int iP1 = o.first;
    int iP2 = o.second;

    if(groupBelongTo[iP1] < 0 && groupBelongTo[iP2] < 0){ // not yet grouped
      int newGroup = intPairsGrouped.size();
      std::vector<int> grp;
      grp.push_back(iP1);
      grp.push_back(iP2);
      std::vector<int> coal;
      coal.push_back(-1);
      coal.push_back(-1);
      intPairsGrouped.push_back(grp);
      intPairsGroupedCoal.push_back(coal);
      groupBelongTo[iP1] = newGroup;
      groupBelongTo[iP2] = newGroup;
    } else if(groupBelongTo[iP1] < 0) {
      int group = groupBelongTo[iP2];
      intPairsGrouped[group].push_back(iP1);
      intPairsGroupedCoal[group].push_back(-1);
      groupBelongTo[iP1] = group;
    } else if(groupBelongTo[iP2] < 0) {
      int group = groupBelongTo[iP1];
      intPairsGrouped[group].push_back(iP2);
      intPairsGroupedCoal[group].push_back(-1);
      groupBelongTo[iP2] = group;
    } else if(groupBelongTo[iP1] != groupBelongTo[iP2]){ // two groups have to be merged
      int ifirst = groupBelongTo[iP1];
      int isecond = groupBelongTo[iP2];
      for(auto& o : intPairsGrouped[isecond]){
        intPairsGrouped[ifirst].push_back(o);
        intPairsGroupedCoal[ifirst].push_back(-1);
        groupBelongTo[o] = ifirst;
      }
      intPairsGrouped[isecond].clear();
      intPairsGroupedCoal[isecond].clear();
    }
  }

  if(! mHsizeGroup){
    mHsizeGroup = new TH1F("mHsizeGroup","",18,1,19);
    mHsizeMerge = new TH1F("mHsizeMerge","",18,1,19);
  }

  for(int ip = 0; ip < intPairsGrouped.size(); ip++){
    int nclusters = 0;
    const auto& o = intPairsGrouped[ip];
    auto& o2 = intPairsGroupedCoal[ip];
    mHsizeGroup->Fill(o.size());
    for(int i1 = 0; i1 < o.size(); i1++){
      particleMC& original1 = part[o[i1]];
      if(original1.daughters.size()){
        continue;
      }
      for(int i2 = i1+1; i2 < o.size(); i2++){
        particleMC& original2 = part[o[i2]];
        if(original2.daughters.size()){
          continue;
        }

/*
        if(o.size() > 8){
          double ks = utils::getKstar(original1,original2);
          printf("%d/%d) k* = %f\n",i1,i2,ks);
        }
*/

        float chargeCoulomb = original1.ColoumbC*original2.ColoumbC;
        float chargeStrong = original1.StrongC*original2.StrongC;
        if(! set(original1,original2, chargeCoulomb, chargeStrong)){
          continue;
        }
        // check and mark for coalescence
        if(o2[i1] < 0 && o2[i2] < 0 && isCoalescence(original1,original2,step)){
          if(o2[i1] < 0 && o2[i2] < 0) { // new cluster
            o2[i1] = nclusters;
            o2[i2] = nclusters;
            nclusters++;
          } else if(o2[i1] < 0) { // add to an already existing cluster
            o2[i1] = o2[i2];
          } else if(o2[i2] < 0) { // add to an already existing cluster
            o2[i2] = o2[i1];
          } else if(o2[i1] != o2[i2]) { // merge clusters
            for(int i3 = 0; i3 < o.size(); i3++){
              if(o2[i3] == o2[i2]) {
                o2[i3] = o2[i1];
              }
            }
          }
        }
        // ------------------------------

        // interact
        if(doScattering){
          doInteract(original1,original2,chargeCoulomb,chargeStrong,0);
        }
        // ------------------------------
      }
    }

    if(! doCoal){
      continue;
    }

    // do coalescence for marked clusters
    std::vector<int> mergeable;
    for(int ic=0; ic < nclusters; ic++){
      mergeable.clear();
      for(int i1 = 0; i1 < o.size(); i1++){
        if(o2[i1] == ic){
          mergeable.push_back(o[i1]);
        }
      }
      // merge all particle in mergeable (if nucleus is allowed otherwise keep the maximum allowed. In ambiguous cases let's decide randomly)

      int newpart = part.size();
      int npro = 0;
      int nneu = 0;
      int ndeu = 0;
      int ntri = 0;
      int nhe = 0;

      for(const auto& o : mergeable){ // count nucleons
        particleMC& p = part[o];
        if(std::abs(p.pdg) == 2212) {
          npro++;
        }
        if(std::abs(p.pdg) == 2112) {
          nneu++;
        }
        if(std::abs(p.pdg) == 4324) {
          ndeu++;
        }
        if(std::abs(p.pdg) == 6436) {
          ntri++;
        }
        if(std::abs(p.pdg) == 6536) {
          nhe++;
        }
      }
      int A = npro + nneu + ndeu*2 + ntri*3 + nhe*3;
      if(1){//A == 2 && npro == 1){ // detueron
        particleMC& p1 = part[mergeable[0]];
        particleMC& p2 = part[mergeable[1]];

        particleMC photon;
        part.push_back(photon);
        p1.daughters.push_back(newpart);
        p2.daughters.push_back(newpart);

        particleMC merged = merge(p1,p2,part[newpart]);

        newpart++;
        p1.daughters.push_back(newpart);
        p2.daughters.push_back(newpart);
        part.push_back(merged);

        mHsizeMerge->Fill(mergeable.size());
      }
    }
  }
}
//_________________________________________________________________________
float vfempto::getCoalProb(const particleMC& p1, const particleMC& p2) {
  double kstar = utils::getKstar(p1,p2);
  return TMath::Exp(-kstar/0.05);
}
//_________________________________________________________________________
bool vfempto::isCoalescence(const particleMC& p1, const particleMC& p2, int step) {
  if(p1.pdg * p2.pdg < 0){ // baryon + antibaryon
    return false;
  }

  int pdgMax = std::abs(p1.pdg);
  if(std::abs(p2.pdg) > pdgMax){
    pdgMax = std::abs(p2.pdg);
  }

  if(step > 0 && pdgMax < 3000){ // baryon-baryon allowed only at step 0
    return false;
  }
  if(step > 1 && pdgMax < 5000){ // D-baryon allowed only at step 1
    return false;
  }
  if(step > 3){ // only 3 steps foreseen
    return false;
  }

  float prob = getCoalProb(p1, p2);

  if(gRandom->Rndm() < prob){
    return true;
  }

  return false;
}
//_________________________________________________________________________
double vfempto::doInteract(particleMC& p1, particleMC& p2, float chargeColoumb, float chargeStrong, float sumRadii, float *pos, float *posLab){
  double kstar = utils::getKstar(p1,p2);

  float *lpos,*lposLab; // 3+3 position vectors

  if(pos){ // return also position in the output (working on it)
    lpos = pos;
  } else { // working on a local variable
    lpos = new float[8];
  }
  if(posLab){ // return also lab position in the output (working on it)
    lposLab = posLab;
  } else { // working on a local variable
    lposLab = new float[8];
  }

  TVector3 impactparam(0,0,0);

  float sumRadiiOr = sumRadii;

  sumRadii = 0.1973 * 0.5 / kstar;

  if(sumRadii < sumRadiiOr){
    sumRadii = sumRadiiOr;
  }

  float extraSmear = gRandom->Gaus(0,mSourceRadius);

  lposLab[0] = sumRadii * 0.5;
  lposLab[4] = -lposLab[0];
  lposLab[1] = extraSmear*0.5;
  lposLab[5] = -lposLab[1];
  lposLab[2] = lposLab[3] = lposLab[6] = lposLab[7] = 0;

  double scaling;

  TLorentzVector pSum = p1.q + p2.q;
  TVector3 b = pSum.BoostVector();
  TVector3 bInv = -b;

  p1.q.Boost(bInv);
  p2.q.Boost(bInv);

  double beta1 = p1.q.Beta();
  double beta2 = p2.q.Beta();

  TLorentzVector v1(lposLab[0],lposLab[1],lposLab[2],lposLab[3]);
  v1.Boost(bInv);
  TLorentzVector v2(lposLab[4],lposLab[5],lposLab[6],lposLab[7]);
  v2.Boost(bInv);
  lpos[0]=lposLab[0];
  lpos[1]=lposLab[1];
  lpos[2]=lposLab[2];
  lpos[3]=lposLab[3];
  lpos[4]=lposLab[4];
  lpos[5]=lposLab[5];
  lpos[6]=lposLab[6];
  lpos[7]=lposLab[7];
  impactparam.SetXYZ(lpos[0]-lpos[4], lpos[1]-lpos[5], lpos[2]-lpos[6]);

  impactparam.SetXYZ(lpos[0]-lpos[4], lpos[1]-lpos[5], lpos[2]-lpos[6]);

  float dist = sqrt(impactparam.Mag2());
  float deltaE = 0;
  deltaE = chargeColoumb * mCoulomb / dist;         // Columb contribution (repulsive/attractive)
  deltaE -= chargeStrong*mStrong*(dist < mStrongR); // Strong potential attractive

  double E0 = p1.q.Energy()+p2.q.Energy();
  double M0 = p1.q.M() + p2.q.M();
  double EK0 = E0 - M0;
  double EKF = EK0 + deltaE;

  if(EKF < 0){
    EKF = 0;
  }
  double EF = M0 + EKF;

  double momFinal = sqrt(EF*EF - M0*M0)*0.5;
  scaling = momFinal / p1.q.P();

  double m1 = p1.q.M();
  double px1= p1.q.Px()*scaling;
  double py1= p1.q.Py()*scaling;
  double pz1= p1.q.Pz()*scaling;
  double e1 = sqrt(px1*px1 + py1*py1 + pz1*pz1 + m1*m1);
  double ptEx = std::abs(momFinal) - std::abs(p1.q.P());
  p1.q.SetPxPyPzE(px1,py1,pz1,e1);
  double m2 = p2.q.M();
  double px2= p2.q.Px()*scaling;
  double py2= p2.q.Py()*scaling;
  double pz2= p2.q.Pz()*scaling;
  double e2 = sqrt(px2*px2 + py2*py2 + pz2*pz2 + m2*m2);
  p2.q.SetPxPyPzE(px2,py2,pz2,e2);

  p1.q.Boost(b);
  p2.q.Boost(b);

  if(!pos){ // clean up local
    delete lpos;
  }
  if(!posLab){ // clean up local
    delete lposLab;
  }

  return ptEx;
}
//_________________________________________________________________________
particleMC vfempto::merge(const particleMC& p1, const particleMC& p2, particleMC& photon){
  particleMC pSum;
  pSum.q = p1.q + p2.q;
  int pdgSum = p1.pdg + p2.pdg;
  double mass = utils::getMass(pdgSum);
  double DeltaE = pSum.q.M() - mass;                      // intial energy - mass consitutents + binding energy
  double energyPhotonRestFrame = DeltaE * (1 + DeltaE*0.5/mass) / (1 + DeltaE/mass);
  double p = pSum.q.P();
  // photon in rest frame
  double theta = acos(1 - 2*gRandom->Rndm());
  double phi = gRandom->Rndm()*2*TMath::Pi();
  photon.q.SetPxPyPzE(energyPhotonRestFrame*sin(theta)*cos(phi),energyPhotonRestFrame*sin(theta)*sin(phi),energyPhotonRestFrame*cos(theta),energyPhotonRestFrame);
  TVector3 b = pSum.q.BoostVector();

  // final detuerons
  pSum.q.SetPxPyPzE(-photon.q.Px(),-photon.q.Py(),-photon.q.Pz(),sqrt(photon.q.P()*photon.q.P() + mass*mass));
  pSum.q.Boost(b); //to lab

  pSum.pdg = p1.pdg + p2.pdg;
  pSum.mother = -1;              // not tracing mothers
  pSum.StrongC = p1.StrongC + p2.StrongC;
  pSum.ColoumbC = p1.ColoumbC + p2.ColoumbC;

  // final photon to lab
  photon.q.Boost(b);
  photon.pdg = 22;

  return pSum;
}
