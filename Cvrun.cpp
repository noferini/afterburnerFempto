#include "Cvrun.h"

void vrun::setEvent(const std::vector<particleMC>& vect){
  mVect.clear();
  for(const auto& o : vect){
    particleCand p;
    p.addOption(o.q.Px(), o.q.Py(), o.q.Pz(), o.pdg); // in this case only one option is set since MC truth is known

    mVect.push_back(p);
  }
}
//_________________________________________________________________________
void vrun::init(){
  initHistos();
}
//_________________________________________________________________________
void vrun::initHistos(){
    mHkstarSE = new TH2D("mHkstarSE","Same Event;k_{T} (GeV/c);k* (MeV/c);N(k*)",10,0,2,100,0,1);
    mHkstarME = new TH2D("mHkstarME","Mixed Events;k_{T} (GeV/c);k* (MeV/c);N(k*)",10,0,2,100,0,1);
}
//_________________________________________________________________________
void vrun::finalize(){
  finalizeHistos();
}
//_________________________________________________________________________
void vrun::finalizeHistos(){
}
//_________________________________________________________________________
void vrun::write(){
  mHkstarSE->Write();
  mHkstarME->Write();
}
//_________________________________________________________________________
void vrun::process(){
  static int nevDone=0; // start to process mixed events only when nevDone >= mNmix

  int lev = nevDone % mNmix;  // event-ID in the current cycle (cycle length defined my mNmix)


  for(int i1=0; i1 < mVect.size(); i1++){ // same event
    const particleCand& p1 = mVect[i1];

    int ipdg1 = selectP1(p1);

    if(ipdg1 < 0){
      continue;
    }

    for(int i2=(i1+1)*mIsSamePart; i2 < mVect.size(); i2++){
      const particleCand& p2 = mVect[i2];

      if(i1 == i2){
        continue;
      }

      int ipdg2 = selectP2(p2);

      if(ipdg2 < 0){
        continue;
      }
      double kt = utils::getKt(p1,p2,ipdg1);
      double kstar = utils::getKstar(p1,p2,ipdg2);

      mHkstarSE->Fill(kt, kstar);
    }
  }

  if(nevDone >= mNmix){ // mix events
    for(int iM=0; iM < mNmix; iM++){

      for(int i1=0; i1 < mVect.size(); i1++){
        const particleCand& p1 = mVect[i1];

        int ipdg1 = selectP1(p1);

        if(ipdg1 < 0){
          continue;
        }

        for(int i2=(i1+1)*mIsSamePart; i2 < mEvPrev[iM].size(); i2++){
           const particleCand& p2 = mEvPrev[iM][i2];

           int ipdg2 = selectP2(p2);

           if(ipdg2 < 0){
             continue;
           }

           double kt = utils::getKt(p1,p2,ipdg1);
           double kstar = utils::getKstar(p1,p2,ipdg2);

           mHkstarME->Fill(kt, kstar);
        }
      }
    }
  }

  // update array of mixed events
  mEvPrev[lev] = mVect;

  nevDone++;
}
//_________________________________________________________________________
int vrun::selectP1(const particleCand& p){
  for(int i=0; i < p.pdgOptions.size(); i++){
    if(p.pdgOptions[i] == mPDG1){
      return i;
    }
  }
  return -1;
}
//_________________________________________________________________________
int vrun::selectP2(const particleCand& p){
  for(int i=0; i < p.pdgOptions.size(); i++){
    if(p.pdgOptions[i] == mPDG2){
      return i;
    }
  }
  return -1;
}
