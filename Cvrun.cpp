#include "Cvrun.h"
#include "TFile.h"

void vrun::setEvent(const std::vector<particleMC>& vect){
  mVect.clear();
  for(const auto& o : vect){
    particleCand p;
    p.addOption(o.q.Px(), o.q.Py(), o.q.Pz(), o.pdg); // in this case only one option is set since MC truth is known

    mVect.push_back(p);
  }
}
void vrun::setIDName(const TString& name){
  mName = name;
}
//_________________________________________________________________________
void vrun::init(){
  initHistos();
}
//_________________________________________________________________________
void vrun::initHistos(){
    mHkstarSE = new TH2D("mHkstarSE" + mName,"Same Event;k_{T} (GeV/c);k* (MeV/c);N(k*)",10,0,2,100,0,1);
    mHkstarME = new TH2D("mHkstarME" + mName,"Mixed Events;k_{T} (GeV/c);k* (MeV/c);N(k*)",10,0,2,100,0,1);
    mDPhiDEtaSE = new TH2D("mDPhiDEtaSE" + mName,"Same Event;#Delta#eta;#Delta#Phi;N(#Delta#Phi)",40,-2,2,20,-TMath::Pi()*0.5,TMath::Pi()*3*0.5);
    mDPhiDEtaME = new TH2D("mDPhiDEtaME" + mName,"Mixed Event;#Delta#eta;#Delta#Phi;N(#Delta#Phi)",40,-2,2,20,-TMath::Pi()*0.5,TMath::Pi()*3*0.5);
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
  gFile->mkdir("vrun" + mName);
  gFile->cd("vrun" + mName);
  mHkstarSE->Write();
  mHkstarME->Write();
  mDPhiDEtaSE->Write();
  mDPhiDEtaME->Write();
  gFile->cd();
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
      double kt = utils::getKt(p1,p2,ipdg1,ipdg2);
      double kstar = utils::getKstar(p1,p2,ipdg1,ipdg2);

      bool conditionEta1 = (p1.q[ipdg1].Eta() > -1. && p1.q[ipdg1].Eta() < 1.);
      bool conditionEta2 = (p2.q[ipdg2].Eta() > -1. && p2.q[ipdg2].Eta() < 1.);
      bool conditionPt1 = (p1.q[ipdg1].Pt() > 0.4 && p1.q[ipdg1].Pt() < 1.);
      bool conditionPt2 = (p2.q[ipdg2].Pt() > 0.4 && p2.q[ipdg2].Pt() < 1.);

      if(conditionEta1 && conditionEta2 && conditionPt1 && conditionPt2){
        double rangeMin = -TMath::Pi()*0.5;
        double rangeMax = TMath::Pi()*3*0.5;
        double shift = TMath::Pi()*2;
        double dPhi = utils::getDPhi(p1,p2,ipdg1,ipdg2,rangeMin,rangeMax,shift);
        double dEta = utils::getDEta(p1,p2,ipdg1,ipdg2);

        mDPhiDEtaSE->Fill(dEta, dPhi);
      }

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

           double kt = utils::getKt(p1,p2,ipdg1,ipdg2);
           double kstar = utils::getKstar(p1,p2,ipdg1,ipdg2);

           bool conditionEta1 = (p1.q[ipdg1].Eta() > -1. && p1.q[ipdg1].Eta() < 1.);
           bool conditionEta2 = (p2.q[ipdg2].Eta() > -1. && p2.q[ipdg2].Eta() < 1.);
           bool conditionPt1 = (p1.q[ipdg1].Pt() > 0.4 && p1.q[ipdg1].Pt() < 1.);
           bool conditionPt2 = (p2.q[ipdg2].Pt() > 0.4 && p2.q[ipdg2].Pt() < 1.);

           if(conditionEta1 && conditionEta2 && conditionPt1 && conditionPt2){
            double rangeMin = -TMath::Pi()*0.5;
            double rangeMax = TMath::Pi()*3*0.5;
            double shift = TMath::Pi()*2;
            double dPhi = utils::getDPhi(p1,p2,ipdg1,ipdg2,rangeMin,rangeMax,shift);
            double dEta = utils::getDEta(p1,p2,ipdg1,ipdg2);

            mDPhiDEtaME->Fill(dEta, dPhi);
           }

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
//_________________________________________________________________________
TH2D* vrun::getKstarSE(){
  return mHkstarSE;
}
//_________________________________________________________________________
TH2D* vrun::getKstarME(){
  return mHkstarME;
}
//_________________________________________________________________________
TH2D* vrun::getDPhiDEtaSE(){
  return mDPhiDEtaSE;
}
//_________________________________________________________________________
TH2D* vrun::getDPhiDEtaME(){
  return mDPhiDEtaME;
}

