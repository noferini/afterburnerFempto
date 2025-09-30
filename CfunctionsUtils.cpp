#include "CfunctionsUtils.h"
#include "TFile.h"
#include<algorithm>
#include "Cvrun.h"

double computeWeight(double pt, TH1D* hWeightTotal){
    int bin = hWeightTotal->FindBin(pt);
    double weight = hWeightTotal->GetBinContent(bin);
    return weight;
}
//_________________________________________________________________________
double returnPt(const std::vector<particleMC>& eventVect, const float minRapidity, const float maxRapidity, const int PDG){
    double maxPt = -1;
    for(const auto& o : eventVect){
        if(o.daughters.size() > 0 || o.pdg != PDG || o.q.Pt() < maxPt){
          continue;
        }
        float y = o.q.Rapidity();
        if (y < minRapidity && y > maxRapidity)
        {
            continue;
        }
        maxPt = o.q.Pt();
    }
    return maxPt;
}