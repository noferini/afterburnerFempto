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
int selectPDGIndex(const particleCand& p, int PDG){
    for (int i = 0; i < p.pdgOptions.size(); ++i)
    {
        if (p.pdgOptions[i] == PDG)
        {
            return i;
        }
    }
    return -1;
}
//_________________________________________________________________________
double returnPt(const std::vector<particleMC>& eventVect, float minRapidity, float maxRapidity, int PDG){
    std::vector<particleCand> mcEventVect;
    for(const auto& o : eventVect){
        if(o.daughters.size() > 0){
          continue;
        }
        particleCand p;
        p.addOption(o.q.Px(), o.q.Py(), o.q.Pz(), o.pdg); // in this case only one option is set since MC truth is known

        mcEventVect.push_back(p);
    }
    std::vector<double> vPt;
    for (int i = 0; i < mcEventVect.size(); ++i)
    {
        const particleCand &p = mcEventVect[i];
        int ipdg = selectPDGIndex(p, PDG);
        if (ipdg == -1)
                continue;
        float y = p.q[ipdg].Rapidity();
        if (y >= minRapidity && y <= maxRapidity)
        {
            if (p.pdgOptions[ipdg] == PDG)
            {
                vPt.push_back(p.q[ipdg].Pt());
            }
        }
    }
    if(vPt.empty()) return std::numeric_limits<double>::quiet_NaN();
    return *std::max_element(vPt.begin(), vPt.end());
}