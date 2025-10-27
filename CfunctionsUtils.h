#ifndef CFUNCTIONSUTILS
#define CFUNCTIONSUTILS

#include <vector>
#include "TH1.h"
#include "Cutils.h"

double computeWeight(double pt, TH1D* hWeightTotal);
int selectPDGIndex(const particleCand& p, int PDG);
double returnPt(const std::vector<particleMC>& eventVect, float minRapidity, float maxRapidity, int PDG);

#endif
