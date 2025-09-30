#ifndef CFUNCTIONSUTILS
#define CFUNCTIONSUTILS

#include <vector>
#include "TH1.h"
#include "Cutils.h"

double computeWeight(double pt, TH1D* hWeightTotal);
double returnPt(const std::vector<particleMC>& eventVect, const float minRapidity, const float maxRapidity, const int PDG);

#endif
