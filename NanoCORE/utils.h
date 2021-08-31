#ifndef UTILS_H
#define UTILS_H
#include "Nano.h"
#include "Math/VectorUtil.h"

double deltaPhi(LorentzVector v1 , LorentzVector v2);
double deltaR(LorentzVector v1 , LorentzVector v2);
double deltaR_v1(LorentzVector *v1 , LorentzVector v2);
double deltaR_v2(LorentzVector *v1 , LorentzVector *v2);
void clear_branches();
vector<vector<int>> categorise( vector<int> electrons, vector<int> muons, vector<int> taus, vector<int> isoTracks );

#endif

