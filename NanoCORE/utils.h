#ifndef UTILS_H
#define UTILS_H
#include "Nano.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "TLorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > svfit_LorentzVector;

float getCosThetaStar_CS_old( LorentzVector gg_p4, LorentzVector tt_p4, float ebeam );
float helicityCosTheta( LorentzVector booster, LorentzVector boosted);
double deltaPhi(float phi1 , float phi2 );
double deltaPhi(LorentzVector v1 , LorentzVector v2);
double deltaPhi(svfit_LorentzVector v1 , LorentzVector v2);
double deltaPhi(LorentzVector v1 , svfit_LorentzVector v2);
double deltaPhi(svfit_LorentzVector v1 , svfit_LorentzVector v2);
double deltaR(LorentzVector v1 , LorentzVector v2);
double deltaR(svfit_LorentzVector v1 , LorentzVector v2);
double deltaR(LorentzVector v1 , svfit_LorentzVector v2);
double deltaR(svfit_LorentzVector v1 , svfit_LorentzVector v2);
double deltaR(LorentzVector *v1 , LorentzVector v2);
double deltaR(LorentzVector *v1 , LorentzVector *v2);
void clear_branches();
vector<vector<int>> categorise( vector<int> electrons, vector<int> muons, vector<int> taus, vector<int> isoTracks );

#endif

