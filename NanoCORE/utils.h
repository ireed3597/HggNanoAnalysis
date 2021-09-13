#ifndef UTILS_H
#define UTILS_H
#include "Nano.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include "TLorentzVector.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > svfit_LorentzVector;

float getCosThetaStar_CS_old( LorentzVector gg_p4, svfit_LorentzVector tt_p4, float ebeam );
float helicityCosTheta( LorentzVector booster, LorentzVector boosted);
float helicityCosTheta( svfit_LorentzVector booster, svfit_LorentzVector boosted);
float helicityCosTheta_phys( LorentzVector particle_1, LorentzVector particle_2);
float helicityCosTheta_phys( svfit_LorentzVector particle_1, svfit_LorentzVector particle_2);
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

