#ifndef UTILS_H
#define UTILS_H
#include "Nano.h"
#include "Math/VectorUtil.h"

double deltaR(LorentzVector v1 , LorentzVector v2);
double deltaR_v1(LorentzVector *v1 , LorentzVector v2);
double deltaR_v2(LorentzVector *v1 , LorentzVector *v2);

#endif

