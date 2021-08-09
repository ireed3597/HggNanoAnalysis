#include "utils.h"

using namespace tas;

double deltaR( LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}

double deltaR_v1( LorentzVector * v1 , LorentzVector v2){
	LorentzVector v1_1;
	v1_1.SetPt( v1->pt() );
	v1_1.SetEta( v1->eta());
	v1_1.SetPhi( v1->phi());
	v1_1.SetM( v1->mass() );
	return ROOT::Math::VectorUtil::DeltaR( v1_1 , v2);
}

double deltaR_v2( LorentzVector * v1 , LorentzVector * v2){
	LorentzVector v1_1;
	v1_1.SetPt( v1->pt() );
	v1_1.SetEta( v1->eta());
	v1_1.SetPhi( v1->phi());
	v1_1.SetM( v1->mass() );
	LorentzVector v2_1;
	v2_1.SetPt( v2->pt() );
	v2_1.SetEta( v2->eta());
	v2_1.SetPhi( v2->phi());
	v2_1.SetM( v2->mass() );
	return ROOT::Math::VectorUtil::DeltaR( v1_1 , v2_1 );
}
