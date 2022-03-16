#include "utils.h"
#include "parameters.h"

using namespace tas;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

float getCosThetaStar_CS_old( LorentzVector gg_p4, svfit_LorentzVector tt_p4, float ebeam=6500 ){
    // cos theta star angle in the Collins Soper frame
    // Copied directly from here: https://github.com/ResonantHbbHgg/Selection/blob/master/selection.h#L3367-L3385
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0,  ebeam, ebeam);
    p2.SetPxPyPzE(0, 0, -ebeam, ebeam);

    LorentzVector hh_lor = gg_p4 + tt_p4;
    TLorentzVector hh;
    hh.SetPxPyPzE(hh_lor.Px(),hh_lor.Py(),hh_lor.Pz(),hh_lor.E()) ;

    TVector3 boost = - hh.BoostVector();
    p1.Boost(boost);
    p2.Boost(boost);
    LorentzVector h1_lor = gg_p4;
    TLorentzVector h_1;
    h_1.SetPxPyPzE(h1_lor.Px(),h1_lor.Py(),h1_lor.Pz(),h1_lor.E()) ; 
    h_1.Boost(boost);

    TVector3 CSaxis = p1.Vect().Unit() - p2.Vect().Unit();
    CSaxis.Unit();
    
    return TMath::Cos(   CSaxis.Angle( h_1.Vect().Unit() )    );
}

float helicityCosTheta( LorentzVector booster, LorentzVector boosted){

    TLorentzVector Booster;
    Booster.SetPxPyPzE(booster.Px(),booster.Py(),booster.Pz(),booster.E()) ; 
    TLorentzVector Boosted;
    Boosted.SetPxPyPzE(boosted.Px(),boosted.Py(),boosted.Pz(),boosted.E()) ; 

    TVector3 BoostVector = Booster.BoostVector();
    Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
    return Boosted.CosTheta();
}

float helicityCosTheta( LorentzVector booster, svfit_LorentzVector boosted){

    TLorentzVector Booster;
    Booster.SetPxPyPzE(booster.Px(),booster.Py(),booster.Pz(),booster.E()) ; 
    TLorentzVector Boosted;
    Boosted.SetPxPyPzE(boosted.Px(),boosted.Py(),boosted.Pz(),boosted.E()) ; 

    TVector3 BoostVector = Booster.BoostVector();
    Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
    return Boosted.CosTheta();
}

float helicityCosTheta( svfit_LorentzVector booster, svfit_LorentzVector boosted){

    TLorentzVector Booster;
    Booster.SetPxPyPzE(booster.Px(),booster.Py(),booster.Pz(),booster.E()) ; 
    TLorentzVector Boosted;
    Boosted.SetPxPyPzE(boosted.Px(),boosted.Py(),boosted.Pz(),boosted.E()) ; 

    TVector3 BoostVector = Booster.BoostVector();
    Boosted.Boost( -BoostVector.x(), -BoostVector.y(), -BoostVector.z() );
    return Boosted.CosTheta();
}

float helicityCosTheta_phys( LorentzVector particle1, LorentzVector particle2 ){

	TLorentzVector particle_1;
	particle_1.SetPxPyPzE(particle1.Px(),particle1.Py(),particle1.Pz(),particle1.E()) ; 
	
	TLorentzVector particle_2;
	particle_2.SetPxPyPzE(particle2.Px(),particle2.Py(),particle2.Pz(),particle2.E()) ; 

	TLorentzVector p1 = particle_1;
	TLorentzVector parent = particle_1 + particle_2;
	
	TVector3 boost_to_parent = -(parent.BoostVector());
	p1.Boost(boost_to_parent);
	
	TVector3 v1 = p1.Vect();
	TVector3 vParent = parent.Vect();
	
	double cos_theta_1 = (v1.Dot(vParent)) / (v1.Mag() * vParent.Mag());
	
	return abs(cos_theta_1);  
}

float helicityCosTheta_phys( svfit_LorentzVector particle1, svfit_LorentzVector particle2 ){

	TLorentzVector particle_1;
	particle_1.SetPxPyPzE(particle1.Px(),particle1.Py(),particle1.Pz(),particle1.E()) ; 
	
	TLorentzVector particle_2;
	particle_2.SetPxPyPzE(particle2.Px(),particle2.Py(),particle2.Pz(),particle2.E()) ; 

	TLorentzVector p1 = particle_1;
	TLorentzVector parent = particle_1 + particle_2;
	
	TVector3 boost_to_parent = -(parent.BoostVector());
	p1.Boost(boost_to_parent);
	
	TVector3 v1 = p1.Vect();
	TVector3 vParent = parent.Vect();
	
	double cos_theta_1 = (v1.Dot(vParent)) / (v1.Mag() * vParent.Mag());
	
	return abs(cos_theta_1);  
}


double deltaPhi( float phi1 , float phi2){
	float dphi = phi1 - phi2;
	if ( dphi > TMath::Pi() ) {
	    dphi -= 2.0*TMath::Pi();
	} else if ( dphi <= -TMath::Pi() ) {
	    dphi += 2.0*TMath::Pi();
	}
	return dphi;
}
double deltaPhi( LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaPhi( v1 , v2);
}
double deltaPhi( svfit_LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaPhi( v1 , v2);
}
double deltaPhi( LorentzVector v1 , svfit_LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaPhi( v1 , v2);
}
double deltaPhi( svfit_LorentzVector v1 , svfit_LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaPhi( v1 , v2);
}

double deltaR( svfit_LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}

double deltaR( svfit_LorentzVector v1 , svfit_LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}
double deltaR( LorentzVector v1 , svfit_LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}
double deltaR( LorentzVector v1 , LorentzVector v2){
	return ROOT::Math::VectorUtil::DeltaR( v1 , v2);
}

double deltaR( LorentzVector * v1 , LorentzVector v2){
	LorentzVector v1_1;
	v1_1.SetPt( v1->pt() );
	v1_1.SetEta( v1->eta());
	v1_1.SetPhi( v1->phi());
	v1_1.SetM( v1->mass() );
	return ROOT::Math::VectorUtil::DeltaR( v1_1 , v2);
}

double deltaR( LorentzVector * v1 , LorentzVector * v2){
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

vector<vector<int>> categorise( vector<int> electrons, vector<int> muons, vector<int> taus, vector<int> isoTracks ){

	vector<vector<int>> output;
	vector<int> idx_tmp(2,0);

	//2tau0lep 
	if ( taus.size() >= 2 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<taus.size(); i++){
			for (unsigned int j=i+1; j<taus.size(); j++){
				double mH = (Tau_p4().at(taus.at(i)) + Tau_p4().at(taus.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Tau_charge().at(taus.at(i)) * Tau_charge().at(taus.at(j)) < 0 ){
					idx_tmp.at(0) = taus.at(i);
					idx_tmp.at(1) = taus.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}
		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	

			output.push_back({ idx_tmp[0], 2 });
			//taus.erase( taus.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 2 });
			//taus.erase( taus.begin() + pos[1] );
	
			return output;
		}
	}

	//1tau1muon
	if ( taus.size() >= 1 && muons.size() >= 1 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<taus.size(); i++){
			for (unsigned int j=0; j<muons.size(); j++){
				double mH = (Tau_p4().at(taus.at(i)) + Muon_p4().at(muons.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Tau_charge().at(taus.at(i)) * Muon_charge().at(muons.at(j)) < 0 ){
					idx_tmp.at(0) = taus.at(i);
					idx_tmp.at(1) = muons.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}

		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	
			
			output.push_back({ idx_tmp[0], 2 });
			//taus.erase( taus.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 1 });
			//muons.erase( muons.begin() + pos[1] );
			return output;
		}
	}

	//1tau1electron
	if ( taus.size() >= 1 && electrons.size() >= 1 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<taus.size(); i++){
			for (unsigned int j=0; j<electrons.size(); j++){
				double mH = (Tau_p4().at(taus.at(i)) + Electron_p4().at(electrons.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Tau_charge().at(taus.at(i)) * Electron_charge().at(electrons.at(j)) < 0 ){
					idx_tmp.at(0) = taus.at(i);
					idx_tmp.at(1) = electrons.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}
		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	

			output.push_back({ idx_tmp[0], 2 });
			//taus.erase( taus.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 0 });
			//electrons.erase( electrons.begin() + pos[1] );
			return output;
		}
	}

	//1muon1electron
	if ( muons.size() >= 1 && electrons.size() >= 1 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<muons.size(); i++){
			for (unsigned int j=0; j<electrons.size(); j++){
				double mH = (Muon_p4().at(muons.at(i)) + Electron_p4().at(electrons.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Muon_charge().at(muons.at(i)) * Electron_charge().at(electrons.at(j)) < 0 ){
					idx_tmp.at(0) = muons.at(i);
					idx_tmp.at(1) = electrons.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}
		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	

			output.push_back({ idx_tmp[0], 1 });
			//muons.erase( muons.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 0 });
			//electrons.erase( electrons.begin() + pos[1] );
			return output;
		}
	}

	//2muons
	if ( muons.size() >= 2 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<muons.size(); i++){
			for (unsigned int j=i+1; j<muons.size(); j++){
				double mH = (Muon_p4().at(muons.at(i)) + Muon_p4().at(muons.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Muon_charge().at(muons.at(i)) * Muon_charge().at(muons.at(j)) < 0 ){
					idx_tmp.at(0) = muons.at(i);
					idx_tmp.at(1) = muons.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}

		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	
			
			output.push_back({ idx_tmp[0], 1 });
			//muons.erase( muons.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 1 });
			//muons.erase( muons.begin() + pos[1] );

			return output;
		}
	}

	//2electrons
	if ( electrons.size() >= 2 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int i=0; i<electrons.size(); i++){
			for (unsigned int j=i+1; j<electrons.size(); j++){
				double mH = (Electron_p4().at(electrons.at(i)) + Electron_p4().at(electrons.at(j))).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Electron_charge().at(electrons.at(i)) * Electron_charge().at(electrons.at(j)) < 0 ){
					idx_tmp.at(0) = electrons.at(i);
					idx_tmp.at(1) = electrons.at(j);
					pos[0] = i;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}
		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	

			output.push_back({ idx_tmp[0], 0 });
			//electrons.erase( electrons.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 0 });
			//electrons.erase( electrons.begin() + pos[1] );

			return output;
		}
	}

	//1tau1isoTrack
	if ( taus.size() >= 1 && isoTracks.size() >= 1 ){
		vector<int> pos(2,-1);
		float mH_best	= -999;
		for (unsigned int k=0; k<taus.size(); k++){
			for (unsigned int j=0; j<isoTracks.size(); j++){
				unsigned int i = isoTracks[j];
				LorentzVector iso_track(IsoTrack_pt()[i], IsoTrack_eta()[i], IsoTrack_phi()[i], 0);
				//iso_track->SetXYZT( IsoTrack_pt().at(i)* TMath::Cos(IsoTrack_phi().at(i)) , IsoTrack_pt().at(i)*TMath::Sin( IsoTrack_phi().at(i)), IsoTrack_pt().at(i)*TMath::SinH( IsoTrack_eta().at(i)),  IsoTrack_pt().at(i)*TMath::CosH( IsoTrack_eta().at(i) ) );
				double mH = (Tau_p4().at(taus.at(k)) + iso_track).M();
				if ( fabs( mH - mHiggs ) < fabs( mH_best - mHiggs ) && Tau_charge().at(taus.at(k)) * IsoTrack_pdgId().at(i) < 0 ){
					idx_tmp.at(0) = taus.at(k);
					idx_tmp.at(1) = i ;
					pos[0] = k;
					pos[1] = j;
					mH_best	= mH;
				}
			}
		}
		if ( pos[0] != -1 && pos[1] != -1 ){
			output.push_back( electrons );	
			output.push_back( muons );	
			output.push_back( taus );	
			output.push_back( isoTracks );	

			output.push_back({ idx_tmp[0], 2 });
			//taus.erase( taus.begin() + pos[0] );
			output.push_back({ idx_tmp[1], 3 });
			//isoTracks.erase( isoTracks.begin() + pos[1] );
			return output;
		}
	}

	//1tau-noIsoTrack
	if ( taus.size() >= 1 ){
		output.push_back( electrons );	
		output.push_back( muons );	
		output.push_back( taus );	
		output.push_back( isoTracks );	

		output.push_back({ taus[0], 2 });
		taus.erase( taus.begin() + 0 );
		output.push_back({ -1 , -1 });
		return output;
	}

	else{
		output.push_back( electrons );	
		output.push_back( muons );	
		output.push_back( taus );	
		output.push_back( isoTracks );	

		output.push_back({ -1 , -1 });
		output.push_back({ -1 , -1 });
		return output;
	}
}
	
void clear_branches(){

	t_run					= -9;
	t_lumiBlock				= -9;
	t_event					= -9;
	t_MET_pt				= -9;
	t_MET_phi				= -9;
	t_weight				= -9;
	//category				= -9; 	//initialised in looper

	lep12_dphi				= -9;
	lep12_deta				= -9;
	lep12_deta_bdt			= -9;
	lep12_dr				= -9;

	cat1					= false;
	cat2					= false;
	cat3					= false;
	cat4					= false;
	cat5					= false;
	cat6					= false;
	cat7					= false;
	cat8					= false;

	n_electrons				= -9;
	n_muons					= -9;
	n_taus 					= -9;
	n_jets 					= -9;
	n_bjets					= -9;
	n_isoTrks				= -9;

	g1_ptmgg				= -9;
	g1_pt					= -9;
	g1_eta					= -9;
	g1_eta_bdt				= -9;
	g1_phi					= -9;
	g1_idmva				= -9;
	g1_pixVeto				= false;

	g2_ptmgg				= -9;
	g2_pt					= -9;
	g2_eta					= -9;
	g2_eta_bdt				= -9;
	g2_phi					= -9;
	g2_idmva				= -9;
	g2_pixVeto				= false;

	gg_pt					= -9;
	gg_ptmgg				= -9;
	gg_eta					= -9;
	gg_eta_bdt				= -9;
	gg_phi					= -9;
	gg_dR					= -9;
	gg_dPhi					= -9;
	gg_hel					= -9;
	gg_hel_phys				= -9;
	gg_tt_CS				= -9;
	gg_tt_hel				= -9;
	gg_tt_hel_phys			= -9;
	mgg						= -9;

	lep1_pt					= -9;
	lep1_eta				= -9;
	lep1_eta_bdt			= -9;
	lep1_phi				= -9;
	lep1_charge				= -9;
	lep1_pdgID				= -9;
	lep1_tightID			= -9;
	lep1_id_vs_e			= -9;
	lep1_id_vs_m			= -9;
	lep1_id_vs_jet			= -9;

	lep2_pt        			= -9;
	lep2_eta       			= -9;
	lep2_eta_bdt			= -9;
	lep2_phi       			= -9;
	lep2_charge    			= -9;
	lep2_pdgID     			= -9;
	lep2_tightID   			= -9;
	lep2_id_vs_e   			= -9;
	lep2_id_vs_m   			= -9;
	lep2_id_vs_jet 			= -9;

	jet1_pt					= -9;
	jet1_eta				= -9;
	jet1_eta_bdt			= -9;
	jet1_phi				= -9;
	jet1_bTag				= -9;
	jet1_id					= -9;

	jet2_pt					= -9;
	jet2_eta				= -9;
	jet2_eta_bdt			= -9;
	jet2_phi				= -9;
	jet2_bTag				= -9;
	jet2_id					= -9;

	max_bTag				= -9;
	tt_hel					= -9;
	tt_hel_phys				= -9;

	m_tautau_vis			= -9;
	pt_tautau_vis			= -9;
	eta_tautau_vis			= -9;
	eta_tautau_vis_bdt		= -9;
	phi_tautau_vis			= -9;

	MET_gg_dPhi				= -9;
	MET_ll_dPhi				= -9;
	dPhi_MET_l				= -9;
	ll_dPhi					= -9;
	ll_dEta					= -9;
	ll_dR					= -9;
	m_Z						= -9;

	dZ						= 0;
	g1_energyErr	= -9;
	g2_energyErr	= -9;
	max_g_ptmgg		= -9;
	min_g_ptmgg		= -9;
	max_g_idmva		= -9;
	min_g_idmva		= -9;
	lep2_pfRelIso03_all	= -9;
	lep2_pfRelIso03_chg	= -9;
	max_lep_pt		= -9;
	min_lep_pt		= -9;

	pt_tautauSVFitLoose		= -9;
	eta_tautauSVFitLoose	= -9;
	eta_tautauSVFitLoose_bdt= -9;
	phi_tautauSVFitLoose	= -9;
	m_tautauSVFitLoose		= -9;
	dR_tautauSVFitLoose		= -9;
	dR_ggtautauSVFitLoose	= -9;
	dPhi_tautauSVFitLoose	= -9;
	dPhi_ggtautauSVFitLoose	= -9;

	tau1_pt_SVFit			= -9;
	tau1_eta_SVFit		= -9;
	tau1_phi_SVFit		= -9;
	tau1_m_SVFit			= -9;
	tau2_pt_SVFit			= -9;
	tau2_eta_SVFit		= -9;
	tau2_phi_SVFit		= -9;
	tau2_m_SVFit			= -9;

	pt_tautau_sntMtt		= -9;
	eta_tautau_sntMtt	= -9;
	eta_tautau_sntMtt_bdt= -9;
	phi_tautau_sntMtt	= -9;
	m_tautau_sntMtt		= -9;
	dR_tautau_sntMtt		= -9;
	dR_ggtautau_sntMtt	= -9;
	dPhi_tautau_sntMtt	= -9;
	dPhi_ggtautau_sntMtt	= -9;

	tau1_pt_sntMtt			= -9;
	tau1_eta_sntMtt		= -9;
	tau1_phi_sntMtt		= -9;
	tau1_m_sntMtt			= -9;
	tau2_pt_sntMtt			= -9;
	tau2_eta_sntMtt		= -9;
	tau2_phi_sntMtt		= -9;
	tau2_m_sntMtt			= -9;

	mX						= -9;
	m_llg_lead			= -9;
	m_llg_subl			= -9;
}
