#include "utils.h"
#include "parameters.h"

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

vector<vector<int>> categorise( vector<int> electrons, vector<int> muons, vector<int> taus, vector<int> isoTracks ){

	vector<vector<int>> output;
	vector<int> idx_tmp(2,0);

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

	cat1					= false;
	cat2					= false;
	cat3					= false;
	cat4					= false;
	cat5					= false;
	cat6					= false;
	cat7					= false;
	cat8					= false;

	mgg						= -9;
	n_electrons				= -9;
	n_muons					= -9;
	n_taus 					= -9;
	n_jets 					= -9;
	n_isoTrks				= -9;

	g1_ptmgg				= -9;
	g1_pt					= -9;
	g1_eta					= -9;
	g1_phi					= -9;
	g1_idmva				= -9;
	g1_pixVeto				= -9;

	g2_ptmgg				= -9;
	g2_pt					= -9;
	g2_eta					= -9;
	g2_phi					= -9;
	g2_idmva				= -9;
	g2_pixVeto				= -9;

	gg_pt					= -9;
	gg_eta					= -9;
	gg_phi					= -9;
	gg_dR					= -9;

	lep1_pt					= -9;
	lep1_eta				= -9;
	lep1_phi				= -9;
	lep1_charge				= -9;
	lep1_pdgID				= -9;
	lep1_tightID			= -9;
	lep1_id_vs_e			= -9;
	lep1_id_vs_m			= -9;
	lep1_id_vs_jet			= -9;

	lep2_pt        			= -9;
	lep2_eta       			= -9;
	lep2_phi       			= -9;
	lep2_charge    			= -9;
	lep2_pdgID     			= -9;
	lep2_tightID   			= -9;
	lep2_id_vs_e   			= -9;
	lep2_id_vs_m   			= -9;
	lep2_id_vs_jet 			= -9;

	jet1_pt					= -9;
	jet1_eta				= -9;
	jet1_bTag				= -9;
	jet1_id					= -9;

	jet2_pt					= -9;
	jet2_eta				= -9;
	jet2_bTag				= -9;
	jet2_id					= -9;

	pt_tautauSVFitLoose		= -9;
	eta_tautauSVFitLoose	= -9;
	phi_tautauSVFitLoose	= -9;
	m_tautauSVFitLoose		= -9;
	dR_tautauSVFitLoose		= -9;
	dR_ggtautauSVFitLoose	= -9;
	dPhi_MET_tau1			= -9;

	m_tautau_vis			= -9;
	pt_tautau_vis			= -9;
	eta_tautau_vis			= -9;
	phi_tautau_vis			= -9;
}
