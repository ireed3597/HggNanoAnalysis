#pragma GCC diagnostic ignored "-Wsign-compare"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TTree.h"
#include "TChain.h"
#include "THStack.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TObjString.h"
#include "TTreeCache.h"
#include "TTreeCacheUnzip.h"
#include "TTreePerfStats.h"

#include "../NanoCORE/Nano.cc"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/utils.cc"
#include "SVfit_utils.cc"

#include <string>
#include <iostream>
#include <iomanip>

#define SUM(vec) std::accumulate((vec).begin(), (vec).end(), 0);
#define SUM_GT(vec,num) std::accumulate((vec).begin(), (vec).end(), 0, [](float x,float y){return ((y > (num)) ? x+y : x); });
#define COUNT_GT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x > (num); });
#define COUNT_LT(vec,num) std::count_if((vec).begin(), (vec).end(), [](float x) { return x < (num); });

#define H1(name,nbins,low,high) TH1F *h_##name = new TH1F(#name,#name,nbins,low,high);

// #define DEBUG

struct debugger { template<typename T> debugger& operator , (const T& v) { cerr<<v<<" "; return *this; } } dbg;
#ifdef DEBUG
    #define debug(args...) do {cerr << #args << ": "; dbg,args; cerr << endl;} while(0)
#else
    #define debug(args...)
#endif

using namespace std;
using namespace tas;

int ScanChain( TChain *ch, string proc, string str_year, float scale_factor = 1, bool resonant = false ) {

	int year;
	if ( str_year == "2016_APV") year = 2016;
	else { year = stoi(str_year); }

	TString file_name = proc + "_ReRecoCheck_" +  str_year;

	TFile* f1 = new TFile("outputs/" + file_name + ".root", "RECREATE");
	H1(mgg, 60, 100 , 180 );
	H1(mgg_1t0l, 60, 100 , 180 );
	H1(mgg_1t1l, 60, 100 , 180 );
	H1(mgg_2t0l, 60, 100 , 180 );
	H1(mgg_0t2l, 60, 100 , 180 );
	H1(mgg_1t0l_iso, 60, 100 , 180 );

	bool ggf_samples = false;
	//define process ids
	if (proc.find(std::string("HH_ggWW_semileptonic")) != std::string::npos)	process_id = -4;
	else if (proc.find(std::string("HH_ggWW_dileptonic")) != std::string::npos) 	process_id = -3;
	//else if (proc.find(std::string("HH_ggZZ")) != std::string::npos) 				process_id = -2;
	else if (proc.find(std::string("HH_ggZZ_2l2q")) != std::string::npos) 				process_id = -6;
	else if (proc.find(std::string("HH_ggZZ_4l")) != std::string::npos) 				process_id = -5;
	else if (proc.find(std::string("HH_ggTauTau")) != std::string::npos) 			process_id = -1;
	else if (proc.find(std::string("Data")) != std::string::npos) 				process_id = 0;
	else if (proc.find(std::string("ZGamma")) != std::string::npos) 				process_id = 2;
	else if (proc.find(std::string("DiPhoton")) != std::string::npos) 			process_id = 3;
	else if (proc.find(std::string("WGamma")) != std::string::npos) 				process_id = 4;
	else if (proc.find(std::string("TTbar")) != std::string::npos) 				process_id = 5;
	else if (proc.find(std::string("TTGamma")) != std::string::npos) 				process_id = 6;
	else if (proc.find(std::string("TTGG")) != std::string::npos) 				process_id = 7;
	else if (proc.find(std::string("GJets")) != std::string::npos)				process_id = 8;
	else if (proc.find(std::string("VH")) != std::string::npos) 					process_id = 9;
	else if (proc.find(std::string("ttH")) != std::string::npos) 					process_id = 10;
	else if (proc.find(std::string("ggH")) != std::string::npos) 					process_id = 11;
	else if (proc.find(std::string("VBFH")) != std::string::npos) 				process_id = 12;
	else if (proc.find(std::string("WW")) != std::string::npos) 				process_id = 14;
	else if (proc.find(std::string("WZ")) != std::string::npos) 				process_id = 15;
	else if (proc.find(std::string("ZZ")) != std::string::npos) 				process_id = 16;
	//non-SM couplings
	else if (proc.find(std::string("ggf_c0_HHggtautau")) != std::string::npos) 		{		process_id = -10;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c1_HHggtautau")) != std::string::npos) 		{		process_id = -11;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c2p45_HHggtautau")) != std::string::npos) {				process_id = -12;		ggf_samples = true;		}
	else if (proc.find(std::string("ggf_c5_HHggtautau")) != std::string::npos) 		{		process_id = -13;		ggf_samples = true;		}
	else if (proc.find(std::string("vbf_cv_0p5_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -14;
	else if (proc.find(std::string("vbf_cv_1p5_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -15;
	else if (proc.find(std::string("vbf_cv_1_c2v_0_c3_1_HHggtautau")) != std::string::npos) 				process_id = -16;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_0_HHggtautau")) != std::string::npos) 				process_id = -17;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_1_HHggtautau")) != std::string::npos) 				process_id = -18;
	else if (proc.find(std::string("vbf_cv_1_c2v_1_c3_2_HHggtautau")) != std::string::npos) 				process_id = -19;
	else if (proc.find(std::string("vbf_cv_1_c2v_2_c3_1_HHggtautau")) != std::string::npos) 				process_id = -20;

    int nEventsTotal = 0;
    int nEventsChain = ch->GetEntries();
    TFile *currentFile = 0;
    TObjArray *listOfFiles = ch->GetListOfFiles();
    TIter fileIter(listOfFiles);
    tqdm bar;

    while ( (currentFile = (TFile*)fileIter.Next()) ) {
        TFile *file = TFile::Open( currentFile->GetTitle() );
        TTree *tree = (TTree*)file->Get("Events");

        tree->SetCacheSize(128*1024*1024);
        tree->SetCacheLearnEntries(100);

        //auto psRead = new TTreePerfStats("readPerf", tree);
				nt.SetYear(year);
        nt.Init(tree);

        for( unsigned int loop_event = 0; loop_event < tree->GetEntriesFast(); ++loop_event) {

         nt.GetEntry(loop_event);

         nEventsTotal++;
         bar.progress(nEventsTotal, nEventsChain);

				clear_branches();

				float weight = 1.;
				if ( proc != "Data" ) weight = genWeight() * scale_factor;

				//////////////////////////////////////////////////////////////////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////////////////////
				///////////////////							photon selection
				vector<int> pho_cands;
				vector<float> pho_pt_cands;
				for (unsigned int i=0; i<nPhoton(); i++){
					if ( Photon_electronVeto().at(i)  >=0.5 &&
						 (	(	Photon_isScEtaEB().at(i)	 && Photon_r9().at(i) > 0.85 )			//pho_EB_highR9
						||	(	Photon_isScEtaEE().at(i)	 && Photon_r9().at(i) > 0.90 )			//pho_EE_highR9
						||	(	Photon_isScEtaEB().at(i)	 && Photon_r9().at(i) < 0.85 && Photon_r9().at(i) > 0.5 && Photon_sieie().at(i) < 0.015 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_photonIso().at(i) - 0.16544*fixedGridRhoAll() ) < 4.0 )			//pho_EB_lowR9
						||	(	Photon_isScEtaEE().at(i)	 && Photon_r9().at(i) < 0.90 && Photon_r9().at(i) > 0.8 && Photon_sieie().at(i) < 0.035 && Photon_trkSumPtHollowConeDR03().at(i) < 6.0  && ( Photon_photonIso().at(i) - 0.13212*fixedGridRhoAll() ) < 4.0 )			/*pho_EE_lowR9 */ )
					 	&& 	Photon_hoe().at(i) < 0.08
						&&	fabs(Photon_eta().at(i)) < 2.5
						&&	(	fabs(Photon_eta().at(i)) < 1.442 || fabs(Photon_eta().at(i)) > 1.566 )
						&&	( Photon_r9().at(i) > 0.8 || Photon_chargedHadronIso().at(i) < 20 || Photon_chargedHadronIso().at(i) / Photon_pt().at(i) < 0.3 )
						&&  ( Photon_isScEtaEB().at(i) || Photon_isScEtaEE().at(i) )
						&&    Photon_mvaID().at(i) > pho_idmva_cut
						){
							pho_cands.push_back(i);					
							pho_pt_cands.push_back( Photon_pt().at(i) );					
					}
				}
				if ( pho_cands.size() < 2 ) continue;

				sort(pho_pt_cands.begin(), pho_pt_cands.end(), greater<float>());
				int gHidx[2] = {-1,-1};
				for (unsigned int i=0; i<pho_cands.size();i++){
					if ( Photon_pt().at(pho_cands[i]) == pho_pt_cands[0] && Photon_pt().at(pho_cands[i]) > 35 ) gHidx[0]	= pho_cands[i];
					if ( Photon_pt().at(pho_cands[i]) == pho_pt_cands[1] && Photon_pt().at(pho_cands[i]) > 25 ) gHidx[1]	= pho_cands[i];
				}
				if ( gHidx[0] < 0 || gHidx[1] < 0 ) continue;

			//di-photon selection
			mgg = (float)(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]) ).M();
			if ( mgg < mgg_lower || mgg > mgg_upper ) continue;
			if ( (proc == "Data" || !resonant) && mgg > mgg_sideband_lower && mgg < mgg_sideband_upper ) continue;

			if ( ( Photon_pt().at(gHidx[0]) < pho_pt_cut 				|| Photon_pt().at(gHidx[1]) < pho_pt_cut ) 
				|| ( Photon_pt().at(gHidx[0]) / mgg < lead_pt_mgg_cut 	|| Photon_pt().at(gHidx[1]) / mgg < sublead_pt_mgg_cut ) 
				|| ( Photon_mvaID().at(gHidx[0]) < pho_idmva_cut 			|| Photon_mvaID().at(gHidx[1]) < pho_idmva_cut ) 
				|| ( Photon_electronVeto().at(gHidx[0]) < pho_eveto_cut 	|| Photon_electronVeto().at(gHidx[1]) < pho_eveto_cut ) 
				|| ( fabs(Photon_eta().at(gHidx[0])) > pho_eta_cut 		|| fabs(Photon_eta().at(gHidx[1])) > pho_eta_cut ) 
				|| ( fabs(Photon_eta().at(gHidx[0])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx[0])) < trans_eta_high ) 
				|| ( fabs(Photon_eta().at(gHidx[1])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx[1])) < trans_eta_high ) 
				){ continue; }
				

			vector<int> sel_eles;
			for(unsigned int i=0; i<nElectron(); i++){
				if (Electron_pt().at(i) > ele_pt && fabs(Electron_eta().at(i)) < ele_eta && ( fabs(Electron_eta().at(i)) < trans_eta_low || fabs(Electron_eta().at(i)) > trans_eta_high ) && fabs(Electron_dxy().at(i)) < ele_dxy && fabs(Electron_dz().at(i)) < ele_dz 
				&& ( (Electron_pfRelIso03_all().at(i) < ele_pfRelIso && Electron_mvaFall17V2noIso_WP90().at(i)) || (Electron_mvaFall17V2Iso_WP90().at(i) ) ) 
				&& deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx[0]) ) > ele_dR_pho && deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx[1]) ) > ele_dR_pho 
				){
					sel_eles.push_back(i);
				}
			}
			
			vector<int> sel_muons;
			for(unsigned int i=0; i<nMuon(); i++){
				if (Muon_pt().at(i) > muon_pt && fabs(Muon_eta().at(i)) < muon_eta && fabs(Muon_dxy().at(i)) < muon_dxy && fabs(Muon_dz().at(i)) < muon_dz 
				&& Muon_pfRelIso03_all().at(i) < muon_pfRelIso 
				&& Muon_isGlobal().at(i) 
				&& deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx[0]) ) > muon_dR_pho && deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx[1]) ) > muon_dR_pho 
				){
					sel_muons.push_back(i);
				}
			}
			
			vector<int> sel_taus;
			for(unsigned int i=0; i<nTau(); i++){
				if (Tau_pt().at(i) > tau_pt && fabs(Tau_eta().at(i)) < tau_eta && /*Tau_idDecayModeNewDMs().at(i) &&*/ fabs(Tau_dz().at(i)) < tau_dz 
				&& Tau_idDeepTau2017v2p1VSe().at(i) >= tau_deepID_e && Tau_idDeepTau2017v2p1VSmu().at(i) >= tau_deepID_m && Tau_idDeepTau2017v2p1VSjet().at(i) >= tau_deepID_j 
				&& deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx[0]) ) > tau_dR_pho && deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx[1]) ) > tau_dR_pho 
				){

					bool overlap = false;
					for (unsigned int j=0; j<sel_eles.size(); j++){
						if ( deltaR( Tau_p4().at(i) , Electron_p4().at(sel_eles.at(j)) ) < tau_dR_lep ) overlap = true;
					}
					for (unsigned int j=0; j<sel_muons.size(); j++){
						if ( deltaR( Tau_p4().at(i) , Muon_p4().at(sel_muons.at(j)) ) < tau_dR_lep ) overlap = true;
					}

					if ( !overlap ) sel_taus.push_back(i);
				}
			}

			//Isolated Tracks selection
            vector<int> sel_isoTracks;
            for (unsigned int i=0; i<nIsoTrack(); i++){
                if ( IsoTrack_isPFcand().at(i) && IsoTrack_fromPV().at(i) 
									&&	fabs(IsoTrack_dxy().at(i)) < 0.2
									&&	fabs(IsoTrack_dz().at(i)) < 0.1
									){
                    LorentzVector *iso_track = new LorentzVector;
                    iso_track->SetXYZT( IsoTrack_pt().at(i)* TMath::Cos(IsoTrack_phi().at(i)) , IsoTrack_pt().at(i)*TMath::Sin( IsoTrack_phi().at(i)), IsoTrack_pt().at(i)*TMath::SinH( IsoTrack_eta().at(i)),  IsoTrack_pt().at(i)*TMath::CosH( IsoTrack_eta().at(i) ) );
                    if ( deltaR( iso_track , Photon_p4().at(gHidx[0]) ) > isoTrk_dR  && deltaR( iso_track , Photon_p4().at(gHidx[1]) ) > isoTrk_dR  ){
						bool iso = true;
						for (unsigned int j=0; j<sel_eles.size(); j++){
                    		if ( deltaR( iso_track , Electron_p4().at(sel_eles.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						for (unsigned int j=0; j<sel_muons.size(); j++){
                    		if ( deltaR( iso_track , Muon_p4().at(sel_muons.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						for (unsigned int j=0; j<sel_taus.size(); j++){
                    		if ( deltaR( iso_track , Tau_p4().at(sel_taus.at(j)) ) < isoTrk_dR ){
								iso = false;
                    		}
                		}
						if ( iso ) sel_isoTracks.push_back( i );
					}
				}
			}
			//Jet Selection			
			vector<int> sel_jets;
			for(unsigned int i=0; i<nJet(); i++){
				if (Jet_pt().at(i) > jet_pt && fabs(Jet_eta().at(i)) < jet_eta && Jet_neEmEF().at(i) < jet_neEmEF && Jet_neHEF().at(i) < jet_neHEF && Jet_chHEF()[i] > jet_chHEF && Jet_chEmEF()[i] < jet_chEmEF && Jet_nConstituents()[i] > jet_nConstituents && deltaR( Jet_p4().at(i) , Photon_p4().at(gHidx[0]) ) > jet_dR_pho && deltaR( Jet_p4().at(i) , Photon_p4().at(gHidx[1]) ) > jet_dR_pho ){

					bool overlap = false;
					for (unsigned int j=0; j<sel_eles.size(); j++){
						if ( !overlap && deltaR( Jet_p4().at(i) , Electron_p4().at(sel_eles.at(j)) ) < jet_dR_lep ){
							overlap = true;
							break;
						}
					}
					for (unsigned int j=0; j<sel_muons.size(); j++){
						if ( !overlap && deltaR( Jet_p4().at(i) , Muon_p4().at(sel_muons.at(j)) ) < jet_dR_lep ){
							overlap = true;
							break;
						}
					}
					for (unsigned int j=0; j<sel_taus.size(); j++){
						if ( !overlap && deltaR( Jet_p4().at(i) , Tau_p4().at(sel_taus.at(j)) ) < jet_dR_tau ){
							overlap = true;
							break;
						}
					}

					if ( !overlap ){
						sel_jets.push_back(i);
						if ( Jet_btagDeepFlavB().at(i) > max_bTag ) max_bTag = Jet_btagDeepFlavB().at(i);
					}
				}
			}

			//bJet Selection			
			vector<int> sel_bJets;
			for(unsigned int i=0; i<sel_jets.size(); i++){
				if ( Jet_btagDeepFlavB().at(sel_jets[i]) > bTag_medium_WP[year - 2016] ) sel_bJets.push_back(sel_jets[i]);
			}


			//Z veto cut
			bool Z_cand = false;
			if ( sel_eles.size() >= 2 ){
				for (unsigned int i=0; i<sel_eles.size(); i++){
					for (unsigned int j=i+1; j<sel_eles.size(); j++){
						if ( (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() > mZ_veto_low  && (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() < mZ_veto_up 
							&& Electron_charge().at(sel_eles[i]) * Electron_charge().at(sel_eles[j]) < 0 ){
							Z_cand = true;
							break;
						}
					}
				}
			}
			if ( sel_muons.size() >= 2 ){
				for (unsigned int i=0; i<sel_muons.size(); i++){
					for (unsigned int j=i+1; j<sel_muons.size(); j++){
						if ( (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() > mZ_veto_low  && (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() < mZ_veto_up 
							&& Muon_charge().at(sel_muons[i]) * Muon_charge().at(sel_muons[j]) < 0 ){
							Z_cand = true;
							break;
						}
					}
				}
			}
			if ( Z_cand ) continue;

			n_electrons		= sel_eles.size();
			n_muons			= sel_muons.size();
			n_taus			= sel_taus.size();
			n_isoTrks		= sel_isoTracks.size();
			n_jets			= sel_jets.size();
			n_bjets			= sel_bJets.size();

			vector<int> h_cand1, h_cand2;
			vector<vector<int>> raw_results = categorise( sel_eles, sel_muons, sel_taus, sel_isoTracks );
			sel_eles		= raw_results[0];
			sel_muons		= raw_results[1];
			sel_taus		= raw_results[2];
			sel_isoTracks	= raw_results[3];
			h_cand1			= raw_results[4];
			h_cand2			= raw_results[5];

			if ( h_cand1[1] == -1 && h_cand2[1] == -1  ) continue;

			if ( h_cand1[1] == 2 && h_cand2[1] == 1  ) cat1 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 0  ) cat2 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 2  ) cat3 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 1  ) cat4 = true;
			if ( h_cand1[1] == 0 && h_cand2[1] == 0  ) cat5 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 0  ) cat6 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 3  ) cat7 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == -1 ) cat8 = true;

			category = 99;
			if ( cat1 ) category = 1;
			if ( cat2 ) category = 2;
			if ( cat3 ) category = 3;
			if ( cat4 ) category = 4;
			if ( cat5 ) category = 5;
			if ( cat6 ) category = 6;
			if ( cat7 ) category = 7;
			if ( cat8 ) category = 8;


			t_run			= run();
			t_lumiBlock		= luminosityBlock();
			t_event			= event();
			t_MET_pt		= MET_pt();
			t_MET_phi		= MET_phi();
			t_weight		= weight;
			MET_gg_dPhi		= deltaPhi( MET_phi() , (Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).phi() );

			gg_pt			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).pt() ;
			gg_ptmgg		=	gg_pt / mgg;
			gg_eta			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).eta() ;
			//gg_eta_bdt		=	gg_eta * sgn( gg_eta ) ;
			gg_phi			=	(Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1])).phi() ;
			gg_dR			=	deltaR(Photon_p4().at(gHidx[0]) , Photon_p4().at(gHidx[1])) ;
			gg_dPhi			=	deltaPhi(Photon_p4().at(gHidx[0]) , Photon_p4().at(gHidx[1])) ;
			gg_hel_phys		= 	fabs(helicityCosTheta_phys( Photon_p4().at(gHidx[0]), Photon_p4().at(gHidx[1]) ) ); 
			bool roll		= 	rand() % 2 == 0;
			if (roll ) 		gg_hel			= 	fabs( helicityCosTheta( Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]), Photon_p4().at(gHidx[0]) ) ) ;
			else	{		gg_hel			= 	fabs( helicityCosTheta( Photon_p4().at(gHidx[0]) + Photon_p4().at(gHidx[1]), Photon_p4().at(gHidx[1]) ) ) ; }

			g1_ptmgg		=	Photon_pt().at(gHidx[0]) / mgg;
			g1_pt			=	Photon_pt().at(gHidx[0]) ;
			g1_eta			=	Photon_eta().at(gHidx[0]) ;
			g1_eta_bdt		=	g1_eta * sgn( gg_eta ) ;
			g1_phi			=	Photon_phi().at(gHidx[0]) ;
			g1_idmva		=	Photon_mvaID().at(gHidx[0]) ;
			g1_pixVeto		=   Photon_pixelSeed().at(gHidx[0]) ;
			g1_energyErr	=   Photon_energyErr().at(gHidx[0]) ;

			g2_ptmgg		=	Photon_pt().at(gHidx[1]) / mgg;
			g2_pt			=	Photon_pt().at(gHidx[1]) ;
			g2_eta			=	Photon_eta().at(gHidx[1]) ;
			g2_eta_bdt		=	g2_eta * sgn( gg_eta ) ;
			g2_phi			=	Photon_phi().at(gHidx[1]) ;
			g2_idmva		=	Photon_mvaID().at(gHidx[1]) ;
			g2_pixVeto		=   Photon_pixelSeed().at(gHidx[1]) ;
			g2_energyErr	=   Photon_energyErr().at(gHidx[1]) ;

			if ( g1_pt > g2_pt ){
				max_g_ptmgg = g1_pt / mgg;
				min_g_ptmgg = g2_pt / mgg;
			}
			else{
				max_g_ptmgg = g2_pt / mgg;
				min_g_ptmgg = g1_pt / mgg;
			}
			if ( g1_idmva > g2_idmva ){
				max_g_idmva = g1_idmva;
				min_g_idmva = g2_idmva;
			}
			else{
				max_g_idmva = g2_idmva;
				min_g_idmva = g1_idmva;
			}

			LorentzVector lep1_p4, lep2_p4;

			if ( cat1 ){
				lep1_p4			=	Muon_p4()[h_cand2[0]];
				lep1_pt			=	Muon_pt()[h_cand2[0]];
				lep1_eta		=	Muon_eta()[h_cand2[0]];
				lep1_eta_bdt	=	Muon_eta()[h_cand2[0]] * sgn( gg_eta );
				lep1_phi		=	Muon_phi()[h_cand2[0]];	
				lep1_charge		=	Muon_charge()[h_cand2[0]];
				lep1_pdgID		=	Muon_pdgId()[h_cand2[0]];
				lep1_tightID	=	Muon_tightId()[h_cand2[0]];

				lep2_p4			=	Tau_p4()[h_cand1[0]];
				lep2_pt			=	Tau_pt()[h_cand1[0]];
				lep2_eta		=	Tau_eta()[h_cand1[0]];
				lep2_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
				lep2_phi		=	Tau_phi()[h_cand1[0]];	
				lep2_charge		=	Tau_charge()[h_cand1[0]];
				lep2_pdgID		=	15 * sgn( lep2_charge );
				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
			
				lep12_dr		= deltaR( Muon_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
			}
			if ( cat2 ){
				lep1_p4			=	Electron_p4()[h_cand2[0]];
				lep1_pt			=	Electron_pt()[h_cand2[0]];
				lep1_eta		=	Electron_eta()[h_cand2[0]];
				lep1_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
				lep1_phi		=	Electron_phi()[h_cand2[0]];	
				lep1_charge		=	Electron_charge()[h_cand2[0]];
				lep1_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
				lep1_pdgID		=	Electron_pdgId()[h_cand2[0]];

				lep2_p4			=	Tau_p4()[h_cand1[0]];
				lep2_pt			=	Tau_pt()[h_cand1[0]];
				lep2_eta		=	Tau_eta()[h_cand1[0]];
				lep2_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
				lep2_phi		=	Tau_phi()[h_cand1[0]];	
				lep2_charge		=	Tau_charge()[h_cand1[0]];
				lep2_pdgID		=	15 * sgn( lep2_charge );
				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];

				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
			}
			if ( cat3 ){
				lep1_p4			=	Tau_p4()[h_cand1[0]];
				lep1_pt			=	Tau_pt()[h_cand1[0]];
				lep1_eta		=	Tau_eta()[h_cand1[0]];
				lep1_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi		=	Tau_phi()[h_cand1[0]];	
				lep1_charge		=	Tau_charge()[h_cand1[0]];
				lep1_pdgID		=	15 * sgn( lep1_charge );
				lep1_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
				lep1_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
				lep1_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];

				lep2_p4			=	Tau_p4()[h_cand2[0]];
				lep2_pt			=	Tau_pt()[h_cand2[0]];
				lep2_eta		=	Tau_eta()[h_cand2[0]];
				lep2_eta_bdt	=	Tau_eta()[h_cand2[0]] * sgn( gg_eta );
				lep2_phi		=	Tau_phi()[h_cand2[0]];	
				lep2_charge		=	Tau_charge()[h_cand2[0]];
				lep2_pdgID		=	15 * sgn( lep2_charge );
				lep2_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand2[0]];	
				lep2_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand2[0]];	
				lep2_id_vs_jet	=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand2[0]];

				lep12_dr		= deltaR( Tau_p4()[h_cand2[0]], Tau_p4()[h_cand1[0]] );
			}
			if ( cat4 ){
				lep1_p4			=	Muon_p4()[h_cand1[0]];
				lep1_pt			=	Muon_pt()[h_cand1[0]];
				lep1_eta		=	Muon_eta()[h_cand1[0]];
				lep1_eta_bdt	=	Muon_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi		=	Muon_phi()[h_cand1[0]];	
				lep1_charge		=	Muon_charge()[h_cand1[0]];
				lep1_pdgID		=	Muon_pdgId()[h_cand1[0]];
				lep1_tightID	=	Muon_tightId()[h_cand1[0]];

				lep2_p4			=	Muon_p4()[h_cand2[0]];
				lep2_pt			=	Muon_pt()[h_cand2[0]];
				lep2_eta		=	Muon_eta()[h_cand2[0]];
				lep2_eta_bdt	=	Muon_eta()[h_cand2[0]] * sgn( gg_eta );
				lep2_phi		=	Muon_phi()[h_cand2[0]];	
				lep2_charge		=	Muon_charge()[h_cand2[0]];
				lep2_pdgID		=	Muon_pdgId()[h_cand2[0]];
				lep2_tightID	=	Muon_tightId()[h_cand2[0]];

				lep12_dr		= deltaR( Muon_p4()[h_cand2[0]], Muon_p4()[h_cand1[0]] );
				m_Z				= ( lep1_p4 + lep2_p4 ).M();
			}
			if ( cat5 ){
				lep1_p4			=	Electron_p4()[h_cand1[0]];
				lep1_pt			=	Electron_pt()[h_cand1[0]];
				lep1_eta		=	Electron_eta()[h_cand1[0]];
				lep1_eta_bdt	=	Electron_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi		=	Electron_phi()[h_cand1[0]];	
				lep1_charge		=	Electron_charge()[h_cand1[0]];
				lep1_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand1[0]];
				lep1_pdgID		=	Electron_pdgId()[h_cand1[0]];

				lep2_p4			=	Electron_p4()[h_cand2[0]];
				lep2_pt			=	Electron_pt()[h_cand2[0]];
				lep2_eta		=	Electron_eta()[h_cand2[0]];
				lep2_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
				lep2_phi		=	Electron_phi()[h_cand2[0]];	
				lep2_charge		=	Electron_charge()[h_cand2[0]];
				lep2_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
				lep2_pdgID		=	Electron_pdgId()[h_cand2[0]];

				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Electron_p4()[h_cand1[0]] );
				m_Z				= ( lep1_p4 + lep2_p4 ).M();
			}
			if ( cat6 ){
				lep1_p4			=	Muon_p4()[h_cand1[0]];
				lep1_pt			=	Muon_pt()[h_cand1[0]];
				lep1_eta		=	Muon_eta()[h_cand1[0]];
				lep1_eta_bdt	=	Muon_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi		=	Muon_phi()[h_cand1[0]];	
				lep1_charge		=	Muon_charge()[h_cand1[0]];
				lep1_pdgID		=	Muon_pdgId()[h_cand1[0]];
				lep1_tightID	=	Muon_tightId()[h_cand1[0]];

				lep2_p4			=	Electron_p4()[h_cand2[0]];
				lep2_pt			=	Electron_pt()[h_cand2[0]];
				lep2_eta		=	Electron_eta()[h_cand2[0]];
				lep2_eta_bdt	=	Electron_eta()[h_cand2[0]] * sgn( gg_eta );
				lep2_phi		=	Electron_phi()[h_cand2[0]];	
				lep2_charge		=	Electron_charge()[h_cand2[0]];
				lep2_tightID	=	Electron_mvaFall17V2Iso_WP90()[h_cand2[0]];
				lep2_pdgID		=	Electron_pdgId()[h_cand2[0]];

				lep12_dr		= deltaR( Electron_p4()[h_cand2[0]], Muon_p4()[h_cand1[0]] );
			}
			if ( cat7 ){
				lep1_p4					=	Tau_p4()[h_cand1[0]];
				lep1_pt					=	Tau_pt()[h_cand1[0]];
				lep1_eta				=	Tau_eta()[h_cand1[0]];
				lep1_eta_bdt			=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi				=	Tau_phi()[h_cand1[0]];	
				lep1_charge				=	Tau_charge()[h_cand1[0]];
				lep1_pdgID				=	15 * sgn( lep2_charge );
				lep1_id_vs_e			=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
				lep1_id_vs_m			=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
				lep1_id_vs_jet			=	 (int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];

				lep2_pt					=	IsoTrack_pt()[h_cand2[0]];
				lep2_eta				=	IsoTrack_eta()[h_cand2[0]];
				lep2_eta_bdt			=	IsoTrack_eta()[h_cand2[0]] * sgn( gg_eta );
				lep2_phi				=	IsoTrack_phi()[h_cand2[0]];	
				lep2_charge				=	IsoTrack_pdgId()[h_cand2[0]]/fabs(IsoTrack_pdgId()[h_cand2[0]]);
				lep2_pdgID				=	IsoTrack_pdgId()[h_cand2[0]];
				lep2_pfRelIso03_all 	= 	IsoTrack_pfRelIso03_all()[h_cand2[0]];;
				lep2_pfRelIso03_chg 	= 	IsoTrack_pfRelIso03_chg()[h_cand2[0]];;

				LorentzVector iso_track(IsoTrack_pt()[h_cand2[0]], IsoTrack_eta()[h_cand2[0]], IsoTrack_phi()[h_cand2[0]], 0);
				lep2_p4 = iso_track;
				lep12_dr		= deltaR( iso_track, Tau_p4()[h_cand1[0]] );
			}

			if ( cat8 ){
				lep1_pt			=	Tau_pt()[h_cand1[0]];
				lep1_eta		=	Tau_eta()[h_cand1[0]];
				lep1_eta_bdt	=	Tau_eta()[h_cand1[0]] * sgn( gg_eta );
				lep1_phi		=	Tau_phi()[h_cand1[0]];	
				lep1_charge		=	Tau_charge()[h_cand1[0]];
				lep1_pdgID		=	15 * sgn( lep1_charge );
				lep1_id_vs_e	=	 (int)Tau_idDeepTau2017v2p1VSe()[h_cand1[0]];	
				lep1_id_vs_m	=	 (int)Tau_idDeepTau2017v2p1VSmu()[h_cand1[0]];	
				lep1_id_vs_jet	=	(int)Tau_idDeepTau2017v2p1VSjet()[h_cand1[0]];
			}

			if ( lep1_pt > lep2_pt ){
				max_lep_pt = lep1_pt;
				min_lep_pt = lep2_pt;
			}
			else{
				max_lep_pt = lep2_pt;
				min_lep_pt = lep1_pt;
			}

			dPhi_MET_l	= deltaPhi( t_MET_phi, lep1_phi );

			//remove main ZGamma bkg in 0tau2lep 
			if ( category >3 && category < 7 ){
				m_llg_lead	= ( lep1_p4 + lep2_p4 + Photon_p4().at(gHidx[0]) ).M();
				m_llg_subl	= ( lep1_p4 + lep2_p4 + Photon_p4().at(gHidx[1]) ).M();
			}
			if ( fabs( m_llg_lead - mZ ) < mllg_window | fabs( m_llg_subl - mZ ) < mllg_window  ) continue;

			if ( sel_jets.size() > 0 ){
				jet1_pt		=	Jet_pt()[sel_jets[0]];
				jet1_eta	=	Jet_eta()[sel_jets[0]];
				jet1_eta_bdt=	jet1_eta* sgn( gg_eta );
				jet1_phi	=	Jet_phi()[sel_jets[0]];
				jet1_bTag	=	Jet_btagDeepFlavB()[sel_jets[0]];
				jet1_id		=	Jet_jetId()[sel_jets[0]];
			}
			if ( sel_jets.size() > 1 ){
				jet2_pt		=	Jet_pt()[sel_jets[1]];
				jet2_eta	=	Jet_eta()[sel_jets[1]];
				jet2_eta_bdt=	jet2_eta* sgn( gg_eta );
				jet2_phi	=	Jet_phi()[sel_jets[1]];
				jet2_bTag	=	Jet_btagDeepFlavB()[sel_jets[1]];
				jet2_id		=	Jet_jetId()[sel_jets[1]];
			}

			//make histograms for yields
			h_mgg->Fill( mgg, weight );
			if ( cat1 || cat2 ) h_mgg_1t1l->Fill( mgg, weight );
			if ( cat3 ) h_mgg_2t0l->Fill( mgg, weight );
			if ( cat4 || cat5 || cat6 ) h_mgg_0t2l->Fill( mgg, weight );
			if ( cat7 ) h_mgg_1t0l_iso->Fill( mgg, weight );
			if ( cat8 ) h_mgg_1t0l->Fill( mgg, weight );

       } // Event loop
       delete file;
    } // File loop
    bar.finish();

	f1->Write();
	f1->Close();
	return 0;
}
