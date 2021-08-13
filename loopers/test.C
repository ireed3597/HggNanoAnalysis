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

#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Base.h"
#include "../NanoCORE/tqdm.h"
#include "../NanoCORE/utils.cc"

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

int ScanChain( TChain *ch, string proc, int year, float scale_factor = 1, bool resonant = false ) {

	TString file_name = proc + "_" +  std::to_string(year);

	fstream syncOut;
	syncOut.open("sync_"+file_name+".txt", ios::out);

	TFile* f1 = new TFile("outputs/" + file_name + ".root", "RECREATE");
	H1(mgg, 60, 100 , 180 );
	H1(mgg_1t0l, 60, 100 , 180 );
	H1(mgg_1t1l, 60, 100 , 180 );
	H1(mgg_2t0l, 60, 100 , 180 );
	H1(mgg_0t2l, 60, 100 , 180 );
	H1(lead_pho_pt, 60, 0 , 70 );
	H1(lead_pho_pt_1t0l, 60, 0 , 70 );
	H1(lead_pho_pt_1t1l, 60, 0 , 70 );
	H1(lead_pho_pt_2t0l, 60, 0 , 70 );
	H1(lead_pho_pt_0t2l, 60, 0 , 70 );

	TTree *out_tree	=	new TTree("Events","output tree");
	out_tree->Branch("year",&year,"year/I");
	out_tree->Branch("run",&t_run,"run/I");
	out_tree->Branch("lumiBlock",&t_lumiBlock,"lumiBlock/I");
	out_tree->Branch("event",&t_event,"event/I");
	out_tree->Branch("MET_pt",&t_MET_pt,"MET_pt/F");
	out_tree->Branch("MET_phi",&t_MET_phi,"MET_phi/F");
	out_tree->Branch("weight",&t_weight,"weight/F");

	out_tree->Branch("mgg",&mgg,"mgg/F");
	out_tree->Branch("n_electrons",&n_electrons,"n_electrons/I");
	out_tree->Branch("n_muons",&n_muons,"n_muons/I");
	out_tree->Branch("n_taus",&n_taus,"n_taus/I");
	out_tree->Branch("n_isoTrks",&n_isoTrks,"n_isoTrks/I");

	out_tree->Branch("g1_ptmgg"		,	&g1_ptmgg	,	"g1_ptmgg/F"		);
	out_tree->Branch("g1_pt"		,	&g1_pt		, 	"g1_pt/F"			);
	out_tree->Branch("g1_eta"		,	&g1_eta		, 	"g1_eta/F"			);	  
	out_tree->Branch("g1_phi"		,	&g1_phi		,  	"g1_phi/F"			);
	out_tree->Branch("g1_idmva"		,	&g1_idmva	,	"g1_idmva/F"		);
	out_tree->Branch("g1_pixVeto"	,	&g1_pixVeto	,	"g1_pixVeto/F"		);
	out_tree->Branch("g2_ptmgg"		,	&g2_ptmgg	,	"g2_ptmgg/F"		);
	out_tree->Branch("g2_pt"		,	&g2_pt		, 	"g2_pt/F"			);
	out_tree->Branch("g2_eta"		,	&g2_eta		, 	"g2_eta/F"			);	  
	out_tree->Branch("g2_phi"		,	&g2_phi		,  	"g2_phi/F"			);
	out_tree->Branch("g2_idmva"		,	&g2_idmva	,	"g2_idmva/F"		);
	out_tree->Branch("g2_pixVeto"	,	&g2_pixVeto	,	"g2_pixVeto/F"		);
	out_tree->Branch("gg_pt"		,	&gg_pt		, 	"gg_pt/F"			);
	out_tree->Branch("gg_eta"		,	&gg_eta		, 	"gg_eta/F"			);	  
	out_tree->Branch("gg_phi"		,	&gg_phi		,  	"gg_phi/F"			);
	out_tree->Branch("gg_dR"		,	&gg_dR		,	"gg_dR/F"			);

	out_tree->Branch("lep1_pt"			,	&lep1_pt				, 	"lep1_pt/F"					);	  
	out_tree->Branch("lep1_eta"			,	&lep1_eta				,  	"lep1_eta/F"				);
	out_tree->Branch("lep1_phi"			,	&lep1_phi				,	"lep1_phi/F"				);
	out_tree->Branch("lep1_charge"		,	&lep1_charge			, 	"lep1_charge/F"				);	  
	out_tree->Branch("lep1_pdgID"		,	&lep1_pdgID				,  	"lep1_pdgID/F"				);
	out_tree->Branch("lep1_tightID"		,	&lep1_tightID			,	"lep1_tightID/F"			);
	out_tree->Branch("lep1_id_vs_e"		,	&lep1_id_vs_e			, 	"lep1_id_vs_e/F"			);	  
	out_tree->Branch("lep1_id_vs_m"		,	&lep1_id_vs_m			,  	"lep1_id_vs_m/F"			);
	out_tree->Branch("lep1_id_vs_jet"	,	&lep1_id_vs_jet			,	"lep1_id_vs_jet/F"			);
	out_tree->Branch("lep2_pt"			,	&lep2_pt				, 	"lep2_pt/F"					);	  
	out_tree->Branch("lep2_eta"			,	&lep2_eta				,  	"lep2_eta/F"				);
	out_tree->Branch("lep2_phi"			,	&lep2_phi				,	"lep2_phi/F"				);
	out_tree->Branch("lep2_charge"		,	&lep2_charge			, 	"lep2_charge/F"				);	  
	out_tree->Branch("lep2_pdgID"		,	&lep2_pdgID				,  	"lep2_pdgID/F"				);
	out_tree->Branch("lep2_tightID"		,	&lep2_tightID			,	"lep2_tightID/F"			);
	out_tree->Branch("lep2_id_vs_e"		,	&lep2_id_vs_e			, 	"lep2_id_vs_e/F"			);	  
	out_tree->Branch("lep2_id_vs_m"		,	&lep2_id_vs_m			,  	"lep2_id_vs_m/F"			);
	out_tree->Branch("lep2_id_vs_jet"	,	&lep2_id_vs_jet			,	"lep2_id_vs_jet/F"			);

	out_tree->Branch("jet1_pt"					,	&jet1_pt		, 	"jet1_pt/F"					);	  
	out_tree->Branch("jet1_eta"					,	&jet1_eta		, 	"jet1_eta/F"				);	  
	out_tree->Branch("jet1_bTag"				,	&jet1_bTag		, 	"jet1_bTag/F"				);	  
	out_tree->Branch("jet1_id"					,	&jet1_id		, 	"jet1_id/F"					);	  
	out_tree->Branch("jet2_pt"					,	&jet2_pt		, 	"jet2_pt/F"					);	  
	out_tree->Branch("jet2_eta"					,	&jet2_eta		, 	"jet2_eta/F"				);	  
	out_tree->Branch("jet2_bTag"				,	&jet2_bTag		, 	"jet2_bTag/F"				);	  
	out_tree->Branch("jet2_id"					,	&jet2_id		, 	"jet2_id/F"					);	  

	out_tree->Branch("pt_tautauSVFitLoose"		,	&pt_tautauSVFitLoose	,	"pt_tautauSVFitLoose/F"			);	  
	out_tree->Branch("eta_tautauSVFitLoose"		,	&eta_tautauSVFitLoose	, 	"eta_tautauSVFitLoose/F"		);	  
	out_tree->Branch("phi_tautauSVFitLoose"		,	&phi_tautauSVFitLoose	, 	"phi_tautauSVFitLoose/F"		);	  
	out_tree->Branch("m_tautauSVFitLoose"		,	&m_tautauSVFitLoose		, 	"m_tautauSVFitLoose/F"			);	  
	out_tree->Branch("dR_tautauSVFitLoose"		,	&dR_tautauSVFitLoose	, 	"dR_tautauSVFitLoose/F"			);	  
	out_tree->Branch("dR_ggtautauSVFitLoose"	,	&dR_ggtautauSVFitLoose	, 	"dR_ggtautauSVFitLoose/F"		);	  
	out_tree->Branch("dPhi_MET_tau1"			,	&dPhi_MET_tau1			, 	"dPhi_MET_tau1/F"				);	  

	out_tree->Branch("m_tautau_vis"				,	&m_tautau_vis  			, 	"m_tautau_vis/F"				);	  
	out_tree->Branch("pt_tautau_vis"			,	&pt_tautau_vis 			, 	"pt_tautau_vis/F"				);	  
	out_tree->Branch("eta_tautau_vis"			,	&eta_tautau_vis			, 	"eta_tautau_vis/F"				);	  
	out_tree->Branch("phi_tautau_vis"			,	&phi_tautau_vis			, 	"phi_tautau_vis/F"				);	  

	int n_mt		= 0;
	int n_et		= 0;
	int n_tt		= 0;
	int n_mm		= 0;
	int n_em		= 0;
	int n_ee		= 0;
	int n_t			= 0;
	int n_it		= 0;
	int n_iso		= 0;

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
            tree->LoadTree(loop_event);

            nEventsTotal++;
            bar.progress(nEventsTotal, nEventsChain);

			//di-photon selection
			mgg = (float)(Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1]) ).M();
			if ( mgg < mgg_lower || mgg > mgg_upper ) continue;
			if ( (proc == "Data" || !resonant) && mgg > mgg_sideband_lower && mgg < mgg_sideband_upper ) continue;

			//photon selection
			if ( Photon_pt().at(gHidx()[0]) < pho_pt_cut 				|| Photon_pt().at(gHidx()[1]) < pho_pt_cut ) continue;
			if ( Photon_pt().at(gHidx()[0]) / mgg < lead_pt_mgg_cut 	|| Photon_pt().at(gHidx()[1]) / mgg < sublead_pt_mgg_cut ) continue;
			if ( Photon_mvaID().at(gHidx()[0]) < pho_idmva_cut 			|| Photon_mvaID().at(gHidx()[1]) < pho_idmva_cut ) continue;
			if ( Photon_electronVeto().at(gHidx()[0]) < pho_eveto_cut 	|| Photon_electronVeto().at(gHidx()[1]) < pho_eveto_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx()[0])) > pho_eta_cut 		|| fabs(Photon_eta().at(gHidx()[1])) > pho_eta_cut ) continue;
			if ( fabs(Photon_eta().at(gHidx()[0])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx()[0])) < trans_eta_high ) continue;
			if ( fabs(Photon_eta().at(gHidx()[1])) > trans_eta_low 		&& fabs(Photon_eta().at(gHidx()[1])) < trans_eta_high ) continue;

			//trigger requirements
			//if ( year == 2016 && !HLT_Diphoton30_18_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;
			//if ( (year == 2017 || year == 2018 ) && !HLT_Diphoton30_22_R9Id_OR_IsoCaloId_AND_HE_R9Id_Mass90() ) continue;

			vector<int> sel_eles;
			for(unsigned int i=0; i<nElectron(); i++){
				if (Electron_pt().at(i) > ele_pt && fabs(Electron_eta().at(i)) < ele_eta && ( fabs(Electron_eta().at(i)) < trans_eta_low || fabs(Electron_eta().at(i)) > trans_eta_high ) && fabs(Electron_dxy().at(i)) < ele_dxy && fabs(Electron_dz().at(i)) < ele_dz 
				&& ( (Electron_pfRelIso03_all().at(i) < ele_pfRelIso && Electron_mvaFall17V2noIso_WP90().at(i)) || (Electron_mvaFall17V2Iso_WP90().at(i) ) ) 
				&& deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > ele_dR_pho && deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > ele_dR_pho 
				){
					sel_eles.push_back(i);
				}
			}
			
			vector<int> sel_muons;
			for(unsigned int i=0; i<nMuon(); i++){
				if (Muon_pt().at(i) > muon_pt && fabs(Muon_eta().at(i)) < muon_eta && fabs(Muon_dxy().at(i)) < muon_dxy && fabs(Muon_dz().at(i)) < muon_dz 
				&& Muon_pfRelIso03_all().at(i) < muon_pfRelIso 
				&& Muon_isGlobal().at(i) 
				&& deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > muon_dR_pho && deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > muon_dR_pho 
				){
					sel_muons.push_back(i);
				}
			}
			
			vector<int> sel_taus;
			for(unsigned int i=0; i<nTau(); i++){
				if (Tau_pt().at(i) > tau_pt && fabs(Tau_eta().at(i)) < tau_eta && Tau_idDecayModeNewDMs().at(i) && fabs(Tau_dz().at(i)) < tau_dz 
				&& Tau_idDeepTau2017v2p1VSe().at(i) >= tau_deepID_e && Tau_idDeepTau2017v2p1VSmu().at(i) >= tau_deepID_m && Tau_idDeepTau2017v2p1VSjet().at(i) >= tau_deepID_j 
				&& deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > tau_dR_pho && deltaR( Tau_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > tau_dR_pho 
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
			if ( sel_taus.size() == 1 && sel_eles.size() + sel_muons.size() == 0 ){
	            for (unsigned int i=0; i<nIsoTrack(); i++){
	                if ( IsoTrack_isPFcand().at(i) && IsoTrack_fromPV().at(i) && IsoTrack_pdgId().at(i)*Tau_charge().at(sel_taus.at(0)) < 0 ){
	                    LorentzVector *iso_track = new LorentzVector;
	                    iso_track->SetXYZT( IsoTrack_pt().at(i)* TMath::Cos(IsoTrack_phi().at(i)) , IsoTrack_pt().at(i)*TMath::Sin( IsoTrack_phi().at(i)), IsoTrack_pt().at(i)*TMath::SinH( IsoTrack_eta().at(i)),  IsoTrack_pt().at(i)*TMath::CosH( IsoTrack_eta().at(i) ) );
	                    if ( deltaR_v1( iso_track , Tau_p4().at(sel_taus.at(0)) ) > isoTrk_dR && deltaR_v1( iso_track , Photon_p4().at(gHidx()[0]) ) > isoTrk_dR  && deltaR_v1( iso_track , Photon_p4().at(gHidx()[1]) ) > isoTrk_dR  ){
	                        sel_isoTracks.push_back( i );
	                    }
	                }
	            }
			}

			//Z veto cut
			bool Z_cand = false;
			if ( sel_eles.size() >= 2 )
				for (unsigned int i=0; i<sel_eles.size(); i++){
					for (unsigned int j=i+1; j<sel_eles.size(); j++){
						if ( (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() > mZ_veto_low  && (Electron_p4().at(sel_eles[i]) + Electron_p4().at(sel_eles[j])).M() < mZ_veto_up ){
							Z_cand = true;
							break;
					}
				}
			if ( sel_muons.size() >= 2 )
				for (unsigned int i=0; i<sel_muons.size(); i++){
					for (unsigned int j=i+1; j<sel_muons.size(); j++){
						if ( (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() > mZ_veto_low  && (Muon_p4().at(sel_muons[i]) + Muon_p4().at(sel_muons[j])).M() < mZ_veto_up ){
							Z_cand = true;
							break;
					}
				}
			}
			if ( Z_cand ) continue;

			n_electrons		= sel_eles.size();
			n_muons			= sel_muons.size();
			n_taus			= sel_taus.size();
			n_isoTrks		= sel_isoTracks.size();

			vector<int> h_cand1, h_cand2;
			vector<vector<int>> raw_results = categorise( sel_eles, sel_muons, sel_taus, sel_isoTracks );
			sel_eles		= raw_results[0];
			sel_muons		= raw_results[1];
			sel_taus		= raw_results[2];
			sel_isoTracks	= raw_results[3];
			h_cand1			= raw_results[4];
			h_cand2			= raw_results[5];

			if ( h_cand1[1] == 2 && h_cand2[1] == 1  ) cat1 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 0  ) cat2 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 2  ) cat3 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 1  ) cat4 = true;
			if ( h_cand1[1] == 0 && h_cand2[1] == 0  ) cat5 = true;
			if ( h_cand1[1] == 1 && h_cand2[1] == 0  ) cat6 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == 3  ) cat7 = true;
			if ( h_cand1[1] == 2 && h_cand2[1] == -1 ) cat8 = true;

			

			float weight = 1.;
			if ( proc != "Data" ) weight = genWeight() * scale_factor;

			h_mgg->Fill( mgg, weight );
			h_lead_pho_pt->Fill( Photon_pt().at(gHidx().at(0)), weight );

			t_run			= run();
			t_lumiBlock		= luminosityBlock();
			t_event			= event();

			g1_ptmgg		=	Photon_pt().at(gHidx()[0]) / mgg;
			g1_pt			=	Photon_pt().at(gHidx()[0]) ;
			g1_eta			=	Photon_eta().at(gHidx()[0]) ;
			g1_phi			=	Photon_phi().at(gHidx()[0]) ;
			g1_idmva		=	Photon_mvaID().at(gHidx()[0]) ;

			g2_ptmgg		=	Photon_pt().at(gHidx()[1]) / mgg;
			g2_pt			=	Photon_pt().at(gHidx()[1]) ;
			g2_eta			=	Photon_eta().at(gHidx()[1]) ;
			g2_phi			=	Photon_phi().at(gHidx()[1]) ;
			g2_idmva		=	Photon_mvaID().at(gHidx()[1]) ;

			gg_pt			=	(Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1])).pt() ;
			gg_eta			=	(Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1])).eta() ;
			gg_phi			=	(Photon_p4().at(gHidx()[0]) + Photon_p4().at(gHidx()[1])).phi() ;
			gg_dR			=	deltaR(Photon_p4().at(gHidx()[0]) , Photon_p4().at(gHidx()[1])) ;

			out_tree->Fill();

			//debugging
			int category;
			if ( cat1 ) category = 1;
			if ( cat2 ) category = 2;
			if ( cat3 ) category = 3;
			if ( cat4 ) category = 4;
			if ( cat5 ) category = 5;
			if ( cat6 ) category = 6;
			if ( cat7 ) category = 7;
			if ( cat8 ) category = 8;
			syncOut<<run()<<","<<luminosityBlock()<<","<<event()<<","<<","<<n_electrons<<","<<n_muons<<","<<n_taus<<","<<category<<std::endl;

			clear_branches();

        } // Event loop
        delete file;
    } // File loop
    bar.finish();

	cout << "category 1 mu-tau: " << n_mt << endl;
	cout << "category 2 ele-tau: " << n_et << endl;
	cout << "category 3 tau-tau: " << n_tt << endl;
	cout << "category 4 mu-mu: " << n_mm << endl;
	cout << "category 5 ele-ele: " << n_ee << endl;
	cout << "category 6 mu-ele: " << n_em << endl;
	cout << "category -1 tau: " << n_t << endl;
	cout << "category 1tau+Iso: " << n_it << endl;
	cout << "number of Iso: " << n_iso << endl;

	f1->Write();
	f1->Close();
	return 0;
}