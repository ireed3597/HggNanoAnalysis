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

vector<unsigned int> event_vec;

int ScanChain( TChain *ch, string proc, int year, float scale_factor = 1, bool resonant = false ) {

	TString file_name = proc + "_" +  std::to_string(year);
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

	int n_ele		= 0;
	int n_mu		= 0;

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

        for( unsigned int event = 0; event < tree->GetEntriesFast(); ++event) {

            nt.GetEntry(event);
            tree->LoadTree(event);

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

			vector<unsigned int> sel_eles;
			for(unsigned int i=0; i<nElectron(); i++){
				if (Electron_pt().at(i) > ele_pt && fabs(Electron_eta().at(i)) < ele_eta && ( fabs(Electron_eta().at(i)) < trans_eta_low || fabs(Electron_eta().at(i)) > trans_eta_high ) && fabs(Electron_dxy().at(i)) < ele_dxy && fabs(Electron_dz().at(i)) < ele_dz 
				&& ( (Electron_pfRelIso03_all().at(i) < ele_pfRelIso && Electron_mvaFall17V2noIso_WP90().at(i)) || (Electron_mvaFall17V2Iso_WP90().at(i) ) ) 
				&& deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > ele_dR_pho && deltaR( Electron_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > ele_dR_pho 
				){
					sel_eles.push_back(i);
				}
			}
			
			vector<unsigned int> sel_muons;
			for(unsigned int i=0; i<nMuon(); i++){
				if (Muon_pt().at(i) > muon_pt && fabs(Muon_eta().at(i)) < muon_eta && fabs(Muon_dxy().at(i)) < muon_dxy && fabs(Muon_dz().at(i)) < muon_dz 
				&& Muon_pfRelIso03_all().at(i) < muon_pfRelIso 
				&& Muon_isGlobal().at(i) 
				&& deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[0]) ) > muon_dR_pho && deltaR( Muon_p4().at(i) , Photon_p4().at(gHidx()[1]) ) > muon_dR_pho 
				){
					sel_muons.push_back(i);
				}
			}
			
			vector<unsigned int> sel_taus;
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
            vector<unsigned int> vec_isoTracks;
			if ( sel_taus.size() == 1 && sel_eles.size() && sel_muons.size() == 0 ){
	            for (unsigned int i=0; i<nIsoTrack(); i++){
	                if ( IsoTrack_isPFcand().at(i) && IsoTrack_fromPV().at(i) ){
	                    LorentzVector *iso_track = new LorentzVector;
	                    iso_track->SetXYZT( IsoTrack_pt().at(i)* TMath::Cos(IsoTrack_phi().at(i)) , IsoTrack_pt().at(i)*TMath::Sin( IsoTrack_phi().at(i)), IsoTrack_pt().at(i)*TMath::SinH( IsoTrack_eta().at(i)),  IsoTrack_pt().at(i)*TMath::CosH( IsoTrack_eta().at(i) ) );
	                    if ( deltaR_v1( iso_track , Tau_p4().at(sel_taus.at(0)) ) > 0.2 && deltaR_v1( iso_track , Photon_p4().at(gHidx()[0]) ) > 0.2  && deltaR_v1( iso_track , Photon_p4().at(gHidx()[1]) ) > 0.2 ){
	                        vec_isoTracks.push_back( i );
	                    }
	                }
	            }
			}

			n_electrons		= sel_eles.size();
			n_muons			= sel_muons.size();
			n_taus			= sel_taus.size();
			n_isoTrks		= vec_isoTracks.size();

			n_ele 	+= n_electrons;
			n_mu	+= n_muons;

			//require opposite sign leptons
			bool os_cut = true;
			//1tau1lep
			if ( n_taus == 1 && (n_electrons + n_muons == 1) && n_electrons == 1 && Tau_charge().at(sel_taus.at(0))*Electron_charge().at(sel_eles.at(0)) > 0 ) 	os_cut = false;
			if ( n_taus == 1 && (n_electrons + n_muons == 1) && n_muons == 1 && Tau_charge().at(sel_taus.at(0))*Muon_charge().at(sel_muons.at(0)) > 0 ) 		os_cut = false;
			//2tau0lep
			if ( n_taus == 2 && (n_electrons + n_muons == 0) && Tau_charge().at(sel_taus.at(0))*Tau_charge().at(sel_taus.at(1)) > 0 ) os_cut = false;
			//0tau2lep
			if ( n_taus == 0 && (n_electrons + n_muons == 2) && n_electrons == 2 && Electron_charge().at(sel_eles.at(0))*Electron_charge().at(sel_eles.at(1)) > 0 ) os_cut = false;
			if ( n_taus == 0 && (n_electrons + n_muons == 2) && n_muons == 2 && Muon_charge().at(sel_muons.at(0))*Muon_charge().at(sel_muons.at(1)) > 0 ) 	os_cut = false;
			if ( n_taus == 0 && (n_electrons + n_muons == 2) && n_electrons == n_muons && Electron_charge().at(sel_eles.at(0))*Muon_charge().at(sel_muons.at(0)) > 0 ) os_cut = false;
			if ( !os_cut ) continue;

			//Z veto cut
			bool mZ_veto = false;
			if ( 
				( n_electrons == 2 && (Electron_p4().at(sel_eles[0]) + Electron_p4().at(sel_eles[1])).M() > mZ_veto_low  && (Electron_p4().at(sel_eles[0]) + Electron_p4().at(sel_eles[1])).M() < mZ_veto_up )
			|| 	( n_muons == 2 && (Muon_p4().at(sel_muons[0]) + Muon_p4().at(sel_muons[1])).M() > mZ_veto_low  && (Muon_p4().at(sel_muons[0]) + Muon_p4().at(sel_muons[1])).M() < mZ_veto_up )
			||  ( n_electrons == 1 && n_muons == 1 && (Electron_p4().at(sel_eles[0]) + Muon_p4().at(sel_muons[0])).M() > mZ_veto_low  && (Electron_p4().at(sel_eles[0]) + Muon_p4().at(sel_muons[0])).M() < mZ_veto_up )
				){
				mZ_veto = true;
			}
			if ( mZ_veto ) continue;

			float weight = 1.;
			if ( proc != "Data" ) weight = genWeight() * scale_factor;

			h_mgg->Fill( mgg, weight );
			h_lead_pho_pt->Fill( Photon_pt().at(gHidx().at(0)), weight );

			if ( n_taus == 1 && (n_electrons + n_muons == 0 ) ){ 
				h_mgg_1t0l->Fill( mgg, weight );
				h_lead_pho_pt_1t0l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( n_taus == 1 && (n_electrons + n_muons == 1 ) ){
				h_mgg_1t1l->Fill( mgg, weight );
				h_lead_pho_pt_1t1l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( n_taus == 2 && (n_electrons + n_muons == 0 ) ){
				h_mgg_2t0l->Fill( mgg, weight );
				h_lead_pho_pt_2t0l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}
			if ( n_taus == 0 && (n_electrons + n_muons == 2 ) ){
				h_mgg_0t2l->Fill( mgg, weight );
				h_lead_pho_pt_0t2l->Fill( Photon_pt().at(gHidx().at(0)), weight );
			}


			t_run			= run();
			t_lumiBlock		= luminosityBlock();
			t_event			= event;

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

        } // Event loop
        delete file;
    } // File loop
    bar.finish();

	cout << "tot number of electrons: " << n_ele << endl;
	cout << "tot number of muons: " << n_mu << endl;

	f1->Write();
	f1->Close();
	return 0;
}
