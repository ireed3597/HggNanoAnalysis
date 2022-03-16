{
    gROOT->ProcessLine(".L NanoCORE/libTauAnalysis_ClassicSVfit.so");
    gROOT->ProcessLine(".L loopers/loop_check_newSkim.C+");

		TChain *ch = new TChain("Events");
		ch->Add("/hadoop/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/GluGluToHHTo2G2W_semileptonic_node_SM_2017_final/skimNano-TestUL_GluGluToHHTo2G2W_semileptonic_node_SM_2017_final_TESTS/220225_220922/0000/*.root");
		ScanChain( ch, "HH_ggWW_semileptonic", "2017", 1, true);
}
