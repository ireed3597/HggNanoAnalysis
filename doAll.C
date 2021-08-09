{
    gROOT->ProcessLine(".L ../NanoCORE/NANO_CORE.so");
    gROOT->ProcessLine(".L test.C+");
    //gROOT->ProcessLine(".L read_samples.C+");
    //TChain *ch = new TChain("Events");
	//TString proc;
	//Int_t year;

    //ch->Add("/hadoop/cms/store/user/legianni/skimNano-Hggselection/HHggtautau_Era2016_private_mc16/skimNano-Hggselection__v5/210216_085328/0000/tree_1.root");
	//proc = "HHggtautau";
	//year = 2016;
    //ScanChain(ch, proc, year );

}
