#!/usr/bin/env python
import os
import sys
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

r.gSystem.Load('NanoCORE/libTauAnalysis_ClassicSVfit.so')
#r.gSystem.Load('loopers/loop_C.so')
r.gSystem.Load('loopers/test_C.so')

lumi = { "2016" : 35.9, "2017" : 41.5, "2018" : 59.8 }
years = ['2016', '2017', '2018']
#years = [ '2016']
samples = {}

corrupted_files = [ '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/EGamma_Run2018A_private_data18/skimNano-HggHtautauselection__v5/210223_101051/0000/tree_298.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_16.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_30.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_31.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_32.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_35.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_47.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016H_private_data16/skimNano-HggHtautauselection__v5/210223_100135/0000/tree_28.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8_private_mc17/skimNano-Hggselection__v5/210215_223313/0000/tree_135.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_16.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_30.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_31.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_32.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_35.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016G_private_data16/skimNano-Hggselection__v5/210215_163513/0000/tree_47.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016H_private_data16/skimNano-Hggselection__v5/210215_163626/0000/tree_28.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018A_private_data18/skimNano-Hggselection__v5/210215_164425/0000/tree_298.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018A_private_data18/skimNano-Hggselection__v5/210215_164425/0000/tree_71.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018B_private_data18/skimNano-Hggselection__v5/210215_164541/0000/tree_23.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018D_private_data18/skimNano-Hggselection__v5/210215_164812/0000/tree_236.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018D_private_data18/skimNano-Hggselection__v5/210215_164812/0000/tree_64.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018D_private_data18/skimNano-Hggselection__v5/210215_164812/0000/tree_95.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018A_private_data18/skimNano-Hggselection__v5/210215_164425/0000/tree_123.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018D_private_data18/skimNano-Hggselection__v5/210215_164812/0000/tree_534.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016E_private_data16/skimNano-Hggselection__v5/210215_163245/0000/tree_3.root',
'/hadoop/cms/store/user/legianni/skimNano-Hggselection/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_private_mc18/skimNano-Hggselection__v5/210215_162125/0000/tree_6.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_private_mc16_2/skimNano-Hggselection__v5/210215_165810/0000/tree_1.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016C_private_data16/skimNano-Hggselection__v5/210215_163017/0000/tree_27.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2017F_private_data17/skimNano-Hggselection__v5/210215_164308/0000/tree_11.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018D_private_data18/skimNano-Hggselection__v5/210215_164812/0000/tree_21.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DoubleEG_Run2016C_private_data16/skimNano-Hggselection__v5/210215_163017/0000/tree_27.root' ]

list_of_files_tmp = [ '/hadoop/cms/store/user/legianni/skimNano-Hggselection/HHggtautau_Era2016_private_mc16/skimNano-Hggselection__v5/210216_085328/0000/tree_9.root' , '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018A_private_data18/skimNano-Hggselection__v5/210215_164425/0000/tree_20.root' , '/hadoop/cms/store/user/legianni/skimNano-Hggselection/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_private_mc18/skimNano-Hggselection__v5/210215_230611/0000/tree_2.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_private_mc17/skimNano-Hggselection__v5/210215_162239/0000/tree_33.root' ]

with open('samples_and_scale1fb.json', "r") as f_in:
	samples = json.load(f_in)

for name, sample in samples.items():
	for year in years:
		if name == "HH_ggZZ_2l2q":
			print 'Start processing ', year, ' ' , str(name)
			ch = r.TChain("Events")
			list_of_files = []
			try:
				for path in sample[year]['paths']:
					list_of_files += glob.glob(path+'/*.root')
				#list_of_files = [ x for x in list_of_files if 'tree' in x and x in list_of_files_tmp ]
				list_of_files = [ x for x in list_of_files if 'tree' in x and x not in corrupted_files ]
				for file_ in list_of_files:
					ch.Add(file_);
				if str(name) != 'Data' and ch.GetEntries() != 0 :
					scale_factor = sample[year]['metadata']['scale1fb'] * lumi[year]
					r.ScanChain(ch, str(name) , int(year) , scale_factor, bool(sample['resonant']) )
				elif ch.GetEntries() != 0 :
					scale_factor = 1
					r.ScanChain(ch, str(name) , int(year) , scale_factor )
			except:
				print name ,' ' , year , ' not available in the samples file or detected corrupted file. '
				print 'Exiting now. '
				sys.exit()

			continue
