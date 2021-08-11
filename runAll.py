#!/usr/bin/env python
import os
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

r.gSystem.Load('NanoCORE/NANO_CORE.so')
r.gSystem.Load('loopers/test_C.so')

lumi = { "2016" : 35.9, "2017" : 41.5, "2018" : 59.8 }
years = ['2016', '2017', '2018']
#years = ['2018']
samples = {}

corrupted_files = [ '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/EGamma_Run2018A_private_data18/skimNano-HggHtautauselection__v5/210223_101051/0000/tree_298.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_16.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_30.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_31.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_32.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_35.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/skimNano-HggHtautauselection__v5/210223_100016/0000/tree_47.root', '/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016H_private_data16/skimNano-HggHtautauselection__v5/210223_100135/0000/tree_28.root' ]

list_of_files_tmp = [ '/hadoop/cms/store/user/legianni/skimNano-Hggselection/HHggtautau_Era2016_private_mc16/skimNano-Hggselection__v5/210216_085328/0000/tree_9.root' , '/hadoop/cms/store/user/legianni/skimNano-Hggselection/EGamma_Run2018A_private_data18/skimNano-Hggselection__v5/210215_164425/0000/tree_20.root' , '/hadoop/cms/store/user/legianni/skimNano-Hggselection/VHToGG_M125_13TeV_amcatnloFXFX_madspin_pythia8_private_mc18/skimNano-Hggselection__v5/210215_230611/0000/tree_2.root', '/hadoop/cms/store/user/legianni/skimNano-Hggselection/DiPhotonJetsBox_MGG-80toInf_13TeV-Sherpa_private_mc17/skimNano-Hggselection__v5/210215_162239/0000/tree_33.root' ]

with open('samples_and_scale1fb.json', "r") as f_in:
	samples = json.load(f_in)

for name, sample in samples.items():
	for year in years:
		print 'Start processing ', year, ' ' , str(name)
		ch = r.TChain("Events")
		list_of_files = []
		for path in sample[year]['paths']:
			list_of_files += glob.glob(path+'/*/*/*/*')
		list_of_files = [ x for x in list_of_files if 'tree' in x and x in list_of_files_tmp ]
		if name == 'Data':
			list_of_files = [ x for x in list_of_files if 'tree' in x and x not in corrupted_files ]
		for file_ in list_of_files:
			ch.Add(file_);
		if str(name) != 'Data':
			scale_factor = sample[year]['metadata']['scale1fb'] * lumi[year]
			r.ScanChain(ch, str(name) , int(year) , scale_factor, bool(sample['resonant']) )
		else:
			scale_factor = 1
			r.ScanChain(ch, str(name) , int(year) , scale_factor )
