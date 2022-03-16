#!/usr/bin/env python
import os
import sys
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

r.gSystem.Load('NanoCORE/libTauAnalysis_ClassicSVfit.so')
#r.gSystem.Load('loopers/loop_add_BDT_vars_test_C.so')
r.gSystem.Load('loopers/loop_add_BDT_vars_C.so')

lumi = { "2016" : 17.95, "2016_APV" : 17.95, "2017" : 41.5, "2018" : 59.8 }
years = ['2016', '2016_APV', '2017', '2018']
#years = [ '2016']
samples = {}

os.system("mkdir -p outputs_UL/add_BDT_vars/")
os.system("mkdir -p outputs_UL/add_BDT_vars/hadded")

corrupted_files = [ "/hadoop/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final_TESTS/220308_012517/0000/tree_1.root", "/hadoop/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final_TESTS/220308_012517/0000/tree_134.root", "/hadoop/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final_TESTS/220308_012517/0000/tree_161.root", "/hadoop/cms/store/user/legianni/skimNano-TestUL__TEST-SamplesV9/DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final/skimNano-TestUL_DoubleEG_Run2016G-UL2016_MiniAODv2-v1_MINIAOD_final_TESTS/220308_012517/0000/tree_165.root" ]

with open('samples_and_scale1fb_UL_nanoAODv9.json', "r") as f_in:
	samples = json.load(f_in)

#switch = False

for name, sample in samples.items()[:]:
	#if 'ggf_' in name or "vbf_" in name:
	#	continue
	#if 'Data' in name :
	#	switch = True
	#	continue

	#if switch == True and "Data" not in name :
	if "HH_ggTauTau" in name :
		os.system("rm -rf outputs_UL/add_BDT_vars/%s"%(name))
		os.system("mkdir -p outputs_UL/add_BDT_vars/%s"%(name))
		for year in years:
			print 'Start processing ', year, ' ' , str(name)
			ch = r.TChain("Events")
			list_of_files = []
			try:
				for path in sample[year]['paths']:
					if name == "HH_ggZZ_2l2q":
						list_of_files += glob.glob(path+'/*.root')
					else :
						list_of_files += glob.glob(path+'/*/*/*/*.root')
				list_of_files = [ x for x in list_of_files if 'tree' in x and x not in corrupted_files ]
				for file_ in list_of_files[:]:
					ch.Add(file_);
				if str(name) != 'Data' and ch.GetEntries() != 0 :
					scale_factor = sample[year]['metadata']['scale1fb'] * lumi[year]
					r.ScanChain(ch, str(name) , year , scale_factor, bool(sample['resonant']) )
				if str(name) == 'Data' and ch.GetEntries() != 0 :
					scale_factor = 1
					r.ScanChain(ch, str(name) , year , scale_factor )
			except:
				print name ,' ' , year , ' not available in the samples file or detected corrupted file. '
				print 'Exiting now. '
				sys.exit()
			if year == '2016':
				continue
			elif year == '2016_APV':
				os.system("hadd -f -k outputs_UL/add_BDT_vars/hadded/%s_%s_hadded.root outputs_UL/add_BDT_vars/%s/%s*10Mar2022*_%s_*.root"%(name,year,name,name,'2016'))
			else:
				os.system("hadd -f -k outputs_UL/add_BDT_vars/hadded/%s_%s_hadded.root outputs_UL/add_BDT_vars/%s/%s*10Mar2022*_%s_*.root"%(name,year,name,name,year))
