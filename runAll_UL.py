#!/usr/bin/env python
import os
import sys
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

r.gSystem.Load('NanoCORE/libTauAnalysis_ClassicSVfit.so')
r.gSystem.Load('loopers/loop_UL_C.so')

lumi = { "2016" : 16.51, "2016_APV" : 19.39, "2017" : 41.5, "2018" : 59.8 }
#years = ['2016', '2016_APV', '2017', '2018']
years = ['2017', '2018']
samples = {}


#with open('samples_and_scale1fb_UL_nanoAODv9.json', "r") as f_in:
with open('samples_and_scale1fb_ul_fullData_v1.json', "r") as f_in:
	samples = json.load(f_in)

for name, sample in samples.items()[:]:
	if 'Data' not in name :
		continue
	for year in years:
		print 'Start processing ', year, ' ' , str(name)
		ch = r.TChain("Events")
		list_of_files = []
		try:
			for path in sample[year]['paths']:
				list_of_files += glob.glob(path+'/*/*/*/*.root')
				list_of_files += glob.glob(path+'/*.root')
			list_of_files = [ x for x in list_of_files if '.root' in x ]
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
