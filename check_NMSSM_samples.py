#!/usr/bin/env python
import os
import sys
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

r.gSystem.Load('NanoCORE/libTauAnalysis_ClassicSVfit.so')
r.gSystem.Load('loopers/check_NMSSM_samples_C.so')

lumi = { "2016" : 35.9, "2017" : 41.5, "2018" : 59.8 }
#years = ['2016', '2017', '2018']
years = [ '2016', '2017']
samples = {}

ch = r.TChain("Events")
list_of_files = []

list_of_files += glob.glob('/hadoop/cms/store/user/fsetti/nanoAOD_runII_UltraLegacy/NMSSM_XYH_Y_tautau_H_gg_MX_300_MY_70_STEP7_v1/*.root')
for file_ in list_of_files:
	ch.Add(file_);
r.ScanChain(ch, 'NMSSM' , 2018 , 1  )
