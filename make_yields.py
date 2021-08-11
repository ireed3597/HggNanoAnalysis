#!/usr/bin/env python
import glob
import json
import ROOT as r
import numpy as np
import pandas as pd

#procs = ['GJets', 'TT', 'ZGamma', 'WGamma', 'DiPhoton', 'VH', 'ggH', 'ttH', 'VBFH', 'Data', 'HH_ggTauTau', 'HH_ggWW_semileptonic', 'HH_ggWW_dileptonic', 'HH_ggZZ']
procs = ['GJets', 'TT', 'ZGamma', 'WGamma', 'DiPhoton', 'VH', 'ggH', 'ttH', 'VBFH', 'Data', 'HH_ggTauTau', 'HH_ggWW_semileptonic', 'HH_ggWW_dileptonic']
yields = {}

for proc in procs:
		yields[proc] = { 'inclusive': 0. , '1tau0lep': 0., '1tau1lep': 0., '2tau0lep': 0., '0tau2lep': 0. }
		files = glob.glob('outputs/'+proc+'*.root')
		for file_ in files:
			file_ = r.TFile(file_)
			h		= file_.Get('mgg')
			h_1t0l	= file_.Get('mgg_1t0l')
			h_1t1l	= file_.Get('mgg_1t1l')
			h_2t0l	= file_.Get('mgg_2t0l')
			h_0t2l	= file_.Get('mgg_0t2l')
			yields[proc]['inclusive'] += h.GetSumOfWeights()
			yields[proc]['1tau0lep'] += h_1t0l.GetSumOfWeights()
			yields[proc]['1tau1lep'] += h_1t1l.GetSumOfWeights()
			yields[proc]['2tau0lep'] += h_2t0l.GetSumOfWeights()
			yields[proc]['0tau2lep'] += h_0t2l.GetSumOfWeights()

df = pd.DataFrame.from_dict( yields, orient='index')
df.sort_index()
print(df.to_latex(index=True))

