#!/usr/bin/env python
import glob
import json
import ROOT as r
import numpy as np
import pandas as pd

procs = ['HH_ggTauTau']
#procs = ['GJets', 'TT', 'ZGamma', 'WGamma', 'DiPhoton', 'VH', 'ggH', 'ttH', 'VBFH', 'Data', 'HH_ggTauTau', 'HH_ggWW_semileptonic', 'HH_ggWW_dileptonic', 'HH_ggZZ_4l', 'HH_ggZZ_2l2q']
#bkgs = ['GJets', 'TT', 'ZGamma', 'WGamma', 'DiPhoton', 'VH', 'ggH', 'ttH', 'VBFH']

skims = [ "*_ReRecoCheck_*" ]
#skims = [ '*old*' , '*new*', '*old*2016*' , '*new*2016*', '*old*2017*' , '*new*2017*', '*old*2018*' , '*new*2018*']
for skim in skims:
	yields = {}
	yields['allBkg'] = { 'inclusive': 0. , '1tau0lep': 0., '1tau0lep_iso':0., '1tau1lep': 0., '2tau0lep': 0., '0tau2lep': 0. }
	for proc in procs:
			yields[proc] = { 'inclusive': 0. , '1tau0lep': 0., '1tau0lep_iso':0., '1tau1lep': 0., '2tau0lep': 0., '0tau2lep': 0. }
			files = glob.glob('outputs/'+proc+skim+'.root')
			for file_ in files:
				file_ = r.TFile(file_)
				h			= file_.Get('mgg')
				h_1t0l		= file_.Get('mgg_1t0l')
				h_1t0l_iso	= file_.Get('mgg_1t0l_iso')
				h_1t1l		= file_.Get('mgg_1t1l')
				h_2t0l		= file_.Get('mgg_2t0l')
				h_0t2l		= file_.Get('mgg_0t2l')
				try:
					yields[proc]['inclusive'] += h.GetSumOfWeights()
					yields[proc]['1tau0lep'] += h_1t0l.GetSumOfWeights()
					yields[proc]['1tau0lep_iso'] += h_1t0l_iso.GetSumOfWeights()
					yields[proc]['1tau1lep'] += h_1t1l.GetSumOfWeights()
					yields[proc]['2tau0lep'] += h_2t0l.GetSumOfWeights()
					yields[proc]['0tau2lep'] += h_0t2l.GetSumOfWeights()
					if proc in bkgs:
						yields['allBkg']['inclusive'] += h.GetSumOfWeights()
						yields['allBkg']['1tau0lep'] += h_1t0l.GetSumOfWeights()
						yields['allBkg']['1tau0lep_iso'] += h_1t0l_iso.GetSumOfWeights()
						yields['allBkg']['1tau1lep'] += h_1t1l.GetSumOfWeights()
						yields['allBkg']['2tau0lep'] += h_2t0l.GetSumOfWeights()
						yields['allBkg']['0tau2lep'] += h_0t2l.GetSumOfWeights()
					
				except:
					print ('File ' , file_.GetName() , ' failed to load weights' )
	
	df = pd.DataFrame.from_dict( yields, orient='index')
	df = df[['inclusive','1tau0lep','1tau0lep_iso','1tau1lep','2tau0lep','0tau2lep']]
	#df.sort_index()
	print(df.to_latex(index=True))

