import os
import uproot
import argparse
import root_pandas

ggHH_proc_ids = { '-10': "ggHHkl0kt1", '-11': "ggHHkl1kt1",'-12': "ggHHkl2p45kt1",'-13': "ggHHkl5kt1"}
#hh_coupl_proc_ids = { '-10': "ggHHkl0kt1", '-11': "ggHHkl1kt1",'-12': "ggHHkl2p45kt1",'-13': "ggHHkl5kt1",'-14': "qqHHCV0p5C2V1kl1",'-16': "qqHHCV1C2V0kl1",'-17': "qqHHCV1C2V1kl0",'-18': "qqHHCV1C2V1kl1",'-19': "qqHHCV1C2V1kl2",'-20': "qqHHCV1C2V2kl1",'-15': "qqHHCV1p5C2V1kl1"}
res_proc_ids = { '9': "VH", '10': "ttH", '11': "ggH", '12': "VBFH" }

#correct xsec of non-SM processes

#xsec_corr_bbbb  = { "ggHHkl0kt1" : 0.023618, "ggHHkl1kt1" : 0.010517, "ggHHkl2p45kt1": 0.004455, "ggHHkl5kt1": 0.031072, "qqHHCV0p5C2V1kl1": 0.003656, "qqHHCV1C2V0kl1": 0.009169, "qqHHCV1C2V1kl0": 0.001558, "qqHHCV1C2V1kl1": 0.000585, "qqHHCV1C2V1kl2": 0.000482, "qqHHCV1C2V2kl1": 0.004823, "qqHHCV1p5C2V1kl1": 0.004823}
xsec_corr_bbbb  = { "ggHHkl0kt1" : 70.38, "ggHHkl1kt1" : 31.05, "ggHHkl2p45kt1": 13.10, "ggHHkl5kt1": 94.82, "qqHHCV0p5C2V1kl1": 0.003656, "qqHHCV1C2V0kl1": 0.009169, "qqHHCV1C2V1kl0": 0.001558, "qqHHCV1C2V1kl1": 0.000585, "qqHHCV1C2V1kl2": 0.000482, "qqHHCV1C2V2kl1": 0.004823, "qqHHCV1p5C2V1kl1": 0.004823}
xsec_corr_ggtt = {}

for key, value in xsec_corr_bbbb.items():
  value_tmp = 0
  if "ggHH" in key:
    value_tmp = value / xsec_corr_bbbb["ggHHkl1kt1"]
  elif "qqHH" in key:
    value_tmp = value / xsec_corr_bbbb["qqHHCV1C2V1kl1"]
  xsec_corr_ggtt[key] = value_tmp

years = [ '2016', '2017', '2018' ]
	
parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input root file",
    type = str,
    default = "/home/users/fsetti/HHggTauTau/HggAnalysisDev/ttH_SR_Optimization/run2_kl1.root"
)
parser.add_argument(
    "--mvas",
	nargs='*',
    help = "mva limits to SRs",
    type = float,
    default = [0.978630,  0.9908]	
)
parser.add_argument(
    "--nSRs",
    help = "number of Signal Regions",
    type = int,
    default = 2
)

args = parser.parse_args()

args.mvas+=[99]
args.mvas.sort(reverse=True)

out_dir = '/home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/'

os.system("rm -rf /home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/*")

os.system("mkdir -p /home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/2016")
os.system("mkdir -p /home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/2017")
os.system("mkdir -p /home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/2018")
os.system("mkdir -p /home/users/fsetti/HHggTauTau/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files/Data")

t = uproot.open(args.input)['t']
df = t.pandas.df()

#df = df.loc[ abs( df.weight ) < 1e-4 ]

#Process HH signal different couplings
for proc_id, proc in ggHH_proc_ids.iteritems():
	for year in years:
		for sr in range(args.nSRs):
			dfs = df.loc[ (df.process_id == int(proc_id) ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) & (df.train_label == 2) ]
			#dfs["weight"] *= xsec_corr_ggtt[proc]			#rescale weights of non-SM couplings
			dfs.to_root(out_dir+year+'/'+proc+'_125_13TeV.root',proc+'_125_13TeV_SR'+str(sr+1), mode='a')


#t = uproot.open("/home/users/fsetti/HHggTauTau/HggAnalysisDev/ttH_SR_Optimization/run2_01Dec2021_final.root")['t']
#t = uproot.open(args.input)['t']
#df = t.pandas.df()

#Process Data
for sr in range(args.nSRs):
	dfs = df.loc[ (df.process_id == 0 ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] )  ]
	dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR'+str(sr+1), mode='a')

#Process Resonant bkg
for proc_id, proc in res_proc_ids.iteritems():
	for year in years:
		for sr in range(args.nSRs):
			dfs = df.loc[ (df.process_id == int(proc_id) ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) & (df.train_label == 2) ]
			dfs.to_root(out_dir+year+'/'+proc+'_125_13TeV.root',proc+'_125_13TeV_SR'+str(sr+1), mode='a')
