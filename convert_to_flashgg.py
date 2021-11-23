import os
import uproot
import argparse
import root_pandas

proc_ids = { '-6': "HH_ggZZ_2l2q",  '-5': "HH_ggZZ_4l", '-4': "HH_ggWW_semileptonic", '-3': "HH_ggWW_dileptonic", '-2': "HH_ggZZ", '-1': "HH_ggTauTau", '0': 'Data', '2': "ZGamma", '3': "DiPhoton", '4': "WGamma", '5': "TTbar", '6': "TTGamma", '7': "TTGG", '8': "GJets", '9': "VH", '10': "ttH", '11': "ggH", '12': "VBFH" }

years = [ '2016', '2017', '2018' ]
	
parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input root file",
    type = str,
    default = "/home/users/fsetti/HHggTauTau/HggAnalysisDev/ttH_SR_Optimization/run2_19Nov2021.root"
)
parser.add_argument(
    "--mvas",
	nargs='*',
    help = "mva limits to SRs",
    type = float,
    default = [0.950229, 0.9800]
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

t = uproot.open(args.input)['t']
df = t.pandas.df()

out_dir = '/home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/'

os.system("rm -rf /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/*")

os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/Data")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2016")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2017")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2018")

for sr in range(args.nSRs):
	dfs = df.loc[ (df.process_id == 0 ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] )  ]
	dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR'+str(sr+1), mode='a')

for year in years:
	for sr in range(args.nSRs):
		dfs = df.loc[ (df.process_id < 0 ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) & (df.train_label == 2) ]
		dfs.to_root(out_dir+year+'/HHggTauTau_125_13TeV.root','ggf_125_13TeV_SR'+str(sr+1), mode='a')

#ggtt = [ 0, 0]
#ggWW = [ 0, 0]
#
#dft = df.loc[ (df.process_id == -1 ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[0] ) & ( df.mva_score >= args.mvas[1] ) & (df.train_label == 2) ]
#dfW = df.loc[ ( (df.process_id == -4 )  | (df.process_id == -3 )) & (df.year == int(year) ) & ( df.mva_score < args.mvas[0] ) & ( df.mva_score >= args.mvas[1] ) & (df.train_label == 2) ]
#
#ggtt[0] = sum(dft.weight)
#ggWW[0] = sum(dfW.weight)
#
#dft = df.loc[ (df.process_id == -1 ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[1] ) & ( df.mva_score >= args.mvas[2] ) & (df.train_label == 2) ]
#dfW = df.loc[ ( (df.process_id == -4 )  | (df.process_id == -3 )) & (df.year == int(year) ) & ( df.mva_score < args.mvas[1] ) & ( df.mva_score >= args.mvas[2] ) & (df.train_label == 2) ]
#
#ggtt[1] = sum(dft.weight)
#ggWW[1] = sum(dfW.weight)
#
#print ggtt
#print ggWW
