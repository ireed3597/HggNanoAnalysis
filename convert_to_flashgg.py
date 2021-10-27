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
    default = "hadded/run2_06Oct2021.root"
)
parser.add_argument(
    "--mvas",
	nargs='*',
    help = "mva limits to SRs",
    type = float,
    default = [0, 1]
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


for proc_id in proc_ids.keys():
	for year in years:
		for sr in range(args.nSRs):
			print 'Now processing: ' , proc_ids[proc_id] , ' for ' , year , ' and SR' , str(sr+1)
			dfs = df.loc[ (df.process_id == int(proc_id) ) & (df.year == int(year) ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) ]
			dfs.to_root('/home/users/fsetti/HggNanoAnalysis/flashgg_files/'+year+'/'+proc_ids[proc_id]+'_125_13TeV_SR'+str(sr+1)+'.root',key=proc_ids[proc_id]+'_125_13TeV_SR'+str(sr+1))
