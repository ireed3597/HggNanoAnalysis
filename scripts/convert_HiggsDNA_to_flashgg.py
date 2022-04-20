import os
import glob
import pandas
import argparse
import root_pandas
import numpy as np

#res_proc_ids = { '9': "VH", '10': "ttH", '11': "ggH", '12': "VBFH" }
sample_id_map = {
    "Data": 0,
    "HHggTauTau": 1,
    "VBFH": 5,
    "VH": 3,
    "ggH": 4,
    "ttH": 2
}


years = [ '2016', '2017', '2018' ]
	
parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input parquet directory",
    type = str,
    default = "test/"
)
args = parser.parse_args()

#get all files including systematic variations
files = glob.glob(str(args.input)+'*.parquet')
	
out_dir = '/home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/'

os.system("rm -rf /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/*")

os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/Data")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2016")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2017")
os.system("mkdir -p /home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files/2018")

for file_ in files:
	df = pandas.read_parquet(file_, engine='pyarrow')

	tag	=	''
	if 'nominal' not in file_.split("/")[-1]:
		if 'up' in file_.split("/")[-1]:
			tag	= file_.split("merged")[-1]
			tag	= tag.split("_up")[0]
			tag	+= 'Up01sigma'
		if 'down' in file_.split("/")[-1]:
			tag	= file_.split("merged")[-1]
			tag	= tag.split("_down")[0]
			tag	+= 'Down01sigma'

	#Define hgg_mass & dZ variable
	df['CMS_hgg_mass'] = df['Diphoton_mass']
	df['dZ']	= np.ones(len(df['Diphoton_mass']))
	df['weight'] = df['weight_central']
	yield_systematics	= [ key for key in df.keys() if ( "weight_" in key ) and ( "_up" in key or "_down" in key )]
	rename_sys	= {}
	for sys in yield_systematics:
		#a bit of gymnastics to get the inputs right for Mr. flashggFinalFit
		sys_central = sys.replace("_up","_central")
		sys_central = sys.replace("_down","_central")
		df[sys] 		= df[sys] / df[sys_central]
		if "_up" in sys:
			rename_sys[sys] = sys.replace("_up","Up01sigma")
		if "_down" in sys:
			rename_sys[sys] = sys.replace("_down","Down01sigma")
	print(rename_sys)
	df = df.rename(columns=rename_sys)

	#Process Data
	dfs = df.loc[ (df.process_id == 0 )  ]
	dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR1'+tag, mode='a')
	
	#Process Signal
	for year in years:
		#dfs = df.loc[ (df.process_id < 0 ) & (df.year == int(year) ) ]
		dfs = df.loc[ (df.process_id == 1. ) & (df.year == int(year) ) ]
		dfs.to_root(out_dir+year+'/HHggTauTau_125_13TeV.root','HHggTauTau_125_13TeV_SR1'+tag, mode='a')
	
	#Process Resonant bkg
	#for proc_id, proc in res_proc_ids.items():
	for proc, proc_id in sample_id_map.items():
		if proc_id <= 1.:
			continue
		for year in years:
			#dfs = df.loc[ (df.process_id == int(proc_id) ) & (df.year == int(year) ) ]
			dfs = df.loc[ (df.process_id == int(proc_id) ) & (df.year == int(year) ) ]
			dfs.to_root(out_dir+year+'/'+proc+'_125_13TeV.root',proc+'_125_13TeV_SR1'+tag, mode='a')
