import os
import glob
import pandas
import argparse
import root_pandas
import numpy as np

#years = [ '2016', '2017', '2018' ]
years = { b'2016UL_pos':'2016', b'2016UL_pre':'2016', b'2017':'2017', b'2018':'2018' }
procs_dict = { "ggH_M125": "ggH", 'HHggTauTau':'HHggTauTau', 'HHggWW_dileptonic':'HHggWWdileptonic', 'HHggWW_semileptonic':'HHggWWsemileptonic', 'ttH_M125':'ttH', 'VBFH_M125':'VBFH', 'VH_M125':'VH', "data":"Data" }

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input parquet directory",
    type = str,
    default = "/ceph/cms/store/user/fsetti/HiggsDNA_output/09Apr2022_sr_resonant_systs/"
)

parser.add_argument(
    "--tag",
    help = "unique tag to identify batch of processed samples",
    type = str,
    default = "10Apr2022"
)
parser.add_argument(
    "--mvas",
	nargs='*',
    help = "mva limits to SRs",
    type = float,
    default = [0.971078, 0.9928]	
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

out_dir = '/home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files_systs/' + str(args.tag) + '/'

#os.system("rm -rf %s"%(out_dir))
os.system("mkdir -p %s"%(out_dir))

os.system("mkdir -p %s/Data"%(out_dir))
os.system("mkdir -p %s/2016"%(out_dir))
#os.system("mkdir -p %s/2016_APV"%(out_dir))
os.system("mkdir -p %s/2017"%(out_dir))
os.system("mkdir -p %s/2018"%(out_dir))

procs = glob.glob(str(args.input)+'/*')

for proc in procs[:]:

	#if "data" not in proc:
	#	continue
	
	#get all files including systematic variations
	files = glob.glob(proc+'/*.parquet')
	
	for file_ in files:
		print ("Now processing: ", file_)
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
		if 'scale' in file_.split("/")[-1]:
				tag = '_MCScale' + tag

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
		#print(rename_sys)
		df = df.rename(columns=rename_sys)

		proc_tag = procs_dict[proc.split("/")[-1].split("_local_")[-1]]

		'''
		Temporary fix to get limits
		'''
		#Process Data
		if "Data" in proc_tag:
			for sr in range(args.nSRs):
				dfs = df.loc[ (df.process_id == 0 ) & ( df.bdt_score < args.mvas[sr] ) & ( df.bdt_score >= args.mvas[sr+1] ) & ( ( df.Diphoton_mass < 120 ) | ( df.Diphoton_mass > 130 ) ) ]
				dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR'+str(sr+1), mode='a')

		else:		
			#Process Signal
			for year,year_str in years.items():
				for sr in range(args.nSRs):
					if '2016' not in year_str :
						dfs = df.loc[ (df.year == year ) & ( df.bdt_score < args.mvas[sr] ) & ( df.bdt_score >= args.mvas[sr+1] )  ]
						dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR'+str(sr+1)+tag, mode='a')
					elif b'2016UL_pos' not in year:
						dfs = df.loc[ ( (df.year == b'2016UL_pre') | (df.year == b'2016UL_pos') ) & ( df.bdt_score < args.mvas[sr] ) & ( df.bdt_score >= args.mvas[sr+1] )  ]
						dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR'+str(sr+1)+tag, mode='a')

		#if "Data" in proc_tag:
		#	#Process Data
		#	dfs = df.loc[ ( df.pass_sr_0 ) ]
		#	dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR1'+tag, mode='a')
		#	dfs = df.loc[ ( df.pass_sr_1 ) ]
		#	dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR2'+tag, mode='a')

		#else:		
		#	#Process Signal- Terrible naming of years :( 
		#	for year,year_str in years.items():
		#		if '2016' not in year_str :
		#			dfs = df.loc[ (df.year == year ) & ( df.pass_sr_0 )  ]
		#			dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR1'+tag, mode='a')
		#			dfs = df.loc[ (df.year == year ) & ( df.pass_sr_1 )  ]
		#			dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR2'+tag, mode='a')
		#		elif b'2016UL_pos' not in year:
		#			dfs = df.loc[ ( (df.year == b'2016UL_pre') | (df.year == b'2016UL_pos') ) & ( df.pass_sr_0 )  ]
		#			dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR1'+tag, mode='a')
		#			dfs = df.loc[ ( (df.year == b'2016UL_pre') | (df.year == b'2016UL_pos') ) & ( df.pass_sr_1 )  ]
		#			dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR2'+tag, mode='a')
