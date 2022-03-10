import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input root file",
    type = str,
    default = "pickles_zipped/run2_20UL_09Mar2022.pkl"
)
parser.add_argument(
    "--mvas",
	nargs='*',
    help = "mva limits to SRs",
    type = float,
    #default = [0.978630,  0.9908]		#gave 19.3000 x SM w/ flashgg , 18.6 - correct
    #default = [0.973863,  0.9901]		#now gives 19.75 x SM with private tools... worse than 18.6 above
    default = [0.975774,  0.9945]		#now gives 19.75 x SM with private tools... worse than 18.6 above
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

df = pd.read_pickle(args.input)

procs = { 'signal': [-6,-1], 'non_res_MC': [1,8], 'non_res_Data': [0,0], 'res_MC': [9,12] }
yields = {}

#Yields in SRs
for sr in range(args.nSRs):
	for proc, bnds in procs.items():
		yields[proc] = { 'inclusive': 0. , '1tau0lep_noIso':0., '1tau0lep_iso':0., '1tau1lep': 0., '2tau0lep': 0., '0tau2lep': 0. }
		dfs = df.loc[  (df.process_id >= bnds[0] ) & (df.process_id <= bnds[1] ) & ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] )  ]
		yields[proc]['inclusive'] = sum(dfs.weight)
		yields[proc]['1tau0lep_noIso'] = sum(dfs.loc[ dfs.Category == 8 ].weight)
		yields[proc]['1tau0lep_iso'] = sum(dfs.loc[ dfs.Category == 7 ].weight)
		yields[proc]['1tau1lep'] = sum(dfs.loc[ (dfs.Category == 1) | (dfs.Category == 2)  ].weight)
		yields[proc]['2tau0lep'] = sum(dfs.loc[ (dfs.Category == 3)  ].weight)
		yields[proc]['0tau2lep'] = sum(dfs.loc[ (dfs.Category == 4) | (dfs.Category == 5) | (dfs.Category == 6)  ].weight)

	y_df = pd.DataFrame.from_dict( yields, orient='index')
	y_df = y_df[['inclusive','2tau0lep','1tau1lep','0tau2lep','1tau0lep_iso','1tau0lep_noIso']]
	print("---------------------------------------------------------------------")
	print("SR", sr+1)
	print(y_df.to_latex(index=True))
	print("---------------------------------------------------------------------")


#Yields at Preselection
for proc, bnds in procs.items():
	yields[proc] = { 'inclusive': 0. , '1tau0lep_noIso':0., '1tau0lep_iso':0., '1tau1lep': 0., '2tau0lep': 0., '0tau2lep': 0. }
	dfs = df.loc[  (df.process_id >= bnds[0] ) & (df.process_id <= bnds[1] )  ]
	yields[proc]['inclusive'] = sum(dfs.weight)
	yields[proc]['1tau0lep_noIso'] = sum(dfs.loc[ dfs.Category == 8 ].weight)
	yields[proc]['1tau0lep_iso'] = sum(dfs.loc[ dfs.Category == 7 ].weight)
	yields[proc]['1tau1lep'] = sum(dfs.loc[ (dfs.Category == 1) | (dfs.Category == 2)  ].weight)
	yields[proc]['2tau0lep'] = sum(dfs.loc[ (dfs.Category == 3)  ].weight)
	yields[proc]['0tau2lep'] = sum(dfs.loc[ (dfs.Category == 4) | (dfs.Category == 5) | (dfs.Category == 6)  ].weight)

y_df = pd.DataFrame.from_dict( yields, orient='index')
y_df = y_df[['inclusive','2tau0lep','1tau1lep','0tau2lep','1tau0lep_iso','1tau0lep_noIso']]
print("---------------------------------------------------------------------")
print("Preselection")
print(y_df.to_latex(index=True))
print("---------------------------------------------------------------------")
