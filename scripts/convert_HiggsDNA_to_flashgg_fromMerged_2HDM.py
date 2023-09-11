import json
import os
import glob
import pandas
import argparse
import root_pandas
import numpy as np
import uproot
from collections import Counter

# def to_root(df,path,key='my_tree',mode='w',store_index=True):

#   column_name_counts = Counter(df.columns)
#   if max(column_name_counts.values()) > 1:
#       raise ValueError('DataFrame contains duplicated column names: ' +
#                         ' '.join({k for k, v in column_name_counts.items() if v > 1}))

#   # We don't want to modify the user's DataFrame here, so we make a shallow copy
#   df_ = df.copy(deep=False)

#   if store_index:
#     name = df_.index.name
#     if name is None:
#         # Handle the case where the index has no name
#         name = ''
#     df_['__index__' + name] = df_.index

#   # Convert categorical columns into something root_numpy can serialise
#   for col in df_.select_dtypes(['category']).columns:
#     name_components = ['__rpCaT', col, str(df_[col].cat.ordered)]
#     name_components.extend(df_[col].cat.categories)
#     if ['*' not in c for c in name_components]:
#         sep = '*'
#     else:
#         raise ValueError('Unable to find suitable separator for columns')
#     df_[col] = df_[col].cat.codes
#     df_.rename(index=str, columns={col: sep.join(name_components)}, inplace=True)

#   arr = df_.to_records (index=False)

#   if mode == 'a':
#     root_file = uproot.update(path)
#   elif mode == 'w':
#     root_file = uproot.recreate(path)
#   else:
#       raise ValueError('Unknown mode: {}. Must be "a" or "w".'.format(mode))

#   if not root_file:
#     raise IOError("cannot open file {0}".format(path))
#   if not root_file.IsWritable():
#     raise IOError("file {0} is not writable".format(path))

#   # Navigate to the requested directory
#   open_dirs = [root_file]
#   for dir_name in key.split('/')[:-1]:
#       current_dir = open_dirs[-1].Get(dir_name)
#       if not current_dir:
#           current_dir = open_dirs[-1].mkdir(dir_name)
#       current_dir.cd()
#       open_dirs.append(current_dir)

#   # The key is now just the top component
#   key = key.split('/')[-1]

#   # If a tree with that name exists, we want to update it \\ You are dealing with uproot now...
#   tree = open_dirs[-1].Get(key)
#   if not tree:
#       tree = None
#   tree = array2tree(arr, name=key, tree=tree)
#   tree.Write(key, ROOT.TFile.kOverwrite)
#   root_file.Close()




# years = [ '2016UL_preVFP', '2017', '2018' ]
years = [ b'2016UL_pre', b'2017', b'2018' ]
procs_dict = {"Data":"Data",
              "DiPhoton":"DiPhoton",
              "DataDrivenGJets":"DataDrivenGJets",
              "GJets_HT-40To100":"GJetsHT40To100",
              "GJets_HT-100To200":"GJetsHT100To200",
              "GJets_HT-200To400":"GJetsHT200To400",
              "GJets_HT-400To600":"GJetsHT400To600",
              "GJets_HT-600ToInf":"GJetsHT600ToInf",
              "TTGG":"TTGG",
              "TTGamma":"TTGamma",
              "TTJets":"TTJets",
              "VBFH_M125":"VBFH",
              "VH_M125":"VH",
              "WGamma":"WGamma",
              "ZGamma":"ZGamma",
              "ggH_M125":"ggH",
              "ttHH_ggTauTau":"ttHHggTauTau",
              "ttHH_ggWW":"ttHHggWW",
              "ttHH_ggbb":"ttHHggbb",
              "ttH_M125":"ttH",
              "HHggTauTau":"HHGGTauTau",
              "HHggWW_dileptonic":"HHGGWWdileptonic",
              "HHggWW_semileptonic":"HHGGWWsemileptonic",
              "HHggbb":"HHGGbb",
              "2HDM_bb_M250":"2HDMbbM250",
              "2HDM_WW_M250":"2HDMWWM250",
              "2HDM_TAUTAU_M250":"2HDMTAUTAUM250"
            }
skip_list = ["DiPhoton",
             "DataDrivenGJets",
             "GJets_HT-40To100",
             "GJets_HT-100To200",
             "GJets_HT-200To400",
             "GJets_HT-400To600",
             "GJets_HT-600ToInf",
             "TTGG",
             "TTGamma",
             "TTJets",
             "WGamma",
             "ZGamma",
            ]

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input parquet directory",
    type = str,
    default = '/home/users/iareed/HiggsDNA/BSM_13Oct22/scored_M250/' 
)

parser.add_argument(
    "--tag",
    help = "unique tag to identify batch of processed samples",
    type = str,
    default = "BSM_13Oct22_M250_div3"
)
parser.add_argument(
    "--mvas",
  nargs="*",
    help = "mva limits to SRs",
    type = float,
    default = [0.957208,0.9897]
)
parser.add_argument(
    "--nSRs",
    help = "number of Signal Regions",
    type = int,
    default = 2
)

#TODO here: we want to open 1 big file, take the all the data (process_id==0) and the rest of the processes divided
# by year. Use summary.json for mapping process name to process_id

args = parser.parse_args()

args.mvas+=[99]
args.mvas.sort(reverse=True)

out_dir = '/home/users/iareed/ttHHggbb/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files_systs/' + str(args.tag) + '/'

os.system("rm -rf %s"%(out_dir))
os.system("mkdir -p %s"%(out_dir))

os.system("mkdir -p %s/Data"%(out_dir))
os.system("mkdir -p %s/2016"%(out_dir))
os.system("mkdir -p %s/2017"%(out_dir))
os.system("mkdir -p %s/2018"%(out_dir))

with open(str(args.input)+'/summary.json',"r") as f_in:
    procs_id_map = json.load(f_in)
procs = procs_id_map["sample_id_map"]
print ("procs {}".format(procs))

#Process Data
file = pandas.read_parquet(str(args.input)+'merged_nominal.parquet', engine='pyarrow') # no systematics on Data
print(file)
df=file[ (file.process_id == procs['Data']) ]
df['CMS_hgg_mass'] = df.Diphoton_mass
for sr in range(args.nSRs):
    print(args.mvas[sr], args.mvas[sr+1])
    dfs = df.loc[ ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) & ( ( df.Diphoton_mass < 120 ) | ( df.Diphoton_mass > 130 ) ) ]
    print("Adding {} events to allData".format(len(dfs)))
    dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR'+str(sr+1), mode='a')

#Process MCs
# I think we need to first open the files, then loop over the process and years
# get all files including systematic variations
#files = glob.glob(str(args.input)+'/merged_nominal.parquet')
files = glob.glob(str(args.input)+'/*.parquet')
print(files)
for file_ in files:
    glob_df = pandas.read_parquet(file_, engine='pyarrow')

    print ("Now processing: ", file_)
    for proc in procs.keys():
        if proc in skip_list:
            continue
        if "Data" in proc:
            continue

        if proc not in procs_dict.keys():
            print("\n\n WARNING Process {} was not found in the process dictionary, Skipping...".format(proc.split("/")[-1].split("_201")[0]) )
            continue

        proc_df = glob_df[glob_df["process_id"]==procs[proc]]
        print("for all years we have {} events for {}".format(len(proc_df),proc))
        for year in years:
            df = proc_df[proc_df["year"]==year]
            print("for year {} we have {} events for {}".format(year,len(df),proc))
            if year == b'2016UL_pre':
                try:
                    df_ext1 = proc_df[proc_df["year"]==b'2016UL_pos']
                    print('for year 2016UL_postVFP we have {} events for {}'.format(len(df_ext1),proc))
                    df = pandas.concat([ df, df_ext1 ], ignore_index=True)
                except:
                    print ("Not finding 2016UL_postVFP for this sample: ", proc )
                    print ("Most likely it is data")

            tag  =  ''
            if 'nominal' not in file_.split("/")[-1]:
                if 'up' in file_.split("/")[-1]:
                    tag  = file_.split("merged")[-1]
                    tag  = tag.split("_up")[0]
                    tag  += 'Up01sigma'
                if 'down' in file_.split("/")[-1]:
                    tag  = file_.split("merged")[-1]
                    tag  = tag.split("_down")[0]
                    tag  += 'Down01sigma'
            if 'scale' in file_.split("/")[-1]:
                tag = '_MCScale' + tag
            if 'smear' in file_.split("/")[-1]:
                tag = '_MCSmear' + tag
            #breakpoint()
            #Define hgg_mass & dZ variable
            df['CMS_hgg_mass'] = df['Diphoton_mass']
            df['dZ'] = np.ones(len(df['Diphoton_mass']))
            df['weight'] = df['weight_central'] * 2
            df['weight_central'] = df['weight']
            yield_systematics = [ key for key in df.keys() if ( "weight_" in key ) and ( "_up" in key or "_down" in key )]
            rename_sys = {}
            for sys in yield_systematics:
                #a bit of gymnastics to get the inputs right for Mr. flashggFinalFit
                if "_up" in sys:
                    sys_central = sys.replace("_up","_central")
                elif "_down" in sys:
                    sys_central = sys.replace("_down","_central")
                rename_sys[sys] = sys
                if 'btag' in sys_central:
                    sys_central = 'weight_btag_deepjet_sf_SelectedJet_central'
                df[sys] = df[sys] / df[sys_central]
                if sys.endswith("_lf"):
                    rename_sys[sys] = rename_sys[sys].replace("_lf","_LF")
                elif sys.endswith("_hf"):
                    rename_sys[sys] = rename_sys[sys].replace("_hf","_HF")
                if "_up" in sys:
                    if 'btag' in sys:
                        rename_sys[sys] = rename_sys[sys].replace("_up","")
                        rename_sys[sys] += "Up01sigma"
                        continue
                    rename_sys[sys] = sys.replace("_up","Up01sigma")
                if "_down" in sys:
                    if 'btag' in sys:
                        rename_sys[sys] = rename_sys[sys].replace("_down","")
                        rename_sys[sys] += "Down01sigma"
                        continue
                    rename_sys[sys] = sys.replace("_down","Down01sigma")
            df = df.rename(columns=rename_sys)
            #Process Signal
            year_str = ''
            if year == b'2016UL_pre':
                year_str = '2016'
            if year == b'2017':
                year_str = '2017'
            if year == b'2018':
                year_str = '2018'

            proc_tag = ''
            # if 'EFT' in proc.split("/")[-1]:
            #   proc_tag = 'ggHHbm' + proc.split("/")[-1].split("node_")[-1].split("_201")[0]
            # else:
            proc_tag = procs_dict[proc.split("/")[-1].split("_201")[0]]

            '''
            Temporary fix to get limits
            '''
            for sr in range(args.nSRs):
                dfs = df.loc[ ( df.mva_score < args.mvas[sr] ) & ( df.mva_score >= args.mvas[sr+1] ) & (df.event % 2 == 1) ]
                print("Adding {} events to {}".format(len(dfs),year_str+"_"+proc_tag))
                dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR'+str(sr+1)+tag, mode='a')

