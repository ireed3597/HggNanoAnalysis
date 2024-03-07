import json
import os
import glob
import pandas
import argparse
import root_pandas
import numpy as np
import uproot
from collections import Counter

import awkward as ak

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
              "ttHH_ggTauTau":"ttHH_ggTautau",
              "ttHH_ggWW":"ttHH_ggWW",
              "ttHH_ggbb":"ttHH_ggbb",
              "ttH_M125":"ttH",
              "HHggTauTau":"ggHH_ggTautau",
              "HHggWW_dileptonic":"ggHH_ggWWdileptonic",
              "HHggWW_semileptonic":"ggHH_ggWWsemileptonic",
              "HHggbb":"ggHH_ggbb",
              "THQ_M125":"tHq",
              "THW_M125":"tHW"
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
    #default = "/home/users/azecchin/Analysis/HiggsDNA/output/heft_presel_FF_syst_condor/"
    default = '/home/users/iareed/HiggsDNA/Full_Samples_21Jun23/fixed_dijet_dummies/'
    #default = '/ceph/cms/store/user/iareed/HiggsDNA_offload/SM_22Sep22/scored_dataframes/' 
)

parser.add_argument(
    "--tag",
    help = "unique tag to identify batch of processed samples",
    type = str,
    #default = "test"
    default = "SM_22Sep23_datacard_name_update"
)
parser.add_argument(
    "--mvas",
    nargs="*",
    help = "mva limits to SRs",
    type = float,
    default = [0.954425,0.9910]
)
parser.add_argument(
    "--mva_name",
    help = "title of mva score to use",
    type = str,
    #default = "mva_score"
    default = "mva_score"
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
print(args.mvas)
out_dir = '/home/users/iareed/CMSSW_10_2_13/src/flashggFinalFit/files_systs/' + str(args.tag) + '/'

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
    dfs = df.loc[ ( df[args.mva_name] < args.mvas[sr] ) & ( df[args.mva_name] >= args.mvas[sr+1] ) & ( ( df.Diphoton_mass < 120 ) | ( df.Diphoton_mass > 130 ) ) ]
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
                dfs = df.loc[ ( df[args.mva_name] < args.mvas[sr] ) & ( df[args.mva_name] >= args.mvas[sr+1] ) & (df.event % 2 == 1) ]
                print("Adding {} events to {}".format(len(dfs),year_str+"_"+proc_tag))
                dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125.38_13TeV_SR'+str(sr+1)+tag, mode='a')
