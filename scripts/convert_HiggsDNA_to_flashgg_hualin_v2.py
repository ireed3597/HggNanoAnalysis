import os
import glob
import pandas
import argparse
import root_pandas
import numpy as np

years = [ '2016UL_preVFP', '2017', '2018' ]
procs_dict = {"Data":"Data",
              "DiPhoton":"DiPhoton",
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
              "ttH_M125":"ttH"
            }

#procs_dict = { "ggH_M125": "ggH", 'HHggTauTau':'HHggTauTau', 'HHggWW_dileptonic':'HHggWWdileptonic', 'HHggWW_semileptonic':'HHggWWsemileptonic', 'ttH_M125':'ttH', 'VBFH_M125':'VBFH', 'VH_M125':'VH', "data":"Data" }

parser = argparse.ArgumentParser()

parser.add_argument(
    "--input",
    help = "path to input parquet directory",
    type = str,
    default = "/home/users/smay/HiggsDNA/tthh_sr_20Apr2022"
    #default = "/ceph/cms/store/user/fsetti/HiggsDNA_output/09Apr2022_sr_resonant_systs/"
)

parser.add_argument(
    "--tag",
    help = "unique tag to identify batch of processed samples",
    type = str,
    default = "debug"
    #default = "10Apr2022"
)
parser.add_argument(
    "--mvas",
        nargs='*',
    help = "mva limits to SRs",
    type = float,
    # low cut, high cut
    default = [0.969307, 0.9883]
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

#out_dir = '/home/users/fsetti/ic_flashgg/CMSSW_10_2_13/src/flashggFinalFit/files_systs/' + str(args.tag) + '/'
out_dir = '/home/users/iareed/ttHHggbb/coupling_scan/CMSSW_10_2_13/src/flashggFinalFit/files_systs/' + str(args.tag) + '/'

os.system("rm -rf %s"%(out_dir))
os.system("mkdir -p %s"%(out_dir))

os.system("mkdir -p %s/Data"%(out_dir))
os.system("mkdir -p %s/2016"%(out_dir))
#os.system("mkdir -p %s/2016_APV"%(out_dir))
os.system("mkdir -p %s/2017"%(out_dir))
os.system("mkdir -p %s/2018"%(out_dir))

procs = glob.glob(str(args.input)+'/*')
doneData=False

for proc in procs[:]:

    if ("Data" in proc) and not doneData:
        files = glob.glob(str(args.input)+'/Data*/*.parquet')
    #Process Data
        for file_ in files:
            df = pandas.read_parquet(file_, engine='pyarrow')
            df['CMS_hgg_mass'] = df['Diphoton_mass']
            for sr in range(args.nSRs):
                dfs = df.loc[ (df.process_id == 0 ) & ( df.bdt_score < args.mvas[sr] ) & ( df.bdt_score >= args.mvas[sr+1] ) & ( ( df.Diphoton_mass < 120 ) | ( df.Diphoton_mass > 130 ) ) ]
                dfs.to_root(out_dir+'/Data/'+'/allData.root',key='Data_13TeV_SR'+str(sr+1), mode='a')
        doneData=True
        continue


    for year in years:
        if year not in proc.split("/")[-1]:
            continue

        #get all files including systematic variations
        files = glob.glob(proc+'/*.parquet')

        for file_ in files:
            print ("Now processing: ", file_)
            df = pandas.read_parquet(file_, engine='pyarrow')
            if year == '2016UL_preVFP':
                df_ext1                                                 = pandas.read_parquet(file_.replace("2016UL_preVFP","2016UL_postVFP"), engine='pyarrow')
                df = pandas.concat([ df, df_ext1 ], ignore_index=True)

            tag     =       ''
            if 'nominal' not in file_.split("/")[-1]:
                if 'up' in file_.split("/")[-1]:
                    tag     = file_.split("merged")[-1]
                    tag     = tag.split("_up")[0]
                    tag     += 'Up01sigma'
                if 'down' in file_.split("/")[-1]:
                    tag     = file_.split("merged")[-1]
                    tag     = tag.split("_down")[0]
                    tag     += 'Down01sigma'
            if 'scale' in file_.split("/")[-1]:
                tag = '_MCScale' + tag

            #Define hgg_mass & dZ variable
            df['CMS_hgg_mass'] = df['Diphoton_mass']
            df['dZ']        = np.ones(len(df['Diphoton_mass']))
            df['weight'] = df['weight_central']
            yield_systematics       = [ key for key in df.keys() if ( "weight_" in key ) and ( "_up" in key or "_down" in key )]
            rename_sys      = {}
            for sys in yield_systematics:
                #a bit of gymnastics to get the inputs right for Mr. flashggFinalFit
                sys_central = sys.replace("_up","_central")
                sys_central = sys.replace("_down","_central")
                if 'btag_deepjet' in sys_central:
                    sys_central = sys_central.replace('_up','_central').split('_central')[0]+'_central'
                    sys_central = sys_central.replace('_down','_central').split('_central')[0]+'_central'
                print(proc, sys, df[sys].mean(), df[sys_central].mean(), (df[sys]/df[sys_central]).mean())
                df[sys]                 = df[sys] / df[sys_central]
                if "_up" in sys:
                    rename_sys[sys] = sys.replace("_up","Up01sigma")
                if "_down" in sys:
                    rename_sys[sys] = sys.replace("_down","Down01sigma")
            #print(rename_sys)
            df = df.rename(columns=rename_sys)

            #Process Signal
            year_str = year
            if '2016' in year:
                year_str = '2016'

            proc_tag = ''
            if 'EFT' in proc.split("/")[-1]:
                proc_tag = 'ggHHbm' + proc.split("/")[-1].split("node_")[-1].split("_201")[0]
            else:
                proc_tag = procs_dict[proc.split("/")[-1].split("_201")[0]]

            '''
            Temporary fix to get limits
            '''
            for sr in range(args.nSRs):
                dfs = df.loc[ ( df.bdt_score < args.mvas[sr] ) & ( df.bdt_score >= args.mvas[sr+1] ) & ( df.weight_central > -9999 ) ]
                dfs.to_root(out_dir+year_str+'/'+proc_tag+'_125_13TeV.root',''+proc_tag+'_125_13TeV_SR'+str(sr+1)+tag, mode='a')
