import pandas as pd
import uproot
import glob

date="15Mar2022_fixIsoTrk"
file_names = glob.glob("../outputs_UL/*"+date+"_201*.root")
file_names += glob.glob("../outputs_UL/Data_04Apr2022_fixIsoTrk*.root")
iter = 0
df = pd.DataFrame()
for file_name in file_names:
	if "vbf_" in file_name:
		continue
	print ( file_name.split("outputs")[-1].replace(".root",".pkl") )
	tree = uproot.open(file_name)["Events"]
	df_tmp = tree.pandas.df()
	if iter == 0:
		df = df_tmp
	else:
		df = pd.concat([df,df_tmp], ignore_index=True)
	iter +=1
#df.to_pickle("../pickles/run2_20UL_"+date+".pkl")
df.to_pickle("../pickles/run2_20UL_06Ape2022_fullData.pkl")
