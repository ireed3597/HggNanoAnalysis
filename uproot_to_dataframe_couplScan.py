import pandas as pd
import uproot
import glob

date="31Jan2022"
file_names = glob.glob("outputs/*"+date+"*.root")
#file_names += glob.glob("outputs/*"+date+"_2017.root")
#file_names += glob.glob("outputs/*"+date+"_2018.root")
iter = 0
df = pd.DataFrame()
for file_name in file_names:
	if "vbf_" in file_name:
		continue
	tree = uproot.open(file_name)["Events"]
	print file_name.split("outputs")[-1].replace(".root",".pkl")
	df_tmp = tree.pandas.df()
	if iter == 0:
		df = df_tmp
	else:
		df = pd.concat([df,df_tmp], ignore_index=True)
	iter +=1
df.to_pickle("pickles/run2_"+date+".pkl")
