import uproot

file_name = "run2_26Jan2022"
tree = uproot.open("hadded/"+file_name+".root")["Events"]
tree.pandas.df().to_pickle("pickles/"+file_name+".pkl")
