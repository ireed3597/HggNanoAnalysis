import pandas
import ROOT as r
import numpy as np
from sys import exit

r.ROOT.EnableImplicitMT()

f = r.TFile.Open("hadded/run2_scan.root")
tree = f.Get("Events")

data, columns = tree.AsMatrix(return_labels=True)
df = pandas.DataFrame(data=data, columns=columns)
#print("Tree converted to a pandas.DataFrame:\n{}".format(df))
df.to_pickle('pickles/run2_scan.pkl')
