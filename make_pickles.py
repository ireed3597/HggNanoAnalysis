import pandas
import ROOT as r
import numpy as np
from sys import exit

r.ROOT.EnableImplicitMT()

f = r.TFile.Open("hadded/hadded_out.root")
tree = f.Get("Events")

data, columns = tree.AsMatrix(return_labels=True)

df = pandas.DataFrame(data=data, columns=columns)
df.to_pickle('pickles/inclusive.pkl')

df1 = df.loc[ df['Category'] == 1 ]
df1.to_pickle('pickles/cat1.pkl')
df2 = df.loc[ df['Category'] == 2 ]
df2.to_pickle('pickles/cat2.pkl')
df3 = df.loc[ df['Category'] == 3 ]
df3.to_pickle('pickles/cat3.pkl')
df4 = df.loc[ df['Category'] == 4 ]
df4.to_pickle('pickles/cat4.pkl')
df5 = df.loc[ df['Category'] == 5 ]
df5.to_pickle('pickles/cat5.pkl')
df6 = df.loc[ df['Category'] == 6 ]
df6.to_pickle('pickles/cat6.pkl')
df7 = df.loc[ df['Category'] == 7 ]
df7.to_pickle('pickles/cat7.pkl')
df8 = df.loc[ df['Category'] == 8 ]
df8.to_pickle('pickles/cat8.pkl')
