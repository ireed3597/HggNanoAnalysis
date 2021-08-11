import ROOT as r
from glob import glob
import os

check_dir = "/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016E_private_data16/"
#check_dir = "/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016F_private_data16/"
#check_dir = "/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016G_private_data16/"
#check_dir = "/hadoop/cms/store/user/legianni/skimNano-HggHtautauselection/DoubleEG_Run2016H_private_data16/"
for ifile in glob(check_dir):
    ifile = r.TFile(ifile)
    ttree = ifile.Get("Events")
    if not ttree:
        print "corrupt",ifile.GetName()
        continue
    ifile.Close()
