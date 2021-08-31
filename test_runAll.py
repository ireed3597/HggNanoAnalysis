#!/usr/bin/env python
import os
import sys
import glob
import json
import ROOT as r
from tqdm import tqdm
from ROOT import gROOT

#sys.path.append('/cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_9/python/TauAnalysis')
#sys.path.append(os.path.dirname(os.path.abspath(__file__).rsplit('/TauAnalysis/ClassicSVfit/',1)[0])+'/cfipython/slc6_amd64_gcc700/TauAnalysis/ClassicSVfit')

r.gSystem.Load('/home/users/fsetti/NanoTools/HggNanoAnalysis/TauAnalysis/ClassicSVfit/lib/libTauAnalysisClassicSVfit.so')
#r.gSystem.Load('libTauAnalysisSVfitTF.so')

