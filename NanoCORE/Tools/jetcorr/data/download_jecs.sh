#!/usr/bin/env sh

# url is like https://raw.githubusercontent.com/cms-jet/JECDatabase/master/textFiles/Autumn18_V19_MC/Autumn18_V19_MC_L1FastJet_AK4PFchs.txt

baseurl="https://raw.githubusercontent.com/cms-jet/JECDatabase/master/textFiles/"

era="Summer16_07Aug2017_V11_MC"

jettypes="
AK4PFchs
"

corrstrs="
L1FastJet
L1RC
L2L3Residual
L2Relative
L2Residual
L3Absolute
Uncertainty
"

mkdir -p $era

cd $era
for jettype in $jettypes; do
    for corrstr in $corrstrs; do
        echo Executing: curl -s -L "$baseurl/$era/${era}_${corrstr}_${jettype}.txt"
        curl -s -O -L "$baseurl/$era/${era}_${corrstr}_${jettype}.txt"
    done
done
cd -
