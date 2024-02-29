#! /usr/bin/bash

input=/home/users/iareed/HiggsDNA/Full_Samples_21Jun23/fixed_dijet_dummies/scored/
declare -A sr1_edges=(
[500]=0.9956
[550]=0.9941
[600]=0.9958
[650]=0.9962
[700]=0.9953
[750]=0.9962
[800]=0.9941
[850]=0.9960
[900]=0.9900
[950]=0.9923
[1000]=0.9904
[1100]=0.9778
[1200]=0.968
[1300]=0.963
[1400]=0.962
[1500]=0.948
)
declare -A sr2_edges=(
[500]=0.985064
[550]=0.983487
[600]=0.989298
[650]=0.983395
[700]=0.981015
[750]=0.983681
[800]=0.979005
[850]=0.984492
[900]=0.978173
[950]=0.955356
[1000]=0.96
[1100]=0.964288
[1200]=0.94971
[1300]=0.938
[1400]=0.932
[1500]=0.91
)
#1000 should be 0.226103
new=500
#for mass in 500 550 600 650 700 750 800 850 900 950 1000 1100 1200 1300 1400 1500
for mass in 1000
do
    tag="Tprime_M${mass}_04Dec23_by_hand_lower_bound"
    mva_name="Tprime_M${mass}_score"
    sr1=${sr1_edges[$mass]}
    sr2=${sr2_edges[$mass]}
    old=$new
    new=$mass
    #echo $tag
    #echo $mva_name
    #echo $sr1
    #echo $sr2
    echo "$mass & [$sr1 - 1] & [$sr2 - $sr1) \\\\"
    python mass_changer.py --old "$old" --new "$new"
    python convert_HiggsDNA_to_flashgg_fromMerged_Tprime.py --tag "$tag" --mvas $sr2 $sr1 --mva_name "$mva_name"
done

#Reset mass to 500
python mass_changer.py --old "$new" --new "500"
