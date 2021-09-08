# HggNanoAnalysis
Area for all code used in studying & developing SNT Hgg analyses, starting from nanoAOD.  
These tools are based on the C++ looper from [NanoTools](https://github.com/cmstas/NanoTools).  

## Installation and setup
**1. Clone repository**
```
git clone https://github.com/cmstas/HggNanoAnalysis -b c++_looper
cd HggNanoAnalysis
```

**2. Set-up**  
These tools include a C++ looper (test.C) that is compiled with ```compile_scripts.C```. The main code is executed in ```runAll.py```.
To run the looper:  
i) Set the NanoTools [environment](https://github.com/cmstas/NanoTools) for both python & CMSSW, by running  
  
```
# download conda installer
curl -O -L https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b 

# add conda to the end of ~/.bashrc, so relogin after executing this line
~/miniconda3/bin/conda init

# stop conda from activating the base environment on login
conda config --set auto_activate_base false
conda config --add channels conda-forge

# install package to tarball environments
conda install --name base conda-pack -y

# create environments with as much stuff from anaconda
conda create --name pyrootenv uproot pandas root matplotlib jupyter

# and then any install residual packages with pip
conda run --name pyrootenv pip install yahist mplhep
```  
  
Then, run 
```
conda activate pyrootenv
cd /cvmfs/cms.cern.ch/slc6_amd64_gcc700/cms/cmssw/CMSSW_10_2_9/ ; cmsenv ; cd -
```  

iii) compile the scripts  with  ```root -l -b -q compile_scripts.C```  and make output directory ```mkdir outputs```  
iv) run the compiled scripts with ```python runAll.py```  
