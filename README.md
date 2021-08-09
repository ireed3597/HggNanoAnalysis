# HggNanoAnalysis
Area for all code used in studying & developing SNT Hgg analyses, starting from nanoAOD.

## Installation and setup
**1. Clone repository**
```
git clone https://github.com/cmstas/HggNanoAnalysis
cd HggNanoAnalysis
```

**2. Set-up**
These tools include a C++ looper (test.C) that is compiled with ```compile_scripts.C```. The main code is executed in ```runAll.py```.
To run the looper:  
i) Set the NanoTools [environment](https://github.com/cmstas/NanoTools) for both python & CMSSW, in that order  
ii) Load the ad-hoc library from this NanoCORE repository (it contains a few additions to the official NanoTools one)  
iii) compile the scripts  
iv) run the compiled scripts with ```python runAll.py```  
