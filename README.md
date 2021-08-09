# HggNanoAnalysis
Area for all code used in studying & developing SNT Hgg analyses, starting from nanoAOD.  
These tools are based on the C++ looper from [NanoTools](https://github.com/cmstas/NanoTools).  

## Installation and setup
**1. Clone repository**
```
git clone https://github.com/cmstas/HggNanoAnalysis
cd HggNanoAnalysis
```

**2. Set-up**  
These tools include a C++ looper (test.C) that is compiled with ```compile_scripts.C```. The main code is executed in ```runAll.py```.
To run the looper:  
i) Set the NanoTools [environment](https://github.com/cmstas/NanoTools) for both python & CMSSW, by running  
```
./py_setup.sh
./setup.sh
```  
iii) compile the scripts  with  ```root -l -b -q compile_scripts.C```  
iv) run the compiled scripts with ```python runAll.py```  
