#### 1. Set Up the CMSSW Area

To set up the CMSSW environment, use the following commands:

```bash
$ cmsrel CMSSW_10_6_29
$ cd CMSSW_10_6_29/src
$ cmsenv
$ scram -b
```

#### 2. used CMSSW_10_6_29 with arch: slc7_amd64_gcc700. Clone the config files from github.

```bash
$ git clone git@github.com:kaizokusk/project.git
```



#### 3. cd into pyfiles folder and cmsRun config files With pythia/flatptgun and mentioned eta ranges 3-5

```bash
 $ cd project/py_files/
```



--(for pythia gun)

```bash
 $ cmsRun GEN_SIM_pythia8Egun.py
```

Or
--(for flatgun)

```bash
$ cmsRun GEN_SIM_flat_random_gun.py
```


 
#### 4. Cmsrun the 2nd config file.

```bash 
$ cmsRun genSimDigiRaw_mcProd.py
```

#### 5. cmsRun the 3rd config file

```bash
 $ cmsRun recoStepUL2018.py
```


#### 6. Then moved the output of “recoStepUL2018.py” step to test folder containing the “ZEE_RecHit_AOD_cfg.py” and then changed the input file in it to the output file name of  “recoStepUL2018.py”.

```bash
$ cp EleGun_HF_AODSIM.root ../HFEleStudy/hfe/test/
```


#### 7. Then move to folder containing Zee_rechit config file and then  cmsrun the Zee_rechit config file

```bash
$ cd ../HFEleStudy/hfe/test/
$ cmsRun  ZEE_RecHit_AOD_cfg.py
```