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

#### 8. Crab file for GEN_SIM step (default = pythia, change pset parameter if flatgun is needed ).Also crab files are working only in cmssw-el7 (singularity)

```bash
$ cmssw-el7
$ cmsRun  MC_gen_crab.py
```


#### 8. Crab file for ZEE_RECHIT_AOD_cfg.py is Crab_config.py

```bash
$ cmssw-el7
$ cmsRun  Crab_config.py
```

#### 9. there are 2 plugin files (one with dressed ele adn other without it) inside HFEleStudy/hfe/plugins . use the one according to your need.

