import CRABClient
from CRABClient.UserUtilities import config 
config = config()

config.General.requestName = 'Gen_pythiagun'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'GEN_SIM_pythia8Egun.py'
config.JobType.numCores=1
config.Data.outputPrimaryDataset = 'Gen_SIM'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1
NJOBS = 1000 # This is not a configuration parameter, but an auxiliary variable that we use in the next line.
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True
config.Data.outputDatasetTag = 'CRAB3_GEN_SIM_PYTHIA'
config.Data.outLFNDirBase = '/store/user/spodem/test_crab'
config.Site.storageSite = 'T3_CH_CERNBOX'

