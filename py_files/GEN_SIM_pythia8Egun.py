# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: SinglePi0E10_pythia8_cfi --conditions auto:run1_mc -n 10 --eventcontent RAWSIM --relval 25000,100 -s GEN,SIM --datatier GEN-SIM --beamspot Realistic8TeVCollision --fileout file:step1.root
import FWCore.ParameterSet.Config as cms



process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.Geometry.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10000)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('SinglePi0E10_pythia8_cfi nevts:30000'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    # compressionAlgorithm = cms.untracked.string('LZMA'),
    # compressionLevel = cms.untracked.int32(1),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM'),
        filterName = cms.untracked.string('')
    ),
    # eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    fileName = cms.untracked.string('file:step1_EleGun_HF.root'),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    # splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v11_L1v1', '')

#Pion Generation

#process.generator = cms.EDProducer("FlatRandomPtGunProducer",
 #  PGunParameters = cms.PSet(
  #    MaxPt = cms.double(50.01),
   #  MinPt = cms.double(14.99),#
   # PartID = cms.vint32(11),
  # MaxEta = cms.double(-3),
 # MaxPhi = cms.double(3.14159265359),
# MinEta = cms.double(-5),
# MinPhi = cms.double(-3.14159265359) ## in radians

 #),
  #  Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts

   # psethack = cms.string('single electron pt 15'),
  # AddAntiParticle = cms.bool(False),
# firstRun = cms.untracked.uint32(1),
# initialSeed = cms.untracked.uint32(123456789),
# engineName = cms.untracked.string('HepJamesRandom')
#)
#process.generator = cms.EDFilter("Pythia8EGun",
 #   PGunParameters = cms.PSet(
  #      AddAntiParticle = cms.bool(False),
   #     MaxE = cms.double(50),
    #    MaxEta = cms.double(5.0),
     #   MaxPhi = cms.double(3.14159265359),
     #   MinE = cms.double(10),
      #  MinEta = cms.double(3),
       # MinPhi = cms.double(-3.14159265359),
     #  ParticleID = cms.vint32(11)
        # ParticleID = cms.vint32(111)
 #  ),
  #  PythiaParameters = cms.PSet(
   #     parameterSets = cms.vstring()
   # ),
  #  Verbosity = cms.untracked.int32(0),
   #  firstRun = cms.untracked.uint32(1),
  #  psethack = cms.string('single pi0 E 10'),
#initialSeed = cms.untracked.uint32(456213))

#process.generatorr = cms.PSet(
 #initialSeed = cms.untracked.uint32(456213),
 #  engineName = cms.untracked.string('HepJamesRandom')
#)

#W generation
process.generator = cms.EDFilter("Pythia8EGun",
	PGunParameters = cms.PSet(
		AddAntiParticle = cms.bool(False),
		MaxEta = cms.double(5),#esit to hit HF
		MaxPhi = cms.double(3.14159265359),
		MaxE = cms.double(250.0),#edit suitably
		MinEta = cms.double(3),#edit to hit HF
		MinPhi = cms.double(-3.14159265359),
		MinE = cms.double(10.0),
	ParticleID = cms.vint32(11)
		), PythiaParameters = cms.PSet(
    parameterSets = cms.vstring()
 ),
  Verbosity = cms.untracked.int32(0),
 firstRun = cms.untracked.uint32(1),
   psethack = cms.string('single pi0 E 10'),


		 initialSeed = cms.untracked.uint32(123456789),
		 engineName = cms.untracked.string('HepJamesRandom')
)



#process.generator = cms.EDFilter("Pythia8PtGun",
 #                        PGunParameters = cms.PSet(
  #      MaxPt = cms.double(100.),
   #    MinPt = cms.double(10.),
#        ParticleID = cms.vint32(11),
 #       AddAntiParticle = cms.bool(False),
  #     MaxEta = cms.double(5),
   #    MaxPhi = cms.double(3.14159265359),
    #    MinEta = cms.double(3),
     #   MinPhi = cms.double(-3.14159265359) ## in radians
      #  ),


        #                 Verbosity = cms.untracked.int32(0), ## set to 1 (or greater)  for printouts
       #                  psethack = cms.string('single electron pt 5 to 100'),
          #              firstRun = cms.untracked.uint32(1),
         #             PythiaParameters = cms.PSet(parameterSets = cms.vstring())
#initialSeed = cms.untracked.uint32(123456789),
#                 engineName = cms.untracked.string('HepJamesRandom'),


           #             )
#













#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)






# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.simulation_step = cms.Path(process.psim)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)
# filter all path with the production filter sequence
for path in process.paths:
#	getattr(process,path).insert(0, process.generator)
        getattr(process,path).insert(0, process.generator)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
