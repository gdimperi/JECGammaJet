from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.workArea = 'WORKINGAREA'
config.General.requestName = 'WORKINGDIR'
config.section_('JobType')
config.JobType.psetName = 'CMSSWCFG'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['OUTFILENAME']
config.section_('Data')
config.Data.inputDataset = 'INPUTDATASET'
config.Data.unitsPerJob = FILESPERJOB
config.Data.splitting = 'FileBased'
config.Data.publication = False
config.Data.outLFNDirBase = '/store/user/gdimperi/JetMET2015/MC/Spring15/'
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_IT_Rome'
