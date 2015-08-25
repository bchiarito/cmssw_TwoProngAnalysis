#!/usr/bin/env python

#
# Support functions for limit setting
#
# Seth I. Cooper, U. Alabama
# November 21 2012
# update to combine and add multiple analysis channels
#    July/August 2015

import string
import os
import sys
import glob
import math
import array
import subprocess
import itertools
import tempfile
import multiprocessing
import copy
import glob

# run root in batch
sys.argv.append('-b')
from ROOT import *
from ModelPoint import *
from MakePlots import *

# systematics for all model points
SigPUSyst = 0.007
SigPDFSyst = 0.05
#SigScaleFactorSystOneGamma = 0.005
SigScaleFactorSystOneGamma = 0.0014
SigPtSFSystOneGamma = 0.03
# background --> take from hists in root files

# globally defined hists
signalHist = TH1F()
signalAccOnlyHistogram = TH1F()
signalTotalEventsHistogram = TH1F()
signalEntriesTotal = 0
signalEntriesTotalUnscaled = 0
signalHistSmeared = TH1F()
signalHistScaleShiftUp = TH1F()
signalHistScaleShiftDown = TH1F()
signalHistPileupShiftUp = TH1F()
signalHistPileupShiftDown = TH1F()
backgroundHist = TH1F()
dataHist = TH1F()
histogramSignalFile = 0
histogramBGFile = 0
histogramDataFile = 0

Farm_Directories = []
optResultsList = []
limitSettingCalled = False
optimizationCouplingsAlreadyDone = []

# rebin the hists to get rid of background=0 at high mass
newBins = [i*10 for i in range(0,150)]
newBins.extend([1500+i*20 for i in range(0,200/20)])
newBins.extend([1700+i*50 for i in range(0,700/50)])
newBins.extend([2500+i*100 for i in range(0,500/100)])
newBins.extend([3000+i*300 for i in range(0,500/300+1)])
newBins.append(4000)
newBinsArr = array.array('d',newBins)


def GetListOfChannels(mpArray):
  channels = []
  for mp in mpArray:
    if not mp.channel in channels:
      channels.append(mp.channel)
  return channels


def CreateFarmDirectoryStructure(FarmDirectory):
  global Farm_Directories
  Farm_Directories = [FarmDirectory+'/',
                        FarmDirectory+'/inputs/',
                        FarmDirectory+'/outputs/',
                        FarmDirectory+'/logs/',
                        FarmDirectory+'/errors/']
  for i in range(0,len(Farm_Directories)):
    if os.path.isdir(Farm_Directories[i]) == False:
      os.system('mkdir -p ' + Farm_Directories[i])


def CreateTheCmdFile(PathCmd):
  cmd_file=open(PathCmd,'w')
  cmd_file.write('#!/bin/bash\n')
  cmd_file.write('# list all bsub commands' + '\n')
  cmd_file.close()


def AddJobToCmdFile(PathCmd,PathShell,PathLog,PathError,doBatch,QueueName):
  cmd_file=open(PathCmd,'a')
  cmd_file.write('\n')
  if doBatch:
    cmd_file.write('bsub -q "%s" '     % QueueName)
    cmd_file.write(' -o %s' % PathLog)
    cmd_file.write(' -e %s ' % PathError)
  cmd_file.write(PathShell)
  if not doBatch:
    cmd_file.write(' >& %s' % PathLog)
  cmd_file.write('\n')
  cmd_file.close()


def CreateShellFileRooStats(modelPoint,lumi,lumiErr,limitsDir):
  global Farm_Directories
  couplingStr = str(modelPoint.coupling).replace('.','p')
  massString = str(int(modelPoint.mass))
  plotFileSamplingName = 'cls_sampling_k_'+couplingStr+"_m"+massString+".pdf"
  plotFileName = 'cls_k_'+couplingStr+"_m"+massString+".pdf"
  outputFile = Farm_Directories[2]+'limits_k_'+couplingStr+'_m'+massString+'.txt'
  if not os.path.isfile(outputFile):
    os.system('touch '+outputFile)
  PathShell = Farm_Directories[1]+'limits_k_'+couplingStr+'_m'+massString+'.sh'
  PathLog = Farm_Directories[3]+'limits_k_'+couplingStr+'_m'+massString+'.log'
  PathError = Farm_Directories[4]+'limits_k_'+couplingStr+'_m'+massString+'.err'
  shell_file=open(PathShell,'w')
  shell_file.write('#!/bin/sh\n')
  shell_file.write('source /afs/cern.ch/cms/cmsset_default.sh\n')
  shell_file.write('cd ' + Farm_Directories[2] + '\n')
  shell_file.write('eval `scram run -sh`\n')
  shell_file.write('cd -' + '\n')
  #void ComputeLimit(float lumi, float lumiError, float totalEff, float totalEffErrStat, float totalEffErrSyst,
  #float nBackground, float nBackgroundErrStat, float nBackgroundErrSyst,
  #float nDataObs, 
  #TODO not needed for computelimit -- save/remove
  #float totalEffMScaleSystUp, float totalEffMScaleSystDown,
  #float totalEffMResSystUp, float totalEffMResSystDown,
  #float totalEffPileupSystUp, float totalEffPileupSystDown,
  #int mass, float coupling, float halfWidth,
  #float totalXSec, int massWindowLow, int massWindowHigh, std::string fileName)
  command = os.getcwd()+'/'
  command+= 'ComputeLimit.C('+str(lumi)+','+str(lumiErr)+','+str(modelPoint.totalEff)+','
  command+=str(modelPoint.totalEffErrStat)+','+str(modelPoint.totalEffErrSyst)+','
  command+=str(modelPoint.totalEffMScaleSystUp)+','+str(modelPoint.totalEffMScaleSystDown)+','
  command+=str(modelPoint.totalEffMResSystUp)+','+str(modelPoint.totalEffMResSystDown)+','
  command+=str(modelPoint.totalEffPileupSystUp)+','+str(modelPoint.totalEffPileupSystDown)+','
  command+=str(modelPoint.nBackground)+','+str(modelPoint.nBackgroundErrStat)+','+str(modelPoint.nBackgroundErrSyst)+','
  command+=str(modelPoint.nDataObs)+','+str(modelPoint.mass)+','
  command+=str(modelPoint.coupling)+','+str(modelPoint.halfWidth)+','+str(modelPoint.totalXSec)+','
  command+=str(modelPoint.optMassWindowLow)+','+str(modelPoint.optMassWindowHigh)+','
  command+='\\"'+outputFile+'\\"'
  command+=')'
  shell_file.write('root -l -b -q "'+command+'"\n')
  shell_file.write('mv cls_sampling.pdf '+limitsDir+'/cls_plots/'+plotFileSamplingName+'\n')
  shell_file.write('mv cls.pdf '+limitsDir+'/cls_plots/'+plotFileName+'\n')
  shell_file.close()
  os.system('chmod 777 '+PathShell)
  return PathShell,PathLog,PathError


# takes an array of model points with the same mass/coupling but different channels
def CreateShellFileCombine(modelPoints):
  global Farm_Directories
  couplingStr = str(modelPoints[0].coupling).replace('.','p')
  massString = str(int(modelPoints[0].mass))
  #ch = modelPoint.channel
  ch = "_".join([mp.channel for mp in modelPoints])
  # create datacard
  pathCard = CreateCombineDatacard(modelPoints)
  #plotFileSamplingName = 'cls_sampling_k_'+couplingStr+"_m"+massString+".pdf"
  #plotFileName = 'cls_k_'+couplingStr+"_m"+massString+".pdf"
  #outputFile = Farm_Directories[2]+'limits_k_'+couplingStr+'_m'+massString+'_'+ch+'.txt'
  #if not os.path.isfile(outputFile):
  #  os.system('touch '+outputFile)
  PathShell = Farm_Directories[1]+'limits_k_'+couplingStr+'_m'+massString+'_'+ch+'.sh'
  PathLog = Farm_Directories[3]+'limits_k_'+couplingStr+'_m'+massString+'_'+ch+'.log'
  PathError = Farm_Directories[4]+'limits_k_'+couplingStr+'_m'+massString+'_'+ch+'.err'
  shell_file=open(PathShell,'w')
  shell_file.write('#!/bin/sh\n')
  shell_file.write('source /afs/cern.ch/cms/cmsset_default.sh\n')
  shell_file.write('cd ' + Farm_Directories[2] + '\n')
  shell_file.write('eval `scram run -sh`\n')
  shell_file.write('cd -' + '\n')
  # -n _k0p01_m750_BB  --> 1st root file will be called higgsCombine_k0p01_m750_BB.Asymptotic.mH120.root
  # -n _k0p01_m750_BB  --> 2nd root file will be called higgsCombine_k0p01_m750_BB.HybridNew.mH120.SEED.root where SEED is the random seed used
  nameString = '_k'+couplingStr+'_m'+massString+'_'+ch
  # check error code and redo again "automatically"
  shell_file.write('combine -n '+nameString+' -M Asymptotic --datacard '+pathCard+'\n')
  shell_file.write('if [ $? -ne 0 ]; then\n')
  shell_file.write('  rm higgsCombine'+nameString+'.Asymptotic.mH120.root\n')
  shell_file.write('  combine -n '+nameString+' -M Asymptotic --datacard '+pathCard+'\n')
  shell_file.write('fi\n')
  shell_file.write('mv higgsCombine'+nameString+'.Asymptotic.mH120.root '+Farm_Directories[2]+'\n')
  shell_file.write('combine -n '+nameString+' -H Asymptotic -M HybridNew --rule CLs --frequentist -s -1 --saveToys --saveHybridResult --datacard '+pathCard+'\n')
  shell_file.write('if [ $? -ne 0 ]; then\n')
  shell_file.write('  rm higgsCombine'+nameString+'.HybridNew.mH120.*.root\n')
  shell_file.write('  combine -n '+nameString+' -H Asymptotic -M HybridNew --rule CLs --frequentist -s -1 --saveToys --saveHybridResult --datacard '+pathCard+'\n')
  shell_file.write('fi\n')
  shell_file.write('mv higgsCombine'+nameString+'.HybridNew.mH120.*.root '+Farm_Directories[2]+'\n')
  # clean up combine's mess
  shell_file.write('rm roostats-*' + '\n')
  shell_file.close()
  os.system('chmod 777 '+PathShell)
  return PathShell,PathLog,PathError
  

def RunSubprocessJob(cmd,logFile):
  print 'Working in process #%d' % os.getpid()
  print '----------> cmd=',cmd
  print '----------> logFile=',logFile
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
  out,err = p.communicate()
  #print 'write output to logFile='+logFile
  #print '----------> out=',out
  with open(logFile,'w') as file:
    file.write(out) # only have stdout since we redirected stderr to stdout
  #return out,err
  

def SubmitJobs(pathCmd, doBatch):
  global limitSettingCalled
  if limitSettingCalled:
    return # already set limits (prevent recursive submission)
  if doBatch:
    os.system('source '+pathCmd)
  else:
    # multiprocess
    pool = multiprocessing.Pool(multiprocessing.cpu_count()) # one job per CPU core (although OS will do the actual schedule on cores)
    #results = []
    #logFiles = []
    with open(pathCmd, 'r') as file:
      for line in file:
        li = line.strip()
        if not li.startswith('#') and len(li):
          #print line
          cmd = line.rstrip().split()[0]
          log = line.rstrip().split()[-1]
          if len(cmd) and len(log):
            #print 'command='+cmd
            #print 'log='+log
            #logFiles.append(log)
            #results.append(pool.apply_async(RunJob,cmd,log))
            pool.apply_async(RunSubprocessJob,args=(cmd,log))
    pool.close()
    pool.join() # wait for completion
    #for res in results:
    #  out,err = result.get()
  limitSettingCalled = True


# takes an array of model points with the same mass/coupling but different channels
def CreateCombineDatacard(modelPoints):
  global Farm_Directories
  global Lumi
  global LumiErr
  couplingStr = str(modelPoints[0].coupling).replace('.','p')
  massStr = str(int(modelPoints[0].mass))
  ch = "_".join([mp.channel for mp in modelPoints])
  numChan = len(modelPoints)
  numBacks = 1

  PathCard = Farm_Directories[1]+'datacard_k'+couplingStr+'_m'+massStr+'_'+ch+'.txt'
  sys.stdout.write('\rwriting card                                                                                                                                             ')
  sys.stdout.write('\rwriting card '+PathCard)
  sys.stdout.flush()
  card_file=open(PathCard,'w')
  card_file.write('#datacard_k'+couplingStr+'_m'+massStr+'_'+ch+'\n')
  card_file.write('imax '+str(numChan)+' number of channels'+'\n')
  card_file.write('jmax '+str(numBacks)+' number of backgrounds'+'\n')
  card_file.write('kmax 3 number of nuisance parameters (sources of systematic uncertainties)'+'\n')
  card_file.write('------------'+'\n')
  #card_file.write('# we have just one channel, in which we observe '+dataObsStr+' events'+'\n')
  card_file.write('bin    ')
  for mp in modelPoints:
    card_file.write(mp.channel+'    ')
  card_file.write('\n')
  card_file.write('observation ')
  for mp in modelPoints:
    card_file.write(str(int(mp.nDataObs))+'  ')
  card_file.write('\n')
  card_file.write('------------'+'\n')
  card_file.write('# now we list the expected events for signal and all backgrounds in that bin'+'\n')
  card_file.write('# the second \'process\' line must have a positive number for backgrounds, and 0 for signal'+'\n')
  card_file.write('# then we list the independent sources of uncertainties, and give their effect (syst. error)'+'\n')
  card_file.write('# on each process and bin'+'\n')
  #card_file.write('bin              1     1'+'\n')
  #card_file.write('process         RSG   Bckg'+'\n')
  #card_file.write('process          0     1'+'\n')
  #card_file.write('rate           '+expSigEvtsStr+'  '+nBackStr+'\n')
  #card_file.write('#expSigEvts = effAcc*thxsec*lumi = '+str(modelPoint.totalEff)+'*'+str(modelPoint.totalXSec)+'*'+str(Lumi)+' = '+expSigEvtsStr+'\n')
  card_file.write('bin    ')
  for mp in modelPoints:
    card_file.write(mp.channel+'    ')
    for i in range(0,numBacks):
      card_file.write(mp.channel+'    ')
  card_file.write('\n')
  card_file.write('process     ')
  for mp in modelPoints:
    card_file.write('RSG    ')
    for i in range(0,numBacks):
      card_file.write('Bckg    ')
  card_file.write('\n')
  card_file.write('process     ')
  for mp in modelPoints:
    card_file.write('0    ')
    for i in range(0,numBacks):
      card_file.write(str(i+1)+'    ')
  card_file.write('\n')
  card_file.write('rate    ')
  for mp in modelPoints:
    expSigEvts = mp.totalEff*mp.totalXSec*Lumi
    expSigEvtsStr = str(expSigEvts)
    card_file.write(expSigEvtsStr+'    ')
    for i in range(0,numBacks):
      nBackStr = str(mp.nBackground)
      card_file.write(nBackStr+'    ')
  card_file.write('\n')
  card_file.write('------------'+'\n')
  #card_file.write('lumi     lnN    '+str(LumiErr)+'   '+str(LumiErr)+'  lumi for both'+'\n')
  #card_file.write('eff_rsg  lnN    '+effAccErrStr+'     -      signal efficiency'+'\n')
  #card_file.write('bg       lnN      -       '+nBackErrStr+'  uncertainty on backgrounds'+'\n')
  card_file.write('lumi    lnN    ')
  for mp in modelPoints:
    card_file.write(str(LumiErr)+'    ')
    for i in range(0,numBacks):
      card_file.write(str(LumiErr)+'    ')
  card_file.write('lumi for all\n')
  card_file.write('eff_rsg    lnN    ')
  for mp in modelPoints:
    effAccErr = 1.0+math.sqrt(math.pow(mp.totalEffErrStat,2)+math.pow(mp.totalEffErrSyst,2))/mp.totalEff
    card_file.write(str(effAccErr)+'    ')
    for i in range(0,numBacks):
      card_file.write('  -  ')
  card_file.write('signal efficiency\n')
  card_file.write('bkg    lnN    ')
  for mp in modelPoints:
    card_file.write('  -  ')
    for i in range(0,numBacks):
      nBackErr = 1.0+math.sqrt(math.pow(mp.nBackgroundErrStat,2)+math.pow(mp.nBackgroundErrSyst,2))/mp.nBackground
      card_file.write(str(nBackErr)+'    ')
  card_file.write('uncertainty on backgrounds\n')
  card_file.close()
  return PathCard


# returns list of list of modelpoints in all different channel combinations
# example:
#  [ (BB), (BE), (EB), (EE), (BB,BE), ..., (BB,BE,EB), ..., (BB,BE,EB,EE) ]
def GetAllChanCombs(modelPointArray,verbose=0):
  allMpsAllChanCombs = []
  #channels = []
  modelPointsCopy = copy.deepcopy(modelPointArray)
  while len(modelPointsCopy) > 0:
    mp = modelPointsCopy[0]
    coupling = mp.coupling
    channel = mp.channel
    #if not channel in channels:
    #  channels.append(channel)
    mass = mp.mass
    thisMpAllChannels = [modelp for modelp in modelPointArray if modelp.coupling==coupling and modelp.mass==mass and modelp.channel != channel]
    thisMpAllChannels.append(mp)
    for mp in thisMpAllChannels:
      #print 'LOOK FOR MP'
      #mp.Print()
      for mp2 in modelPointsCopy:
        if mp2.coupling==mp.coupling and mp2.mass==mp.mass and mp2.channel==mp.channel:
          modelPointsCopy.remove(mp2)
    thisMpAllChanCombs = []
    for L in range(1, len(thisMpAllChannels)+1):
      els = [list(x) for x in itertools.combinations(thisMpAllChannels, L)]
      thisMpAllChanCombs.extend(els)
    allMpsAllChanCombs.extend(thisMpAllChanCombs)
    if verbose > 1:
      for x in thisMpAllChanCombs:
        print '[',
        for mp in x:
          #mp.Print()
          print 'k_'+str(mp.coupling)+'_m'+str(int(mp.mass))+'_'+mp.channel+',',
        print ']'
  if verbose > 2:
    print '[ ',
    for mpList in allMpsAllChanCombs:
        print '[',
        for mp in mpList:
          #mp.Print()
          print 'k_'+str(mp.coupling)+'_m'+str(int(mp.mass))+'_'+mp.channel+',',
        print '] ',
    print ' ]'
  return allMpsAllChanCombs


def ComputeLimits(cl95MacroPath,doBatch,queueName,lumi,lumiErr,modelPointArray,limitsDir,useCombine):
  global Farm_Directories
  global Lumi
  global LumiErr
  Lumi = lumi
  # lumiErr here is the absolute value, error*lumi. combine needs relative. so divide this by lumi and add 1
  LumiErr = 1+lumiErr/lumi
  if not os.path.isdir(limitsDir+'/cls_plots'):
    os.mkdir(limitsDir+'/cls_plots')
  farmDirBase = limitsDir+'/Farm'
  CreateFarmDirectoryStructure(farmDirBase)
  cmdPath = farmDirBase+'/submit.sh'
  CreateTheCmdFile(cmdPath)
  ## find all channels
  #channels = []
  #for mp in modelPointArray:
  #  if not mp.channel in channels:
  #    channels.append(mp.channel)
  ## make all combinations of channels from the channel list
  #chanCombs = []
  #for L in range(1, len(stuff)+1):
  #  els = [list(x) for x in itertools.combinations(stuff, L)]
  #  chanCombs.extend(els)
  ###
  ###str(decimal.Decimal('0.1'))
  ### now make lists of modelPoints with the same but different channels
  #modelPointsByChannel = {}
  #for ch in channels:
  #  modelPointsByChannel[ch] = []
  #for mp in modelPointArray:
  #  modelPointsByChannel[mp.channel].append(mp)
  ###
  #modelPointsCopy = copy.deepcopy(modelPointArray)
  #while len(modelPointsCopy) > 0:
  #  mp = modelPointsCopy.pop(0)
  #  thisMpAllChannels = []
  #  thisMpAllChannels.append(mp)
  if not useCombine:
    for modelPoint in modelPointArray:
      pathShell,pathLog,pathError = CreateShellFileRooStats(modelPoint,lumi,lumiErr,limitsDir)
      AddJobToCmdFile(cmdPath,pathShell,pathLog,pathError,doBatch,queueName)
  else:
    # must pass CreateShellFileCombine all combinations of channels for which to calculate the limits
    allMpsAllChanCombs = GetAllChanCombs(modelPointArray)
    for modelPointList in allMpsAllChanCombs:
      pathShell,pathLog,pathError = CreateShellFileCombine(modelPointList)
      AddJobToCmdFile(cmdPath,pathShell,pathLog,pathError,doBatch,queueName)
  print
  print 'Submit all jobs'
  SubmitJobs(cmdPath,doBatch)
  #FIXME make git-compatible somehow
  #print 'Used StatisticalTools/RooStatsRoutine CVS tag:',GetRooStatsMacroCVSTag(cl95MacroPath)


def MergeLimitJobs(modelPointArray,limitsDir,limitsFileNameBase, useCombine):
  if not useCombine:
    MergeLimitsRooStats(modelPointArray,limitsDir,limitsFileNameBase)
    return True
  else:
    retVal = MergeLimitsCombine(modelPointArray,limitsDir,limitsFileNameBase)
    return retVal


# FIXME update for channel support
def MergeLimitsRooStats(modelPointsArray,limitsDir,limitsFileNameBase):
  global Farm_Directories
  farmDirBase = limitsDir+'/Farm'
  CreateFarmDirectoryStructure(farmDirBase)
  newModelPointArray = []
  for modelPoint in modelPointArray:
    # FIXME save original info; eventually do all non-limit-related stuff like this
    # TODO removing from pass-through in ComputeLimit.C macro
    thisMPFileName = modelPoint.fileName
    thisMPKFactor = modelPoint.kFactor
    thisMPAcceptance = modelPoint.acceptance
    thisMPpreMWEff = modelPoint.preMWEff
    thisMPsignalEvents = modelPoint.totalSignalEvents
    # read this model point (now with limits) from output file
    couplingStr = str(modelPoint.coupling).replace('.','p')
    massString = str(int(modelPoint.mass))
    outputFile = Farm_Directories[2]+'limits_k_'+couplingStr+'_m'+massString+'.txt'
    with open(outputFile, 'r') as file:
      lines = file.readlines()
    mpFileArr = ReadFromLines(lines)
    if len(mpFileArr) < 1:
      print 'ERROR: Did not find model point in file:',outputFile
      print 'Perhaps the limit-setting job failed? Quitting.'
      exit(-1)
    thisModelPoint = mpFileArr[0]
    if thisModelPoint.expLimit==-1:
      print 'ERROR: Did not find exp. limit result for k=',modelPoint.coupling,'mass=',str(modelPoint.mass),'in file:',outputFile
      print 'Perhaps the limit-setting job failed? Quitting.'
      exit(-1)
    thisModelPoint.fileName = thisMPFileName
    thisModelPoint.kFactor = thisMPKFactor
    thisModelPoint.acceptance = thisMPAcceptance
    thisModelPoint.preMWEff = thisMPpreMWEff
    thisModelPoint.totalSignalEvents = thisMPsignalEvents
    newModelPointArray.append(thisModelPoint)
    ## replace this model point in the array
    #for index,mp in enumerate(modelPointArray):
    #  if mp.mass==thisModelPoint.mass and mp.coupling==thisModelPoint.coupling:
    #    modelPointArray[index] = thisModelPoint
    ## write out the result file
    #if os.path.isfile(fileName):
    #  os.remove(fileName)
  for mp in newModelPointArray:
    couplingString = str(mp.coupling).replace('.','p')
    fileName=limitsFileNameBase+couplingString+'.txt'
    with open(fileName,'a') as file:
      mp.Write(file)


def GetTFileContents(f):
  returnStr = ''
  f.cd()
  keyInfoTuple = tuple((key.GetClassName(),key.GetName()) for key in gDirectory.GetListOfKeys())
  for className,name in keyInfoTuple:
    returnStr+=className+'          '+name+'\n'
  return returnStr


def MergeLimitsCombine(modelPointArray,limitsDir,limitsFileNameBase):
  # the input model points have been read from the optimization txt file
  global Farm_Directories
  farmDirBase = limitsDir+'/Farm'
  CreateFarmDirectoryStructure(farmDirBase)
  newModelPointArray = []
  errorMsg = ''
  jobsToRerun = ''
  filesToRemove = ''
  allMpsAllChanCombs = GetAllChanCombs(modelPointArray)
  for modelPointList in allMpsAllChanCombs:
    #print 'modelPointList: channels =',
    #print [mp.channel for mp in modelPointList],
    #print ' objects:',
    #print [mp for mp in modelPointList]
    if len(modelPointList)==1:
      mp = copy.deepcopy(modelPointList[0])
    else:
      # don't write the efficiencies, etc. for multi-channel points
      mp = ModelPoint(coupling=modelPointList[0].coupling,mass=modelPointList[0].mass,totalXSec=modelPointList[0].totalXSec)
    #print 'modifying',mp
    couplingStr = str(mp.coupling).replace('.','p')
    massStr = str(int(mp.mass))
    chan = "_".join([mp2.channel for mp2 in modelPointList])
    mp.channel = chan
    #print 'mp:',mp,'channel is now:',mp.channel
    shellFileName = Farm_Directories[1]+'limits_k_'+couplingStr+'_m'+massStr+'_'+chan+'.sh'
    # read the limit results and update the mp
    # 1st root file will be called, e.g., higgsCombine_k0p01_m750_BB.Asymptotic.mH120.root
    # 2nd root file will be called, e.g., higgsCombine_k0p01_m750_BB.HybridNew.mH120.SEED.root where SEED is the random seed used
    fileNameAsymp = Farm_Directories[2]+'higgsCombine_k'+couplingStr+'_m'+massStr+'_'+chan+'.Asymptotic.mH120.root'
    fileAsymp = TFile(fileNameAsymp)
    #tree = fileAsymp.Get('limit')
    tree = MakeNullPointer(TTree)
    try:
      fileAsymp.GetObject('limit',tree)
    except:
      errorMsg+='No TTree named "limit" in file: '+fileNameAsymp+'\n'
      errorMsg+='TFileContents:\n'
      errorMsg+='\t'+GetTFileContents(fileAsymp)
      errorMsg+='Please re-run the limit job: '+shellFileName+'\n\n'
      jobsToRerun+='\t'+shellFileName+'\n'
      filesToRemove+='\trm '+fileNameAsymp+'\n'
      continue
    tol = 1e-3 # tolerance for float subtraction/comparison
    # limits in tree are the ratio of cross sections (limit/exp) so convert to xsec limits
    treeContents = []
    for lim in tree:
      if math.fabs(lim.quantileExpected-0.025) < tol:
        mp.expLimitTwoSigmaLow = lim.limit*mp.totalXSec
      elif math.fabs(lim.quantileExpected-0.16) < tol:
        mp.expLimitOneSigmaLow = lim.limit*mp.totalXSec
      elif math.fabs(lim.quantileExpected-0.5) < tol:
        mp.expLimit = lim.limit*mp.totalXSec
      elif math.fabs(lim.quantileExpected-0.84) < tol:
        mp.expLimitOneSigmaHigh = lim.limit*mp.totalXSec
      elif math.fabs(lim.quantileExpected-0.975) < tol:
        mp.expLimitTwoSigmaHigh = lim.limit*mp.totalXSec
      treeContents.append([str(lim.quantileExpected),str(lim.limit)])
    fileAsymp.Close()
    # check extracted limits
    if mp.expLimitTwoSigmaLow < 0 or mp.expLimitTwoSigmaHigh < 0 or mp.expLimitOneSigmaLow < 0 or mp.expLimitOneSigmaHigh < 0 or mp.expLimit < 0:
      errorMsg+='Problem in limit tree--one of the obtained expected limits was negative from tree in file: '+fileNameAsymp+'\n'
      errorMsg+='TTree limit contents (should have 6 expLim rows + 1 obsLim row):\n'
      errorMsg+='quantileExp\tlimit\n'
      for listPair in treeContents:
        errorMsg+='\t'+listPair[0]+'\t'+listPair[1]+'\n'
      if len(treeContents) < 1: errorMsg+='\t[empty]\n'
      errorMsg+='Please re-run the limit job: '+shellFileName+'\n\n'
      jobsToRerun+='\t'+shellFileName+'\n'
      filesToRemove+='\trm '+fileNameAsymp+'\n'
      continue
    # now for the observed limit
    # find the root file
    fileNameHNglob = Farm_Directories[2]+'higgsCombine_k'+couplingStr+'_m'+massStr+'_'+chan+'.HybridNew.mH120.*.root'
    fileNamesHN = glob.glob(fileNameHNglob)
    if len(fileNamesHN) > 1 or len(fileNamesHN) < 1:
      errorMsg+='Multiple (or zero) output files for this modelpoint found; list of filenames for glob '+fileNameHNglob+': \n'+'\n'.join(fileNamesHN)
      if len(fileNamesHN) > 1:
        errorMsg+='\nNot sure why this happened. Try removing the files and re-running the limit job: '+shellFileName+'\n\n'
        jobsToRerun+='\t'+shellFileName+'\n'
        filesToRemove+='\trm '+'\n\trm '.join(fileNamesHN)+'\n'
      else:
        errorMsg+='\n\tPlease re-run the limit job: '+shellFileName+'\n\n'
        jobsToRerun+='\t'+shellFileName+'\n'
      continue
    fileNameHN = fileNamesHN[0]
    #fileNameHN = max(glob.iglob(fileNameHNglob), key=os.path.getctime)
    ##old issue: we have multiple outputs files now--why? just take newest match for now
    fileHN = TFile(fileNameHN)
    tree = MakeNullPointer(TTree)
    try:
      fileHN.GetObject('limit',tree)
    except:
      errorMsg+='No TTree named "limit" in file: '+fileNameHN+'\n'
      errorMsg+='TFileContents:\n'
      errorMsg+='\t'+GetTFileContents(fileHN)
      errorMsg+='Please re-run the limit job: '+shellFileName+'\n\n'
      jobsToRerun+='\t'+shellFileName+'\n'
      filesToRemove+='\trm '+fileNameHN+'\n'
      continue
    treeContents = []
    for lim in tree:
      # should only be one row here, but check just in case
      if math.fabs(lim.quantileExpected+1) < tol:
        mp.obsLimit = lim.limit*mp.totalXSec
        mp.obsLimitErr = lim.limitErr*mp.totalXSec
      treeContents.append([str(lim.quantileExpected),str(lim.limit),str(lim.limitErr)])
    fileHN.Close()
    # check extracted limits
    if mp.obsLimit < 0 or mp.obsLimitErr < 0:
      errorMsg+='Problem in limit tree--one of the obtained observed limits was negative from tree in file: '+fileNameHN+'\n'
      errorMsg+='TTree limit contents:\n'
      errorMsg+='quantileExp\tlimit\tlimitErr\n'
      for listPair in treeContents:
        errorMsg+='\t'+listPair[0]+'\t'+listPair[1]+'\t'+listPair[2]+'\n'
      if len(treeContents) < 1: errorMsg+='\t[empty]\n'
      errorMsg+='Please re-run the limit job: '+shellFileName+'\n\n'
      jobsToRerun+='\t'+shellFileName+'\n'
      filesToRemove+='\trm '+fileNameAsymp+'\n'
      continue
    # write out the result
    fileName=limitsFileNameBase+couplingStr+'.txt'
    with open(fileName,'a') as file:
      mp.Write(file)
  if len(errorMsg) > 0:
    print
    print '----------> ERRORS FOUND:'
    print errorMsg
    if len(filesToRemove) > 0:
      print 'Files to remove:'
      print filesToRemove
    print 'Jobs to re-run:'
    print jobsToRerun
    return False
  else:
    return True


# FIXME: make git-compatible somehow
def GetRooStatsMacroCVSTag(cl95MacroPath):
  cl95Split = cl95MacroPath.split('/')
  cl95MacroName = cl95Split[len(cl95Split)-1]
  cl95MacroDir = cl95MacroPath.rstrip(cl95Split[len(cl95Split)-1])
  proc = subprocess.Popen(['cvs','status',cl95MacroName],cwd=cl95MacroDir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  proc.wait()
  #proc = subprocess.call(['cvs','status',cl95MacroName],cwd=cl95MacroDir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out,err = proc.communicate()
  split = out.split()
  tagVersion = 'unknown'
  for i in range(0,len(split)):
    if 'Tag' in split[i]:
      tagVersion = split[i+1]
      break
  return tagVersion


def PrintInfo(name, listEntries, listErrors):
  print name,"#entries = {",
  for entry in listEntries:
    print "%.5f,"%entry,
  print "}"
  print name,"stat = {",
  for entry in listEntries:
    if entry >= 0:
      print "%.5f,"%math.sqrt(entry),
    else:
      print "NaN,",
  print "}"
  print name,"syst = {",
  for entry in listErrors:
    print "%.5f,"%entry,
  print "}"


def PrintEntries(name, listEntries):
  print name,"#entries = {",
  for entry in listEntries:
    print "%.5f,"%entry,
  print "}"


def OpenFilesAndGetHists(modelPoint):
  global histogramSignalFile
  global histogramBGFile
  global histogramDataFile
  global signalHist
  global signalAccOnlyHistogram
  global signalTotalEventsHistogram
  global signalEntriesTotal
  global signalEntriesTotalUnscaled
  global signalHistSmeared
  global signalHistScaleShiftUp
  global signalHistScaleShiftDown
  global signalHistPileupShiftUp
  global signalHistPileupShiftDown
  global backgroundHist
  global dataHist
  #
  histogramSignalFile = TFile.Open(modelPoint.fileName)
  if not histogramSignalFile:
    print 'file:',modelPoint.fileName,'not found; quitting'
    return
  couplStr = str(modelPoint.coupling).replace('.','')
  massStr = str(modelPoint.mass)
  signalHist = histogramSignalFile.Get("h_kMpl-{coupl}_M-{mass}_{ch}_passKinAndLooseID".format(coupl=couplStr,mass=massStr,ch=modelPoint.channel))
  signalAccOnlyHistogram = histogramSignalFile.Get("h_kMpl-{coupl}_M-{mass}_{ch}_passKin".format(coupl=couplStr,mass=massStr,ch=modelPoint.channel))
  signalNoCutsHistogram = histogramSignalFile.Get("h_kMpl-{coupl}_M-{mass}_noCuts".format(coupl=couplStr,mass=massStr))
  #FIXME: to be added back later?
  #signalTotalEventsHistogram = histogramSignalFile.Get("h_nEvents")
  signalTotalEventsHistogram = TH1F()
  signalEntriesTotal = signalNoCutsHistogram.Integral()
  signalEntriesTotalUnscaled = signalNoCutsHistogram.GetEntries()
  #TODO: get these back later
  #signalHistSmeared = histogramSignalFile.Get("h_Diphoton_Minv_Smeared_FineBinning")
  #signalHistScaleShiftUp = histogramSignalFile.Get("h_Diphoton_Minv_ScaleShiftedUp_FineBinning")
  #signalHistScaleShiftDown = histogramSignalFile.Get("h_Diphoton_Minv_ScaleShiftedDown_FineBinning")
  #signalHistPileupShiftUp = histogramSignalFile.Get("h_Diphoton_Minv_PileupShiftedUp_FineBinning")
  #signalHistPileupShiftDown = histogramSignalFile.Get("h_Diphoton_Minv_PileupShiftedDown_FineBinning")
  signalHistSmeared = signalHist
  signalHistScaleShiftUp = signalHist
  signalHistScaleShiftDown = signalHist
  signalHistPileupShiftUp = signalHist
  signalHistPileupShiftDown = signalHist
  # we've hacked them for now ^^^^^^^^^^^^^^^^^^^^
  # background
  histogramBGFile = TFile.Open(modelPoint.bgFileName)
  if not histogramBGFile:
    print 'file:',modelPoint.bgFileName,'not found; quitting'
    return
  histoBG = histogramBGFile.Get("h_Diphoton_Minv_FineBinning")
  backgroundHist = histoBG.Clone()
  backgroundHist.SetName('backgroundHist')
  # open data file
  histogramDataFile = TFile.Open(modelPoint.dFileName)
  if not histogramDataFile:
    print 'file:',modelPoint.dFileName,'not found; quitting'
    return
  histoData = histogramDataFile.Get("h_Diphoton_Minv_FineBinning")
  dataHist = histoData.Clone()
  dataHist.SetName('dataHist')
  #rebin all
  signalHist = signalHist.Rebin(len(newBins)-1,'',newBinsArr)
  signalAccOnlyHistogram = signalAccOnlyHistogram.Rebin(len(newBins)-1,'',newBinsArr)
  signalNoCutsHistogram = signalNoCutsHistogram.Rebin(len(newBins)-1,'',newBinsArr)
  signalHistSmeared = signalHistSmeared.Rebin(len(newBins)-1,'',newBinsArr)
  signalHistScaleShiftUp = signalHistScaleShiftUp.Rebin(len(newBins)-1,'',newBinsArr)
  signalHistScaleShiftDown = signalHistScaleShiftDown.Rebin(len(newBins)-1,'',newBinsArr)
  signalHistPileupShiftUp = signalHistPileupShiftUp.Rebin(len(newBins)-1,'',newBinsArr)
  signalHistPileupShiftDown = signalHistPileupShiftDown.Rebin(len(newBins)-1,'',newBinsArr)
  backgroundHist = backgroundHist.Rebin(len(newBins)-1,'',newBinsArr)
  dataHist = dataHist.Rebin(len(newBins)-1,'',newBinsArr)


def CloseFilesAndDeleteHists():
  global histogramSignalFile
  global histogramBGFile
  global histogramDataFile
  global signalHist
  global signalAccOnlyHistogram
  global signalTotalEventsHistogram
  global signalEntriesTotal
  global signalHistSmeared
  global signalHistScaleShiftUp
  global signalHistScaleShiftDown
  global signalHistPileupShiftUp
  global signalHistPileupShiftDown
  global backgroundHist
  global dataHist
  histogramSignalFile.Close()
  histogramBGFile.Close()
  histogramDataFile.Close()
  del signalHist
  del signalAccOnlyHistogram
  del signalTotalEventsHistogram
  del signalHistSmeared
  del signalHistScaleShiftUp
  del signalHistScaleShiftDown
  del signalHistPileupShiftUp
  del signalHistPileupShiftDown
  del backgroundHist
  del dataHist


def FillModelPointInfoForWindow(modelPoint,minBin,maxBin):
  optMassRangeLow = dataHist.GetBinLowEdge(minBin)
  optMassRangeHigh = dataHist.GetBinLowEdge(maxBin)+dataHist.GetBinWidth(maxBin)
  entriesData = dataHist.Integral(minBin,maxBin)
  errorStatBG = Double(0)
  entryBG = backgroundHist.IntegralAndError(minBin,maxBin,errorStatBG)
  #if entryBG <= 0:
  #  print 'ERROR: entryBG=',entryBG
  #  print 'backgroundHist.IntegralAndError(minBin,maxBin,errorStatBG)=',backgroundHist.IntegralAndError(minBin,maxBin,errorStatBG)
  #  print 'coupling=',modelPoint.coupling,'mass=',modelPoint.mass,'channel=',modelPoint.channel
  #  print 'optMassRangeLow:',optMassRangeLow,'optMassRangeHigh:',optMassRangeHigh
  #  print 'Opt bin range:',minBin,'-',maxBin
  #  print 'bail out!'
  #  exit(-1)
  #FIXME
  #errSystMC = histosmcUpperError.Integral(minBin,maxBin)
  modelPoint.nDataObs = entriesData
  modelPoint.nBackground = entryBG
  # FIXME: take the upper limit...in case of zero background--1.14?
  modelPoint.nBackgroundErrStat = float(errorStatBG)
  #FIXME
  #modelPoint.nBackgroundErrSyst = math.sqrt(errSystGamJet*errSystGamJet+errSystJetJet*errSystJetJet+errSystMC*errSystMC)
  totalBGErr = math.sqrt(pow(modelPoint.nBackgroundErrStat,2)+pow(modelPoint.nBackgroundErrSyst,2))
  modelPoint.totalSignalEvents = 1.0*signalEntriesTotal
  modelPoint.preMWEff = 1.0*signalHist.Integral()/signalEntriesTotal
  modelPoint.totalEff = 1.0*signalHist.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffErrStat = math.sqrt(modelPoint.totalEff*(1-modelPoint.totalEff)/signalEntriesTotalUnscaled)
  modelPoint.totalEffMScaleSystUp = 1.0*signalHistScaleShiftUp.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMScaleSystDown = 1.0*signalHistScaleShiftDown.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMResSystUp = 1.0*signalHistSmeared.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMResSystDown = 1.0*signalHistSmeared.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffPileupSystUp = 1.0*signalHistPileupShiftUp.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffPileupSystDown = 1.0*signalHistPileupShiftDown.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.optMassWindowLow = optMassRangeLow
  modelPoint.optMassWindowHigh = optMassRangeHigh
  sigMScaleSyst = max(math.fabs(modelPoint.totalEffMScaleSystUp-modelPoint.totalEff)/modelPoint.totalEff,math.fabs(modelPoint.totalEffMScaleSystDown-modelPoint.totalEff)/modelPoint.totalEff)
  sigMResSyst = math.fabs(modelPoint.totalEffMResSystUp-modelPoint.totalEff)/modelPoint.totalEff
  acceptance = 1.0*signalAccOnlyHistogram.Integral()/signalEntriesTotal
  acceptanceMassWindow = 1.0*signalAccOnlyHistogram.Integral(minBin,maxBin)/signalEntriesTotal
  singlePhotonEff = math.sqrt(modelPoint.totalEff/acceptance) # singlePhEff =~ sqrt(diPhotonEff)
  sigEffSystSingleGammaToTwoGamma = (2*singlePhotonEff*math.sqrt(pow(SigScaleFactorSystOneGamma,2)+pow(SigPtSFSystOneGamma,2)))/modelPoint.totalEff
  # above without /totalEff is absolute error on sigma_eff_diphoton; we need % to add with other sigma_eff_diphoton effects here (next line)
  sigEffSyst = math.sqrt(pow(sigMScaleSyst,2)+pow(sigMResSyst,2)+pow(SigPUSyst,2)+pow(SigPDFSyst,2)+pow(sigEffSystSingleGammaToTwoGamma,2))
  # now we make it into absolute error, starting from %
  modelPoint.totalEffErrSyst = sigEffSyst*modelPoint.totalEff
  modelPoint.acceptance = acceptance
  modelPoint.acceptanceMassWindow = acceptanceMassWindow
  totalEffErr = math.sqrt(pow(modelPoint.totalEffErrStat,2)+pow(modelPoint.totalEffErrSyst,2))
  optMassLow = modelPoint.optMassWindowLow
  optMassHigh = modelPoint.optMassWindowHigh
  peakBin = signalHist.GetMaximumBin()
  peakMass = signalHist.GetBinLowEdge(peakBin)
  modelPoint.massPeak = peakMass
  print
  print 'Fill ModelPoint: coupling=',modelPoint.coupling,'mass=',modelPoint.mass,'channel=',modelPoint.channel
  print 'optMassRangeLow:',optMassRangeLow,'optMassRangeHigh:',optMassRangeHigh
  print 'totalSignalEvents:',modelPoint.totalSignalEvents
  print 'Mass peak:',peakMass
  print 'Opt bin range:',minBin,'-',maxBin
  print 'Opt mass window:',modelPoint.optMassWindowLow,'-',modelPoint.optMassWindowHigh
  print 'totalEff =',modelPoint.totalEff,'+/-',modelPoint.totalEffErrStat,'(stat) +/-',modelPoint.totalEffErrSyst,'(syst) = ',totalEffErr
  print 'acceptance =',modelPoint.acceptance
  print 'acceptanceMassWindow =',modelPoint.acceptanceMassWindow
  print 'totalEffMScaleSystDown =',modelPoint.totalEffMScaleSystDown
  print 'totalEffMScaleSystUp =',modelPoint.totalEffMScaleSystUp
  print 'totalEffMResSystDown =',modelPoint.totalEffMResSystDown
  print 'totalEffMResSystUp =',modelPoint.totalEffMResSystUp
  print 'totalEffPileupSystDown =',modelPoint.totalEffPileupSystDown
  print 'totalEffPileupSystUp =',modelPoint.totalEffPileupSystUp
  print 'background = ',modelPoint.nBackground,'+/-',modelPoint.nBackgroundErrStat,'(stat) +/-',modelPoint.nBackgroundErrSyst,'(syst) = ',totalBGErr
  #print 'binRange mScaleDown=',minBinMassScaleSystDown,'-',maxBinMassScaleSystDown
  #print 'binRange mScaleUp=',minBinMassScaleSystUp,'-',maxBinMassScaleSystUp
  print 'entriesNominal=',signalHist.Integral(minBin,maxBin)
  print 'entriesMScaleUp=',signalHistScaleShiftUp.Integral(minBin,maxBin)
  print 'entriesMScaleDown=',signalHistScaleShiftDown.Integral(minBin,maxBin)
  print 'entriesMResSmeared=',signalHistSmeared.Integral(minBin,maxBin)
  print 'entriesPUShiftUp=',signalHistPileupShiftUp.Integral(minBin,maxBin)
  print 'entriesPUShiftDown=',signalHistPileupShiftDown.Integral(minBin,maxBin)


def OptimizeWindow(modelPoint, lumi, maxWindowRange, useAsymmWindow, useSSB, rootPlotFile, imageDir, extraMargin):
  #print 'OptimizeWindow: coupling=',modelPoint.coupling,'mass=',modelPoint.mass,'channel=',modelPoint.channel
  # compute optimization variable for various window sizes
  peakBin = signalHist.GetMaximumBin()
  useAsymmetricWindow = useAsymmWindow
  # maxWindowRange in GeV; get the maximum low bin range in bins
  maxLowBinRangeToUse = peakBin-signalHist.FindBin(signalHist.GetBinCenter(peakBin)-maxWindowRange)
  maxHighBinRangeToUse = signalHist.GetNbinsX()-peakBin+1
  minHalfWindowSize = 10 # GeV
  massRangesUsedForWindow = []
  sOverSqrtBForWindow = []
  ssbForWindow = []
  maxSOverSqrtB = -1
  indexMaxSOverSqrtB = -1
  maxSsb = -1
  indexMaxSsb = -1
  for nBinsMin in xrange(0,maxLowBinRangeToUse):
    minBin = peakBin-nBinsMin
    if minBin < 1:
      break # don't look at underflow (or less) bin
    for nBinsMax in xrange(0,maxHighBinRangeToUse):
      if useAsymmetricWindow:
        maxBin = peakBin+nBinsMax
      else:
        maxBin = peakBin+nBinsMin
      background = backgroundHist.Integral(minBin,maxBin)
      ## for this bin selection, if background is negative or zero, skip this selection
      ##   we could just set it to zero, but it could be the optimal point, which would lead to issues later
      ##   and in any case, this is probably just due to stat. fluctuations
      #if background <= 0:
      #  continue
      # we now bin in variable bins to get around this
      # signal events = (eff*acc)*lumi*crossSec
      #signal = (signalHist.Integral(minBin,maxBin)/signalEntriesTotal)*lumi*modelPoint.totalXSec
      signal = signalHist.Integral(minBin,maxBin)
      if background > 0:
        sOverRootB = signal/math.sqrt(background)
      else:
        sOverRootB = -1
      #print 'minBin=',minBin,'maxBin=',maxBin
      #print 'signal=',signal,'background=',background
      #print 'massRangeUsedForWindow:',signalHist.GetBinLowEdge(minBin),'-',signalHist.GetBinLowEdge(maxBin)+signalHist.GetBinWidth(maxBin)
      #print 'ssb=',signal/math.sqrt(signal+background)
      #modelPoint.Print()
      ssb = signal/math.sqrt(signal+background)
      massRangesUsedForWindow.append((signalHist.GetBinLowEdge(minBin),signalHist.GetBinLowEdge(maxBin)+signalHist.GetBinWidth(maxBin)))
      sOverSqrtBForWindow.append(sOverRootB)
      ssbForWindow.append(ssb)
      fullWindowWidthGeV = signalHist.GetBinLowEdge(maxBin)+signalHist.GetBinWidth(maxBin)-signalHist.GetBinLowEdge(minBin)
      #if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and sOverRootB > maxSOverSqrtB:
      if fullWindowWidthGeV/2 >= minHalfWindowSize and sOverRootB > maxSOverSqrtB:
        maxSOverSqrtB = sOverRootB
        indexMaxSOverSqrtB = len(massRangesUsedForWindow)-1
      #if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and ssb > maxSsb:
      if fullWindowWidthGeV/2 >= minHalfWindowSize and ssb > maxSsb:
        #print '---> SET MAX SSB here'
        maxSsb = ssb
        indexMaxSsb = len(massRangesUsedForWindow)-1
      if maxBin >= signalHist.GetNbinsX():
        break # don't let it go to overflow bin or beyond
      if not useAsymmetricWindow:
        break # break out of the loop at one iteration if we're not using the asymm window
  # add margin to mass window limits
  if useSSB:
    optMassRangeLow = (1-extraMargin)*massRangesUsedForWindow[indexMaxSsb][0]
    optMassRangeHigh = (1+extraMargin)*massRangesUsedForWindow[indexMaxSsb][1]
  else:
    optMassRangeLow = (1-extraMargin)*massRangesUsedForWindow[indexMaxSOverSqrtB][0]
    optMassRangeHigh = (1+extraMargin)*massRangesUsedForWindow[indexMaxSOverSqrtB][1]
  minBin = dataHist.FindBin(optMassRangeLow)
  maxBin = dataHist.FindBin(optMassRangeHigh)-1 # will take the next bin without -1
  print 'opt massRange=',optMassRangeLow,'-',optMassRangeHigh,'(added extra margin of',extraMargin,')'
  print 'opt +extraMargin(if any) binRange =',minBin,'-',maxBin
  # Fill the model point
  #print 'fill ssbForWindow'
  modelPoint.optSSBValue = ssbForWindow[indexMaxSsb]
  #print 'OptimizeWindow: FillModelPointInfoForWindow(): coupling=',modelPoint.coupling,'mass=',modelPoint.mass,'channel=',modelPoint.channel
  FillModelPointInfoForWindow(modelPoint,minBin,maxBin)
  peakMass = signalHist.GetBinLowEdge(peakBin)
  backgroundHist.SetName('backgroundHist')
  rootPlotFile.cd()
  if not rootPlotFile.Get('backgroundHist'):
    backgroundHist.Write()
  signalHist.SetName('diPhotonMinv_'+modelPoint.channel+'_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
  signalHist.Write()
  return peakMass, indexMaxSsb, massRangesUsedForWindow, ssbForWindow, indexMaxSOverSqrtB, sOverSqrtBForWindow


def CollectResults(mp,ssbOpt,peakMass):
  optResultsList.append(mp)


def RunOptimizeJob(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFileNameTemplate,imageDir,extraWindowMargin):
  print 'RunOptimizeJob: Working in process #%d' % os.getpid()
  print '-----> ModelPoint: coupling=',mp.coupling,'mass=',mp.mass,'channel=',mp.channel
  rootFile = TFile(rootFileNameTemplate.format(chan=mp.channel,mass=mp.mass,coup=str(mp.coupling).replace('.','p')),'recreate')
  OpenFilesAndGetHists(mp)
  peakMass, optSsbIndex, massRangesTried, ssbTried, optSRootBIndex, sRootBTried = OptimizeWindow(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFile,imageDir,extraWindowMargin)
  print 'Done with OptimizeWindow() -----> ModelPoint: coupling=',mp.coupling,'mass=',mp.mass,'channel=',mp.channel
  print '-----> opt s/sqrt(s+b)=',ssbTried[optSsbIndex]
  print '-----> opt s/sqrt(b)=',sRootBTried[optSRootBIndex]
  optMassLow = mp.optMassWindowLow
  optMassHigh = mp.optMassWindowHigh
  #optMinMasses.append(optMassLow)
  #optMaxMasses.append(optMassHigh)
  #peakMasses.append(peakMass)
  #optSSBValues.append(ssbTried[optSsbIndex])
  #masses.append(mp.mass)
  rootFile.cd()
  minMassTried = min(zip(*massRangesTried)[0])
  maxMassTried = max(zip(*massRangesTried)[1])
  MakeOptimizationGraph(peakMass,mp,minMassTried,maxMassTried,massRangesTried,ssbTried,optSsbIndex,useAsymmWindow,rootFile)
  MakeOptimizationGraphSRootB(peakMass,mp,minMassTried,maxMassTried,massRangesTried,sRootBTried,optSRootBIndex,useAsymmWindow,rootFile)
  CloseFilesAndDeleteHists()
  rootFile.Close()
  #from pprint import pprint
  #pprint (vars(mp))
  #print mp.__dict__
  #for i in mp.__dict__.itervalues():
  #  print type(i)
  print '-----> Return the result for ModelPoint: coupling=',mp.coupling,'mass=',mp.mass,'channel=',mp.channel
  return mp,float(ssbTried[optSsbIndex]),float(peakMass)


# here we assume if we call this function, it's for a different coupling
def OptimizeSignalMassWindows(modelPointArray,lumi,useAsymmWindow,useSSB,maxWindowRange,txtFile,rootFileNameTemplate,optPlotFile,colorIndex,imageDir,extraWindowMargin):
  pool = multiprocessing.Pool(multiprocessing.cpu_count()) # one job per CPU core (although OS will do the actual schedule on cores)
  asyncResultList = []
  global optimizationCouplingsAlreadyDone
  if modelPointArray[0].coupling in optimizationCouplingsAlreadyDone:
    return
  # loop over model points
  for mp in modelPointArray:
    print
    print 'Optimize:',
    if useSSB:
      print ' Using s/sqrt(s+b),',
    else:
      print ' Using s/sqrt(b),',
    print 'ModelPoint: coupling=',mp.coupling,'mass=',mp.mass,'channel=',mp.channel
    # execute
    #RunOptimizeJob(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFileNameTemplate,imageDir,extraWindowMargin)
    #pool.apply_async(RunOptimizeJob,args=(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFileNameTemplate,imageDir,extraWindowMargin),callback=CollectResults)
    asyncResultList.append(pool.apply_async(RunOptimizeJob,args=(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFileNameTemplate,imageDir,extraWindowMargin)))
  pool.close()
  pool.join() # wait for completion
  optimizationCouplingsAlreadyDone.append(modelPointArray[0].coupling)
  newMPs = []
  for res in asyncResultList:
    #print 'res=',res
    #print '----> res.get()=',res.get()
    mp = res.get()[0]
    #mp.Print()
    #print '-----> Got the result for ModelPoint: coupling=',mp.coupling,'mass=',mp.mass,'channel=',mp.channel
    newMPs.append(mp)
    #print 'added',mp,' from res.get()[0] to newMPs'
    print
  # plotting/smoothing
  # replace our list of modelPoints with the one from the optimize jobs
  modelPointArray = newMPs
  print 'Done with optimize jobs for all model points. Plotting and smoothing.'
  rootFile = optPlotFile
  rootFile.cd()
  # find channels
  channels = GetListOfChannels(modelPointArray);
  #mp.Print()
  # remake arrays above for each channel
  #print 'channels:',channels
  for ch in channels:
    masses = [mp.mass for mp in modelPointArray if mp.channel==ch]
    optMinMasses = [mp.optMassWindowLow for mp in modelPointArray if mp.channel==ch]
    optMaxMasses = [mp.optMassWindowHigh for mp in modelPointArray if mp.channel==ch]
    optSSBValues = [mp.optSSBValue for mp in modelPointArray if mp.channel==ch]
    peakMasses = [mp.massPeak for mp in modelPointArray if mp.channel==ch]
    #print 'masses:',masses
    #print 'optMinMasses:',optMinMasses
    #print 'optMaxMasses:',optMaxMasses
    #print 'peakMasses:',peakMasses
    # make graphs and save to file
    if(len(modelPointArray) > 0):
      coupling = modelPointArray[0].coupling
    else:
      coupling = 0.0
    MakeOptSSBValuesGraph(coupling,masses,ch,optSSBValues,colorIndex,rootFile)
    print 'MakeOptMassWindowGraphs: coupling=',coupling,'masses=',masses,'ch=',ch
    param0,param1 = MakeOptMassWindowGraphs(coupling,masses,ch,optMinMasses,optMaxMasses,peakMasses,colorIndex,rootFile)
    myFunc = TF1("myFunc","pol1")
    myFunc.SetParameters(param0,param1)
    optMinMassesSmoothed = []
    optMaxMassesSmoothed = []
    mpsDone = []
    for mp in modelPointArray:
      if not mp.channel==ch:
        continue
      if mp in mpsDone:
        continue
      mpsDone.append(mp)
      OpenFilesAndGetHists(mp)
      peakBin = signalHist.GetMaximumBin()
      peakMass = signalHist.GetBinLowEdge(peakBin)
      # SIC FIXME: smoothing gets a bit screwed up with variable bins at high masses
      #windowHalfWidth = myFunc.Eval(mp.mass)/2
      #minBin = dataHist.FindBin(peakMass-windowHalfWidth)
      #maxBin = dataHist.FindBin(peakMass+windowHalfWidth)-1 # will take the next bin without -1
      ## Fill the model point
      #print 'smoothed mass window'
      #FillModelPointInfoForWindow(mp,minBin,maxBin)
      #print 'opt smoothed massRange=',mp.optMassWindowLow,'-',mp.optMassWindowHigh
      #print 'opt smoothed binRange =',minBin,'-',maxBin
      mp.Write(txtFile)
      optMinMassesSmoothed.append(mp.optMassWindowLow)
      optMaxMassesSmoothed.append(mp.optMassWindowHigh)
      MakeOptMassWindowSignalBackgroundPlot(rootFile,signalHist,backgroundHist,mp.optMassWindowLow,mp.optMassWindowHigh,mp,imageDir)
      CloseFilesAndDeleteHists()
    # make new half window vs mass
    print 'MakeSmoothedWindowVsMassPlot(coupling='+str(coupling)+','+'masses='+str(masses)+','+'channel='+ch+',...)'
    MakeSmoothedWindowVsMassPlot(coupling,masses,ch,optMinMassesSmoothed,optMaxMassesSmoothed,colorIndex,rootFile)
  return modelPointArray
  


def GetMinMaxMassArrays(modelPointArray, numsigmas = 3):
  # make list of lower/upper edges of mass windows
  minMasses = [(modelPoint.mass - numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  maxMasses = [(modelPoint.mass + numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  return minMasses,maxMasses


# NB: no longer maintained in 2015
## import of yield-calculating code from ExoDiPhotonAnalyzer/test/PlottingCode/CreateHistogramFiles.C
#def CalculateYieldsForMassRanges(HistogramFileLocationData, HistogramFileLocationMC, modelPointArray, lumi, numsigmas, txtFile):

