#!/usr/bin/env python

#
# Support functions for limit setting
#
# Seth I. Cooper, U. Alabama
# November 21 2012

import string
import os
import sys
import glob
import math
import array
import subprocess
import itertools
import tempfile

# run root in batch
sys.argv.append('-b')
from ROOT import *
from ModelPoint import *
from MakePlots import *

# systematics for all model points
SigPUSyst = 0.007
SigPDFSyst = 0.05
SigScaleFactorSyst = 0.005
SigPtSFSyst = 0.03
# background
BGOverallSyst = 0.15

# globally defined hists
histosdata = TH1F()
histosmc = TH1F()
histosJetJet = TH1F()
histosJetJetUpperError = TH1F()
histosJetJetLowerError = TH1F()
histosGammaJet = TH1F()
histosGammaJetUpperError = TH1F()
histosGammaJetLowerError = TH1F()
signalHistogram = TH1F()
signalTotalEventsHistogram = TH1F()
signalEntriesTotal = 0
signalHistogramSmeared = TH1F()
signalHistogramScaleShiftUp = TH1F()
signalHistogramScaleShiftDown = TH1F()
signalHistogramPileupShiftUp = TH1F()
signalHistogramPileupShiftDown = TH1F()
histogramSignalFile = 0


def ComputeLimits(cl95MacroPath,lumi,lumiErr,modelPointArray,fileName,useKFactor):
  if not os.path.isdir('cls_plots'):
    os.mkdir('cls_plots')
  # Get setup script path
  setupScriptPath = cl95MacroPath.split('root')[0]
  setupScriptPath+='setup/lxplus_standalone_setup.sh'
  for modelPoint in modelPointArray:
    # save some info
    thisMPFileName = modelPoint.fileName
    thisMPKFactor = modelPoint.kFactor
    handle,tempFileName = tempfile.mkstemp()
    #void ComputeLimit(float lumi, float lumiError, float totalEff, float totalEffErrStat, float totalEffErrSyst,
    #float nBackground, float nBackgroundErrStat, float nBackgroundErrSyst,
    #float nDataObs, 
    #TODO not needed for computelimit -- save/remove
    #float totalEffMScaleSystUp, float totalEffMScaleSystDown,
    #float totalEffMResSystUp, float totalEffMResSystDown,
    #float totalEffPileupSystUp, float totalEffPileupSystDown,
    #int mass, float coupling, float halfWidth,
    #float totalXSec, int massWindowLow, int massWindowHigh, std::string fileName)
    if useKFactor:
      modelPoint.totalEff*=modelPoint.kFactor
      modelPoint.totalEffErrStat*=modelPoint.kFactor
      modelPoint.totalEffErrSyst*=modelPoint.kFactor
      modelPoint.totalEffMScaleSystUp*=modelPoint.kFactor
      modelPoint.totalEffMScaleSystDown*=modelPoint.kFactor
      modelPoint.totalEffMResSystUp*=modelPoint.kFactor
      modelPoint.totalEffMResSystDown*=modelPoint.kFactor
      modelPoint.totalEffPileupSystUp*=modelPoint.kFactor
      modelPoint.totalEffPileupSystDown*=modelPoint.kFactor
    command = 'ComputeLimit.C('+str(lumi)+','+str(lumiErr)+','+str(modelPoint.totalEff)+','
    command+=str(modelPoint.totalEffErrStat)+','+str(modelPoint.totalEffErrSyst)+','
    command+=str(modelPoint.totalEffMScaleSystUp)+','+str(modelPoint.totalEffMScaleSystDown)+','
    command+=str(modelPoint.totalEffMResSystUp)+','+str(modelPoint.totalEffMResSystDown)+','
    command+=str(modelPoint.totalEffPileupSystUp)+','+str(modelPoint.totalEffPileupSystDown)+','
    command+=str(modelPoint.nBackground)+','+str(modelPoint.nBackgroundErrStat)+','+str(modelPoint.nBackgroundErrSyst)+','
    command+=str(modelPoint.nDataObs)+','+str(modelPoint.mass)+','
    command+=str(modelPoint.coupling)+','+str(modelPoint.halfWidth)+','+str(modelPoint.totalXSec)+','
    command+=str(modelPoint.optMassWindowLow)+','+str(modelPoint.optMassWindowHigh)+','
    command+='\"'+tempFileName+'\"'
    command+=')'
    subprocess.call(['root','-b','-q',command]) # wait for each one to finish
    os.close(handle)
    couplingStr = str(modelPoint.coupling).replace('.','p')
    plotFileSamplingName = 'cls_sampling_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    plotFileName = 'cls_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    subprocess.call(['mv','cls_sampling.pdf','cls_plots/'+plotFileSamplingName])
    subprocess.call(['mv','cls.pdf','cls_plots/'+plotFileName])
    # read this model point (now with limits) from temp file
    file = open(tempFileName, 'r')
    lines = file.readlines()
    file.close()
    os.remove(tempFileName)
    thisModelPoint = ReadFromLines(lines)[0]
    thisModelPoint.fileName = thisMPFileName
    thisModelPoint.kFactor = thisMPKFactor
    #FIXME/TODO: store all info besides limit info like this and don't use ComputeLimit.C to write it out...
    # replace this model point in the array
    for index,mp in enumerate(modelPointArray):
      if mp.mass==thisModelPoint.mass and mp.coupling==thisModelPoint.coupling:
        modelPointArray[index] = thisModelPoint
    # write out the result file
    if os.path.isfile(fileName):
      os.remove(fileName)
    with open(fileName,'w') as file:
      for mp in modelPointArray:
        mp.Write(file)
  print 'Used StatisticalTools/RooStatsRoutine CVS tag:',GetRooStatsMacroCVSTag(cl95MacroPath)


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


def OpenSignalFilesAndGetHists(modelPoint):
  global histogramSignalFile
  global signalHistogram
  global signalTotalEventsHistogram
  global signalEntriesTotal
  global signalHistogramSmeared
  global signalHistogramScaleShiftUp
  global signalHistogramScaleShiftDown
  global signalHistogramPileupShiftUp
  global signalHistogramPileupShiftDown
  #
  histogramSignalFile = TFile.Open(modelPoint.fileName)
  if not histogramSignalFile:
    print 'file:',modelPoint.fileName,'not found; quitting'
    return
  signalHistogram = histogramSignalFile.Get("h_Diphoton_Minv_FineBinning")
  signalTotalEventsHistogram = histogramSignalFile.Get("h_nEvents")
  signalEntriesTotal = signalTotalEventsHistogram.Integral()
  signalHistogramSmeared = histogramSignalFile.Get("h_Diphoton_Minv_Smeared_FineBinning")
  signalHistogramScaleShiftUp = histogramSignalFile.Get("h_Diphoton_Minv_ScaleShiftedUp_FineBinning")
  signalHistogramScaleShiftDown = histogramSignalFile.Get("h_Diphoton_Minv_ScaleShiftedDown_FineBinning")
  signalHistogramPileupShiftUp = histogramSignalFile.Get("h_Diphoton_Minv_PileupShiftedUp_FineBinning")
  signalHistogramPileupShiftDown = histogramSignalFile.Get("h_Diphoton_Minv_PileupShiftedDown_FineBinning")


def CloseSignalFilesAndDeleteHists():
  global histogramSignalFile
  global signalHistogram
  global signalTotalEventsHistogram
  global signalHistogramSmeared
  global signalHistogramScaleShiftUp
  global signalHistogramScaleShiftDown
  global signalHistogramPileupShiftUp
  global signalHistogramPileupShiftDown
  histogramSignalFile.Close()
  del signalHistogram
  del signalTotalEventsHistogram
  del signalHistogramSmeared
  del signalHistogramScaleShiftUp
  del signalHistogramScaleShiftDown
  del signalHistogramPileupShiftUp
  del signalHistogramPileupShiftDown


def FillModelPointInfoForWindow(modelPoint,minBin,maxBin):
  optMassRangeLow = histosdata.GetBinLowEdge(minBin)
  optMassRangeHigh = histosdata.GetBinLowEdge(maxBin)+histosdata.GetBinWidth(maxBin)
  entriesData = histosdata.Integral(minBin,maxBin)
  errGamJet = Double(0) #ROOT.Double()
  entGamJet = histosGammaJet.IntegralAndError(minBin,maxBin,errGamJet)
  errJetJet = Double(0)
  entJetJet = histosJetJet.IntegralAndError(minBin,maxBin,errJetJet)
  errMC = Double(0)
  entMC = histosmc.IntegralAndError(minBin,maxBin,errMC)
  entryBG = entGamJet+entJetJet+entMC
  errorBG = math.sqrt(errGamJet*errGamJet+errJetJet*errJetJet+errMC*errMC)
  modelPoint.nDataObs = entriesData
  modelPoint.nBackground = entryBG
  # FIXME: take the upper limit...in case of zero background--1.14?
  modelPoint.nBackgroundErrStat = errorBG
  modelPoint.totalEff = 1.0*signalHistogram.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffErrStat = math.sqrt(modelPoint.totalEff*(1-modelPoint.totalEff)/signalEntriesTotal)
  modelPoint.totalEffMScaleSystUp = 1.0*signalHistogramScaleShiftUp.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMScaleSystDown = 1.0*signalHistogramScaleShiftDown.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMResSystUp = 1.0*signalHistogramSmeared.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffMResSystDown = 1.0*signalHistogramSmeared.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffPileupSystUp = 1.0*signalHistogramPileupShiftUp.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.totalEffPileupSystDown = 1.0*signalHistogramPileupShiftDown.Integral(minBin,maxBin)/signalEntriesTotal
  modelPoint.optMassWindowLow = optMassRangeLow
  modelPoint.optMassWindowHigh = optMassRangeHigh
  sigMScaleSyst = max(math.fabs(modelPoint.totalEffMScaleSystUp-modelPoint.totalEff)/modelPoint.totalEff,math.fabs(modelPoint.totalEffMScaleSystDown-modelPoint.totalEff)/modelPoint.totalEff)
  sigMResSyst = math.fabs(modelPoint.totalEffMResSystUp-modelPoint.totalEff)/modelPoint.totalEff
  sigEffSyst = math.sqrt(pow(sigMScaleSyst,2)+pow(sigMResSyst,2)+pow(SigPUSyst,2)+pow(SigPDFSyst,2)+pow(SigScaleFactorSyst,2)+pow(SigPtSFSyst,2))
  modelPoint.totalEffErrSyst = sigEffSyst*modelPoint.totalEff # make into number of events, not %
  totalEffErr = math.sqrt(pow(modelPoint.totalEffErrStat,2)+pow(modelPoint.totalEffErrSyst,2))
  bgErrSyst = BGOverallSyst
  modelPoint.nBackgroundErrSyst = bgErrSyst*modelPoint.nBackground # make into number of events, not %
  totalBGErr = math.sqrt(pow(modelPoint.nBackgroundErrStat,2)+pow(modelPoint.nBackgroundErrSyst,2))
  optMassLow = modelPoint.optMassWindowLow
  optMassHigh = modelPoint.optMassWindowHigh
  peakBin = signalHistogram.GetMaximumBin()
  peakMass = signalHistogram.GetBinLowEdge(peakBin)
  print 'Fill model point coupling=',modelPoint.coupling,'mass=',modelPoint.mass
  print 'Mass peak:',peakMass
  print 'totalEff =',modelPoint.totalEff,'+/-',modelPoint.totalEffErrStat,'(stat) +/-',modelPoint.totalEffErrSyst,'(syst) = ',totalEffErr
  print 'totalEffMScaleSystDown =',modelPoint.totalEffMScaleSystDown
  print 'totalEffMScaleSystUp =',modelPoint.totalEffMScaleSystUp
  print 'totalEffMResSystDown =',modelPoint.totalEffMResSystDown
  print 'totalEffMResSystUp =',modelPoint.totalEffMResSystUp
  print 'totalEffPileupSystDown =',modelPoint.totalEffPileupSystDown
  print 'totalEffPileupSystUp =',modelPoint.totalEffPileupSystUp
  print 'background = ',modelPoint.nBackground,'+/-',modelPoint.nBackgroundErrStat,'(stat) +/-',modelPoint.nBackgroundErrSyst,'(syst) = ',totalBGErr
  #print 'binRange mScaleDown=',minBinMassScaleSystDown,'-',maxBinMassScaleSystDown
  #print 'binRange mScaleUp=',minBinMassScaleSystUp,'-',maxBinMassScaleSystUp
  print 'entriesNominal=',signalHistogram.Integral(minBin,maxBin)
  print 'entriesMScaleUp=',signalHistogramScaleShiftUp.Integral(minBin,maxBin)
  print 'entriesMScaleDown=',signalHistogramScaleShiftDown.Integral(minBin,maxBin)
  print 'entriesMResSmeared=',signalHistogramSmeared.Integral(minBin,maxBin)
  print 'entriesPUShiftUp=',signalHistogramPileupShiftUp.Integral(minBin,maxBin)
  print 'entriesPUShiftDown=',signalHistogramPileupShiftDown.Integral(minBin,maxBin)
  print


def OptimizeWindow(modelPoint, lumi, maxWindowRange, useAsymmWindow, useSSB, rootPlotFile, imageDir, extraMargin):
  # compute optimization variable for various window sizes
  peakBin = signalHistogram.GetMaximumBin()
  useAsymmetricWindow = useAsymmWindow
  maxLowBinRangeToUse = maxWindowRange # bins/GeV
  maxHighBinRangeToUse = signalHistogram.GetNbinsX()-peakBin+1
  minHalfWindowSize = 10 # bins
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
        # loop is redundant here
      entGamJet = histosGammaJet.Integral(minBin,maxBin)
      entJetJet = histosJetJet.Integral(minBin,maxBin)
      entMC = histosmc.Integral(minBin,maxBin)
      background = entGamJet+entJetJet+entMC
      # signal events = (eff*acc)*lumi*crossSec
      signal = (signalHistogram.Integral(minBin,maxBin)/signalEntriesTotal)*lumi*modelPoint.totalXSec
      if background > 0:
        sOverRootB = signal/math.sqrt(background)
      else:
        sOverRootB = -1
      ssb = signal/math.sqrt(signal+background)
      massRangesUsedForWindow.append((signalHistogram.GetBinLowEdge(minBin),signalHistogram.GetBinLowEdge(maxBin)+signalHistogram.GetBinWidth(maxBin)))
      sOverSqrtBForWindow.append(sOverRootB)
      ssbForWindow.append(ssb)
      if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and sOverRootB > maxSOverSqrtB:
        maxSOverSqrtB = sOverRootB
        indexMaxSOverSqrtB = len(massRangesUsedForWindow)-1
      if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and ssb > maxSsb:
        maxSsb = ssb
        indexMaxSsb = len(massRangesUsedForWindow)-1
      if maxBin >= signalHistogram.GetNbinsX():
        break # don't let it go to overflow bin or beyond
  # add margin to mass window limits
  if useSSB:
    optMassRangeLow = (1-extraMargin)*massRangesUsedForWindow[indexMaxSsb][0]
    optMassRangeHigh = (1+extraMargin)*massRangesUsedForWindow[indexMaxSsb][1]
  else:
    optMassRangeLow = (1-extraMargin)*massRangesUsedForWindow[indexMaxSOverSqrtB][0]
    optMassRangeHigh = (1+extraMargin)*massRangesUsedForWindow[indexMaxSOverSqrtB][1]
  minBin = histosdata.FindBin(optMassRangeLow)
  maxBin = histosdata.FindBin(optMassRangeHigh)-1 # will take the next bin without -1
  # Fill the model point
  FillModelPointInfoForWindow(modelPoint,minBin,maxBin)
  print 'opt massRange=',optMassRangeLow,'-',optMassRangeHigh,'(added extra margin of',extraMargin,')'
  print 'opt +extraMargin(if any) binRange =',minBin,'-',maxBin
  peakMass = signalHistogram.GetBinLowEdge(peakBin)
  backgroundHist = histosGammaJet.Clone()
  backgroundHist.Add(histosJetJet)
  backgroundHist.Add(histosmc)
  backgroundHist.SetName('backgroundHist')
  rootPlotFile.cd()
  if not rootPlotFile.Get('backgroundHist'):
    backgroundHist.Write()
  signalHistogram.SetName('diPhotonMinv_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
  signalHistogram.Write()
  return peakMass, indexMaxSsb, massRangesUsedForWindow, ssbForWindow, indexMaxSOverSqrtB, sOverSqrtBForWindow


def OptimizeSignalMassWindows(HistogramFileLocationData,HistogramFileLocationMC,modelPointArray,lumi,useAsymmWindow,useSSB,maxWindowRange,txtFile,rootFile,colorIndex,imageDir,extraWindowMargin,DataSample):
  optMinMasses = []
  optMaxMasses = []
  peakMasses = []
  optSSBValues = []
  optSRootBValues = []
  masses = []

  # open MC files
  global fJetJethists
  global fGammaJethists
  global fMChists
  global histosmc
  global histosJetJet
  global histosJetJetUpperError
  global histosJetJetLowerError
  global histosGammaJet
  global histosGammaJetUpperError
  global histosGammaJetLowerError
  global fdatahists
  global histosdata
  # FIXME hardcoded names/locations
  # background
  print 'using data sample:',DataSample+';','histogram file location for fakes/data:',HistogramFileLocationData
  print 'histogram file for backgroundMC:',HistogramFileLocationMC
  histogramFileJetJet = HistogramFileLocationData+"/histograms_"+DataSample+"_JetJet.root"
  histogramFileGammaJet = HistogramFileLocationData+"/histograms_"+DataSample+"_GammaJet.root"	
  histogramFileMC = HistogramFileLocationMC
  fJetJethists = TFile.Open(histogramFileJetJet)
  if not fJetJethists:
    print 'file:',histogramFileJetJet,'not found; quitting'
    return
  fGammaJethists = TFile.Open(histogramFileGammaJet)
  if not fGammaJethists:
    print 'file:',histogramFileGammaJet,'not found; quitting'
    return
  fMChists = TFile.Open(histogramFileMC)
  if not fMChists:
    print 'file:',histogramFileMC,'not found; quitting'
    return
  histosmc = fMChists.Get("h_Diphoton_Minv_FineBinning")
  histosmc.Scale(lumi)
  histosJetJet = fJetJethists.Get("h_JetJet_minv_FineBinning")
  histosJetJetUpperError = fJetJethists.Get("h_JetJet_minv_UpperError")
  histosJetJetLowerError = fJetJethists.Get("h_JetJet_minv_LowerError")
  histosGammaJet = fGammaJethists.Get("h_GammaJet_minv_FineBinning")
  histosGammaJetUpperError = fGammaJethists.Get("h_GammaJet_minv_UpperError")
  histosGammaJetLowerError = fGammaJethists.Get("h_GammaJet_minv_LowerError")
  backgroundHist = histosGammaJet.Clone()
  backgroundHist.Add(histosJetJet)
  backgroundHist.Add(histosmc)
  backgroundHist.SetName('backgroundHist')
  histogramFileData = HistogramFileLocationData+"/histograms_"+DataSample+".root"
  # open data file
  fdatahists = TFile.Open(histogramFileData)
  if not fdatahists:
    print 'file:',histogramFileData,'not found; quitting'
    return
  histosdata = fdatahists.Get("h_Diphoton_Minv_FineBinning")

  # loop over model points
  for mp in modelPointArray:
    print
    print 'Optimize:',
    if useSSB:
      print ' Using s/sqrt(s+b),',
    else:
      print ' Using s/sqrt(b),',
    print 'ModelPoint: coupling=',mp.coupling,'mass=',mp.mass
    OpenSignalFilesAndGetHists(mp)
    peakMass, optSsbIndex, massRangesTried, ssbTried, optSRootBIndex, sRootBTried = OptimizeWindow(mp,lumi,maxWindowRange,useAsymmWindow,useSSB,rootFile,imageDir,extraWindowMargin)
    print 'opt s/sqrt(s+b)=',ssbTried[optSsbIndex]
    print 'opt s/sqrt(b)=',sRootBTried[optSRootBIndex]
    optMassLow = mp.optMassWindowLow
    optMassHigh = mp.optMassWindowHigh
    optMinMasses.append(optMassLow)
    optMaxMasses.append(optMassHigh)
    peakMasses.append(peakMass)
    optSSBValues.append(ssbTried[optSsbIndex])
    masses.append(mp.mass)
    rootFile.cd()
    minMassTried = min(zip(*massRangesTried)[0])
    maxMassTried = max(zip(*massRangesTried)[1])
    MakeOptimizationGraph(peakMass,mp,minMassTried,maxMassTried,massRangesTried,ssbTried,optSsbIndex,useAsymmWindow,rootFile)
    MakeOptimizationGraphSRootB(peakMass,mp,minMassTried,maxMassTried,massRangesTried,sRootBTried,optSRootBIndex,useAsymmWindow,rootFile)
    CloseSignalFilesAndDeleteHists()
  # make graphs and save to file
  if(len(modelPointArray) > 0):
    coupling = modelPointArray[0].coupling
  else:
    coupling = 0.0
  MakeOptSSBValuesGraph(coupling,masses,optSSBValues,colorIndex,rootFile)
  param0,param1 = MakeOptMassWindowGraphs(coupling,masses,optMinMasses,optMaxMasses,peakMasses,colorIndex,rootFile)
  myFunc = TF1("myFunc","pol1")
  myFunc.SetParameters(param0,param1)
  # now update mass windows with smoothing from linear fit
  optMinMassesSmoothed = []
  optMaxMassesSmoothed = []
  for mp in modelPointArray:
    OpenSignalFilesAndGetHists(mp)
    peakBin = signalHistogram.GetMaximumBin()
    peakMass = signalHistogram.GetBinLowEdge(peakBin)
    windowHalfWidth = myFunc.Eval(mp.mass)/2
    minBin = histosdata.FindBin(peakMass-windowHalfWidth)
    maxBin = histosdata.FindBin(peakMass+windowHalfWidth)-1 # will take the next bin without -1
    # Fill the model point
    print 'smoothed mass window'
    FillModelPointInfoForWindow(mp,minBin,maxBin)
    print 'opt smoothed massRange=',mp.optMassWindowLow,'-',mp.optMassWindowHigh
    print 'opt smoothed binRange =',minBin,'-',maxBin
    mp.Write(txtFile)
    optMinMassesSmoothed.append(mp.optMassWindowLow)
    optMaxMassesSmoothed.append(mp.optMassWindowHigh)
    MakeOptMassWindowSignalBackgroundPlot(rootFile,signalHistogram,backgroundHist,mp.optMassWindowLow,mp.optMassWindowHigh,mp,imageDir)
    CloseSignalFilesAndDeleteHists()
  # make new half window vs mass
  MakeSmoothedWindowVsMassPlot(coupling,masses,optMinMassesSmoothed,optMaxMassesSmoothed,colorIndex,rootFile)
  
  # close mc/data
  fJetJethists.Close()
  fGammaJethists.Close()
  fMChists.Close()
  fdatahists.Close()



def GetMinMaxMassArrays(modelPointArray, numsigmas = 3):
  # make list of lower/upper edges of mass windows
  minMasses = [(modelPoint.mass - numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  maxMasses = [(modelPoint.mass + numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  return minMasses,maxMasses


# import of yield-calculating code from ExoDiPhotonAnalyzer/test/PlottingCode/CreateHistogramFiles.C
def CalculateYieldsForMassRanges(HistogramFileLocationData, HistogramFileLocationMC, modelPointArray, lumi, numsigmas, txtFile, DataSample):
  gStyle.SetOptStat("ourmei")

  # FIXME hardcoded names/locations
  # background
  histogramFileData = HistogramFileLocationData+"/histograms_"+DataSample+".root"
  histogramFileJetJet = HistogramFileLocationData+"/histograms_"+DataSample+"_JetJet.root"
  histogramFileGammaJet = HistogramFileLocationData+"/histograms_"+DataSample+"_GammaJet.root"	
  histogramFileMC = HistogramFileLocationMC
  fdatahists = TFile.Open(histogramFileData)
  if not fdatahists:
    print 'file:',histogramFileData,'not found; quitting'
    return
  fJetJethists = TFile.Open(histogramFileJetJet)
  if not fJetJethists:
    print 'file:',histogramFileJetJet,'not found; quitting'
    return
  fGammaJethists = TFile.Open(histogramFileGammaJet)
  if not fGammaJethists:
    print 'file:',histogramFileGammaJet,'not found; quitting'
    return
  fMChists = TFile.Open(histogramFileMC)
  if not fMChists:
    print 'file:',histogramFileMC,'not found; quitting'
    return

  #print "Getting data histogram"
  histosdata = fdatahists.Get("h_Diphoton_Minv_FineBinning")
  #print "histo",histosdata.GetName(),histosdata.GetEntries(),"entries",histosdata.Integral(),"(integral)"

  #print "Getting MC histogram"
  histosmc = fMChists.Get("h_Diphoton_Minv_FineBinning")
  histosmc.Scale(lumi)
  #print "histo",histosmc.GetName(),histosmc.GetEntries(),"entries",histosmc.Integral(),"(integral)"

  #print "Getting Jet+Jet histograms"
  histosJetJet = fJetJethists.Get("h_JetJet_minv_FineBinning")
  histosJetJetUpperError = fJetJethists.Get("h_JetJet_minv_UpperError")
  histosJetJetLowerError = fJetJethists.Get("h_JetJet_minv_LowerError")
  #print "histo",histosJetJet.GetName(),histosJetJet.GetEntries(),"entries",histosJetJet.Integral(),"(integral)"
  #print "histo",histosJetJetUpperError.GetName(),histosJetJetUpperError.GetEntries(),"entries",histosJetJetUpperError.Integral(),"(integral)"
  #print "histo",histosJetJetLowerError.GetName(),histosJetJetLowerError.GetEntries(),"entries",histosJetJetLowerError.Integral(),"(integral)"

  #print "Getting Gamma+Jet histograms"
  histosGammaJet = fGammaJethists.Get("h_GammaJet_minv_FineBinning")
  histosGammaJetUpperError = fGammaJethists.Get("h_GammaJet_minv_UpperError")
  histosGammaJetLowerError = fGammaJethists.Get("h_GammaJet_minv_LowerError")
  #print "histo",histosGammaJet.GetName(),histosGammaJet.GetEntries(),"entries",histosGammaJet.Integral(),"(integral)"
  #print "histo",histosGammaJetUpperError.GetName(),histosGammaJetUpperError.GetEntries(),"entries",histosGammaJetUpperError.Integral(),"(integral)"
  #print "histo",histosGammaJetLowerError.GetName(),histosGammaJetLowerError.GetEntries(),"entries",histosGammaJetLowerError.Integral(),"(integral)"

  ## make list of lower/upper edges of mass windows
  minMasses,maxMasses = GetMinMaxMassArrays(modelPointArray, numsigmas)

  # calculate integrals in different mass ranges
  # They all have the same binning
  # Let's take the data one
  nbinsX = histosdata.GetNbinsX()
  binnr = nbinsX
  binminnumbers = []
  binmaxnumbers = []

  print "---------------Mass intervals to seek----------------"
# debug
  for minMass,maxMass in itertools.izip(minMasses,maxMasses):
    print "interval [",minMass,",",maxMass,"]"
    binminnumbers.append(histosdata.FindBin(minMass)+1)
    maxMassBin = histosdata.FindBin(maxMass)
    if maxMassBin == histosdata.GetNbinsX()+1: # is it the overflow bin?
      maxMassBin = histosdata.GetNbinsX()
    binmaxnumbers.append(maxMassBin)

  #Once we have the bin numbers
  #we compute the integrals
  #and their stat uncertainties
  #and their syst uncertainties
  #For each contribution, we take nbinsX+1 to have the overflow

  for binMin,binMax in itertools.izip(binminnumbers,binmaxnumbers):
    print 'binNumber range:',binMin,"---",binMax

  entriesJetJet = [histosJetJet.Integral(binMin,binMax) for binMin, binMax in zip(binminnumbers,binmaxnumbers)]
  errorsJetJet = [histosJetJetUpperError.Integral(binMin,binMax) for binMin, binMax in zip(binminnumbers,binmaxnumbers)]
  PrintInfo('JetJet',entriesJetJet,errorsJetJet)

  entriesGammaJet = [histosGammaJet.Integral(binMin,binMax) for binMin, binMax in zip(binminnumbers,binmaxnumbers)]
  errorsGammaJet = [histosGammaJetUpperError.Integral(binMin,binMax) for binMin, binMax in zip(binminnumbers,binmaxnumbers)]
  PrintInfo('GammaJet',entriesGammaJet,errorsGammaJet)

  entriesFake = [gammaJet+jetJet for gammaJet, jetJet in zip(entriesGammaJet,entriesJetJet)]
  errorsFake = [math.sqrt((errGammaJet*errGammaJet) + (errJetJet*errJetJet)) for errGammaJet, errJetJet in zip(errorsGammaJet,errorsJetJet)]
  PrintInfo('Fake',entriesFake,errorsFake)

  entriesMC = [histosmc.Integral(binMin,binMax) for binMin, binMax in zip(binminnumbers,binmaxnumbers)]
  errorsMC = [0]*len(binminnumbers)
  PrintInfo('MC',entriesMC,errorsMC)

  entriesBackground = [entGamJet+entJetJet+entMC for entGamJet,entJetJet,entMC in zip(entriesGammaJet,entriesJetJet,entriesMC)]
  errorsBackground = [sqrt((errGamJet*errGamJet) + (errJetJet*errJetJet) + (errMC*errMC)) for errGamJet,errJetJet,errMC in zip(errorsGammaJet,errorsJetJet,errorsMC)]
  PrintInfo('Background',entriesBackground,errorsBackground)

  entriesData = [histosdata.Integral(binMin,binMax) for binMin,binMax in zip(binminnumbers,binmaxnumbers)]
  errorsData = [0]*len(binminnumbers)
  PrintInfo('Data',entriesData,errorsData)

  entriesDataOverBackground = [entData/entBg if entBg > 0 else 'bg=0' for entData,entBg in zip(entriesData,entriesBackground)]
  errorsDataOverBackground = []
  for entData,entBg,errData,errBg in itertools.izip(entriesData,entriesBackground,errorsData,errorsBackground):
    if entData != 0:
      errorsDataOverBackground.append(entData/entBg * sqrt( (errData/entData)*(errData/entData) + (errBg/entBg)*(errBg/entBg) ))
    else:
      errorsDataOverBackground.append(-1)
  # note: printinfo doesn't work here since 'bg=0' isn't a float
  #PrintInfo('Data/Background',entriesDataOverBackground,errorsDataOverBackground)
  
  entriesSignalInWindow = []
  entriesSignalTotal = []
  signalCouplings = []
  signalMasses = []

  # Now put the info into the ModelPoints
  for mp, entryData, entryBG, errorBG, binMin, binMax, minMass, maxMass in itertools.izip(modelPointArray,entriesData,entriesBackground,errorsBackground,binminnumbers,binmaxnumbers,minMasses,maxMasses):
    mp.nDataObs = entryData
    mp.nBackground = entryBG
    mp.nBackgroundErrStat = math.sqrt(entryBG) #FIXME
    mp.optMassWindowLow = minMass
    mp.optMassWindowHigh = maxMass
    # signal
    histogramSignalFile = TFile.Open(mp.fileName)
    if not histogramSignalFile:
      print 'file:',mp.fileName,'not found; quitting'
      return
    signalHistogram = histogramSignalFile.Get("h_Diphoton_Minv_FineBinning")
    signalEntriesInWindow = signalHistogram.Integral(binMin,binMax)
    signalTotalEventsHistogram = histogramSignalFile.Get("h_nEvents")
    signalEntriesTotal = signalTotalEventsHistogram.Integral()
    histogramSignalFile.Close()
    mp.totalEff = 1.0*signalEntriesInWindow/signalEntriesTotal
    entriesSignalInWindow.append(signalEntriesInWindow)
    entriesSignalTotal.append(signalEntriesTotal)
    signalCouplings.append(mp.coupling)
    signalMasses.append(mp.mass)
    mp.Write(txtFile)

  PrintEntries('Signal couplings',signalCouplings)
  PrintEntries('Signal Masses',signalMasses)
  PrintEntries('signalEntriesInWindow',entriesSignalInWindow)
  PrintEntries('signalEntriesTotal',entriesSignalTotal)
  signalEffs = [sigEvtsWindow/sigEvtsTotal for sigEvtsWindow,sigEvtsTotal in itertools.izip(entriesSignalInWindow,entriesSignalTotal)]
  PrintEntries('Signal efficiencies*acceptances',signalEffs)


