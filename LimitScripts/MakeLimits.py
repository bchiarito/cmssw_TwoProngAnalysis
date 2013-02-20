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

# run root in batch
sys.argv.append('-b')
from ROOT import *
from ModelPoint import *


def ComputeLimits(cl95MacroPath,lumi,lumiErr,modelPointArray,fileName):
  # check to see if results txt file exists; if so, remove
  if os.path.isfile(fileName):
    os.remove(fileName)
  # Get setup script path
  setupScriptPath = cl95MacroPath.split('root')[0]
  setupScriptPath+='setup/lxplus_standalone_setup.sh'
  for modelPoint in modelPointArray:
    # void ComputeLimit(float lumi, float lumiError, float totalEff, float nBackground, float nBackgroundErr,
    # float nDataObs, int mass, float coupling, float halfWidth, float totalXSec, std::string fileName)
    command = 'ComputeLimit.C('+str(lumi)+','+str(lumiErr)+','+str(modelPoint.totalEff)+','+str(modelPoint.nBackground)+','
    command+=str(modelPoint.nBackgroundErr)+','+str(modelPoint.nDataObs)+','+str(modelPoint.mass)+','
    command+=str(modelPoint.coupling)+','+str(modelPoint.halfWidth)+','+str(modelPoint.totalXSec)+','
    command+='\"'+fileName+'\"'
    command+=')'
    subprocess.call(['root','-b','-q',command]) # wait for each one to finish
    couplingStr = str(modelPoint.coupling).replace('.','p')
    plotFileSamplingName = 'cls_sampling_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    plotFileName = 'cls_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    if not os.path.isdir('cls_plots'):
      os.mkdir('cls_plots')
    subprocess.call(['mv','cls_sampling.pdf','cls_plots/'+plotFileSamplingName])
    subprocess.call(['mv','cls.pdf','cls_plots/'+plotFileName])
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


def OptimizeWindow(HistogramFileLocation, modelPoint, lumi):
  #Sample = "ExoDiPhotonAnalyzer_DataABC"
  Sample = "ExoDiPhotonAnalyzer_PFDec14th_DataABCD"

  # FIXME hardcoded names/locations
  # background
  histogramFileJetJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root"
  histogramFileGammaJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root"	
  histogramFileMC = HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root"
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

  # signal
  histogramSignalFile = TFile.Open(modelPoint.fileName)
  if not histogramSignalFile:
    print 'file:',modelPoint.fileName,'not found; quitting'
    return
  signalHistogram = histogramSignalFile.Get("h_Diphoton_Minv_FineBinning")
  signalTotalEventsHistogram = histogramSignalFile.Get("h_nEvents")
  signalEntriesTotal = signalTotalEventsHistogram.Integral()

  # compute optimization variable for various window sizes
  useAsymmetricWindow = False
  minHalfWindowSize = 10 # bins
  maxBinRangeToUse = 400 # bins
  nBinsUsedForWindow = []
  sOverSqrtBForWindow = []
  ssbForWindow = []
  peakBin = signalHistogram.FindBin(signalHistogram.GetMaximumBin())
  maxSOverSqrtB = -1
  indexMaxSOverSqrtB = -1
  maxSsb = -1
  indexMaxSsb = -1
  for nBinsMin in xrange(0,maxBinRangeToUse):
    minBin = peakBin-nBinsMin
    for nBinsMax in xrange(0,maxBinRangeToUse):
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
      nBinsUsedForWindow.append((nBinsMin,maxBin-peakBin))
      sOverSqrtBForWindow.append(sOverRootB)
      ssbForWindow.append(ssb)
      if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and sOverRootB > maxSOverSqrtB:
        maxSOverSqrtB = sOverRootB
        indexMaxSOverSqrtB = len(nBinsUsedForWindow)-1
      if (nBinsMin+maxBin-peakBin)/2 >= minHalfWindowSize and ssb > maxSsb:
        maxSsb = ssb
        indexMaxSsb = len(nBinsUsedForWindow)-1
      # INFO
      ##print 'nBins=',nBins,'binMin=',minBin,'binMax=',maxBin
      ##print 'mass range=',histosmc.FindBin(minBin),'-',histosmc.FindBin(maxBin)
      #print 'entGamJet=',entGamJet,'entJetJet=',entJetJet,'entMC',entMC
      #print 'background=',background,'signal=',signal,'s/sqrt(b)=',sOverRootB
      ##print 'background=',background,'signal=',signal,'Z_Bi=',zbi
      #if modelPoint.mass==3000:
      #  print 'background=',background,'signal=',signal,'s/sqrt(s+b)=',ssb

  # INFO
  #print 'modelPoint k=',modelPoint.coupling,'mass=',modelPoint.mass
  #print 'Opt. window halfsize =',nBinsForMaxSsb
  #print 'Opt. mass range',histosmc.GetBinLowEdge(histosmc.FindBin(modelPoint.mass)-nBinsForMaxSOverSqrtB),'-',histosmc.GetBinLowEdge(histosmc.FindBin(modelPoint.mass)+nBinsForMaxSOverSqrtB)
  #print 'Opt. mass range',histosmc.GetBinLowEdge(peakBin-nBinsForMaxSsb),'-',histosmc.GetBinLowEdge(peakBin+nBinsForMaxSsb)

  # now get data yields
  histogramFileData = HistogramFileLocation+Sample+"/histograms_"+Sample+".root"
  fdatahists = TFile.Open(histogramFileData)
  if not fdatahists:
    print 'file:',histogramFileData,'not found; quitting'
    return
  #print "Getting data histogram"
  histosdata = fdatahists.Get("h_Diphoton_Minv_FineBinning")
  #print "histo",histosdata.GetName(),histosdata.GetEntries(),"entries",histosdata.Integral(),"(integral)"
  minBin = peakBin-nBinsUsedForWindow[indexMaxSsb][0]
  maxBin = peakBin+nBinsUsedForWindow[indexMaxSsb][1]
  entriesData = histosdata.Integral(minBin,maxBin)
  entGamJet = histosGammaJet.Integral(minBin,maxBin)
  entJetJet = histosJetJet.Integral(minBin,maxBin)
  entMC = histosmc.Integral(minBin,maxBin)
  entryBG = entGamJet+entJetJet+entMC
  # Fill the model point
  modelPoint.nDataObs = entriesData
  modelPoint.nBackground = entryBG
  modelPoint.nBackgroundErr = math.sqrt(entryBG) 
  modelPoint.totalEff = 1.0*signalHistogram.Integral(minBin,maxBin)/signalEntriesTotal

  massRangeLow = signalHistogram.GetBinLowEdge(minBin)
  massRangeHigh = signalHistogram.GetBinLowEdge(maxBin)+signalHistogram.GetBinWidth(maxBin)
  #return nBinsForMaxSOverSqrtB, maxSOverSqrtB, nBinsUsedForWindow, sOverSqrtBForWindow
  return massRangeLow, massRangeHigh, indexMaxSsb, nBinsUsedForWindow, ssbForWindow


def GetMinMaxMassArrays(modelPointArray, numsigmas = 3):
  # make list of lower/upper edges of mass windows
  minMasses = [(modelPoint.mass - numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  maxMasses = [(modelPoint.mass + numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  return minMasses,maxMasses


# import of yield-calculating code from ExoDiPhotonAnalyzer/test/PlottingCode/CreateHistogramFiles.C
def CalculateYieldsForMassRanges(HistogramFileLocation, modelPointArray, lumi, numsigmas = 3):
  gStyle.SetOptStat("ourmei")
  #Sample = "ExoDiPhotonAnalyzer_DataABC"
  Sample = "ExoDiPhotonAnalyzer_PFDec14th_DataABCD"

  # FIXME hardcoded names/locations
  # background
  histogramFileData = HistogramFileLocation+Sample+"/histograms_"+Sample+".root"
  histogramFileJetJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root"
  histogramFileGammaJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root"	
  histogramFileMC = HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root"
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
    mp.nBackgroundErr = math.sqrt(entryBG) 
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

  PrintEntries('Signal couplings',signalCouplings)
  PrintEntries('Signal Masses',signalMasses)
  PrintEntries('signalEntriesInWindow',entriesSignalInWindow)
  PrintEntries('signalEntriesTotal',entriesSignalTotal)
  signalEffs = [sigEvtsWindow/sigEvtsTotal for sigEvtsWindow,sigEvtsTotal in itertools.izip(entriesSignalInWindow,entriesSignalTotal)]
  PrintEntries('Signal efficiencies*acceptances',signalEffs)


