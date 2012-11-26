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

from ROOT import *
from ModelPoint import *


def ComputeLimits(cl95MacroPath,lumi,lumiErr,modelPointArray,file):
  gROOT.ProcessLine('.L '+cl95MacroPath+'+')
  ROOT.SetParameter("Optimize",False)
  ROOT.SetParameter("MakePlot",True)
  ROOT.SetParameter("WriteResult",True)
  ROOT.SetParameter("PlotHypoTestResult",True)
  for modelPoint in modelPointArray:
    # new interface
    #limit = ROOT.GetClsLimit(lumi, lumiErr, modelPoint.GetTotalEff(), 
    #                         math.sqrt((modelPoint.GetTotalEff()*(1.-modelPoint.GetTotalEff()))/25000.),
    #                         modelPoint.GetNBackground(), modelPoint.GetNBackgroundErr(),
    #                         modelPoint.GetNDataObs());
  
    # legacy interface
    limit = ROOT.roostats_limit(lumi, lumiErr, modelPoint.totalEff,
                                       math.sqrt((modelPoint.totalEff*(1.-modelPoint.totalEff))/25000.),
                                       modelPoint.nBackground, modelPoint.nBackgroundErr,
                                       modelPoint.nDataObs, False, 0, "cls",
                                       "plots_cl95_"+str(modelPoint.coupling)+"_"+str(modelPoint.mass)+".pdf",0)
                                       #"",0)
                                       #"",123456) # seed for testing
    couplingStr = str(modelPoint.coupling).replace('.','p')
    plotFileSamplingName = 'cls_sampling_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    plotFileName = 'cls_k_'+couplingStr+"_m"+str(modelPoint.mass)+".pdf"
    if os.path.isdir('cls_plots') == False:
      subprocess.Popen(['mkdir','cls_plots'])
    subprocess.Popen(['mv','cls_sampling.pdf','cls_plots/'+plotFileSamplingName])
    subprocess.Popen(['mv','cls.pdf','cls_plots/'+plotFileName])
    modelPoint.AddLimitResult(limit);
    modelPoint.Write(file)
  print 'Used StatisticalTools/RooStatsRoutine CVS tag:',GetRooStatsMacroCVSTag(cl95MacroPath)


def GetRooStatsMacroCVSTag(cl95MacroPath):
  cl95Split = cl95MacroPath.split('/')
  cl95MacroName = cl95Split[len(cl95Split)-1]
  cl95MacroDir = cl95MacroPath.rstrip(cl95Split[len(cl95Split)-1])
  proc = subprocess.Popen(['cvs','status',cl95MacroName],cwd=cl95MacroDir,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
  out,err = proc.communicate()
  split = out.split()
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
    print "%.5f,"%math.sqrt(entry),
  print "}"
  print name,"syst = {",
  for entry in listErrors:
    print "%.5f,"%entry,
  print "}"


# import of yield-calculating code from ExoDiPhotonAnalyzer/test/PlottingCode/CreateHistogramFiles.C
def MakeYieldsTableForMassRanges(HistogramFileLocation, modelPointArray, lumi, numsigmas = 2):
  gStyle.SetOptStat("ourmei")
  Sample = "ExoDiPhotonAnalyzer_DataABC"

  # FIXME hardcoded names/locations
  histogramdata = HistogramFileLocation+Sample+"/histograms_"+Sample+".root"
  histogramJetJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_JetJet.root"
  histogramGammaJet = HistogramFileLocation+Sample+"/histograms_"+Sample+"_GammaJet.root"	
  histogramMC = HistogramFileLocation+"/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root"
  fdatahists = TFile.Open(histogramdata)
  fJetJethists = TFile.Open(histogramJetJet)
  fGammaJethists = TFile.Open(histogramGammaJet)
  fMChists = TFile.Open(histogramMC)

  print "Getting data histogram"
  histosdata = fdatahists.Get("h_Diphoton_Minv")
  print "histo",histosdata.GetName(),histosdata.GetEntries(),"entries",histosdata.Integral(),"(integral)"

  print "Getting MC histogram"
  histosmc = fMChists.Get("h_Diphoton_Minv")
  histosmc.Scale(lumi)
  print "histo",histosmc.GetName(),histosmc.GetEntries(),"entries",histosmc.Integral(),"(integral)"

  print "Getting Jet+Jet histograms"
  histosJetJet = fJetJethists.Get("h_JetJet_minv")
  histosJetJetUpperError = fJetJethists.Get("h_JetJet_minv_UpperError")
  histosJetJetLowerError = fJetJethists.Get("h_JetJet_minv_LowerError")
  print "histo",histosJetJet.GetName(),histosJetJet.GetEntries(),"entries",histosJetJet.Integral(),"(integral)"
  print "histo",histosJetJetUpperError.GetName(),histosJetJetUpperError.GetEntries(),"entries",histosJetJetUpperError.Integral(),"(integral)"
  print "histo",histosJetJetLowerError.GetName(),histosJetJetLowerError.GetEntries(),"entries",histosJetJetLowerError.Integral(),"(integral)"

  print "Getting Gamma+Jet histograms"
  histosGammaJet = fGammaJethists.Get("h_GammaJet_minv")
  histosGammaJetUpperError = fGammaJethists.Get("h_GammaJet_minv_UpperError")
  histosGammaJetLowerError = fGammaJethists.Get("h_GammaJet_minv_LowerError")
  print "histo",histosGammaJet.GetName(),histosGammaJet.GetEntries(),"entries",histosGammaJet.Integral(),"(integral)"
  print "histo",histosGammaJetUpperError.GetName(),histosGammaJetUpperError.GetEntries(),"entries",histosGammaJetUpperError.Integral(),"(integral)"
  print "histo",histosGammaJetLowerError.GetName(),histosGammaJetLowerError.GetEntries(),"entries",histosGammaJetLowerError.Integral(),"(integral)"

  ## make list of lower/upper edges of mass windows
  minMasses = [(modelPoint.mass - numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]
  maxMasses = [(modelPoint.mass + numsigmas * modelPoint.halfWidth) for modelPoint in modelPointArray]

  # calculate integrals in different mass ranges
  # They all have the same binning
  # Let's take the data one
  nbinsX = histosdata.GetNbinsX()
  binnr = nbinsX
  binminnumbers = []
  binmaxnumbers = []

  print "---------------Mass intervals to seek----------------"
# debug
#  print 'nBins',histosdata.GetNbinsX()
  for minMass,maxMass in itertools.izip(minMasses,maxMasses):
    print "interval [",minMass,",",maxMass,"]"
#    print 'minMass',minMass,'binLowEdge',histosdata.GetBinLowEdge(histosdata.FindBin(minMass)),'binNr',histosdata.FindBin(minMass)
    binminnumbers.append(histosdata.FindBin(minMass)+1)
    maxMassBin = histosdata.FindBin(maxMass)
    if maxMassBin == histosdata.GetNbinsX()+1: # is it the overflow bin?
      maxMassBin = histosdata.GetNbinsX()
#    print 'maxMass',maxMass,'binLowEdge',histosdata.GetBinLowEdge(maxMassBin),'binNr',maxMassBin
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

  entriesDataOverBackground = [entData/entBg for entData,entBg in zip(entriesData,entriesBackground)]
  errorsDataOverBackground = []
  for entData,entBg,errData,errBg in itertools.izip(entriesData,entriesBackground,errorsData,errorsBackground):
    if entData != 0:
      errorsDataOverBackground.append(entData/entBg * sqrt( (errData/entData)*(errData/entData) + (errBg/entBg)*(errBg/entBg) ))
    else:
      errorsDataOverBackground.append(-1)
  PrintInfo('Data/Background',entriesDataOverBackground,errorsDataOverBackground)



