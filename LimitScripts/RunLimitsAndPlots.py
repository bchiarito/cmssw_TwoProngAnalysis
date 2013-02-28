#!/usr/bin/env python

#
# Script to compute limits and make plots
# Runs support functions in MakePlots.py and MakeLimits.py
# Seth I. Cooper, U. Alabama
# November 2012
#
#   * Must set up CMSSW area and do cmsenv
#   * Must check out StatisticalTools/RooStatsRoutines package and then:
#     source $CMSSW_BASE/src/StatisticalTools/RooStatsRoutines/setup/setup/lxplus_standalone_setup.sh
#   * More instructions/notes can be found on the twiki page:
#     https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaDiphotonResonance2012#Limit_Setting_Code_Instructions


import string
import os
import sys
import glob
import math
import array
import subprocess

# run root in batch
sys.argv.append('-b')
#import ROOT
from ROOT import *
from MakePlots import *
from MakeLimits import *


# mapping of couplings/masses to halfWidths --> taken as input from fits
halfWidths0p01 = dict()
halfWidths0p01[750]  = 5.2041
halfWidths0p01[1000] = 6.50326
halfWidths0p01[1250] = 8.53237
halfWidths0p01[1500] = 10.0917
halfWidths0p01[1750] = 11.5976
halfWidths0p01[2000] = 13.6845
halfWidths0p01[3000] = 22.4852
# 0.05
halfWidths0p05 = dict()
halfWidths0p05[1750] = 13.7476
halfWidths0p05[2000] = 16.795
halfWidths0p05[2500] = 20.4539
halfWidths0p05[2750] = 22.7374
halfWidths0p05[3000] = 25.2982
# 0.1
halfWidths0p1 = dict()
halfWidths0p1[2250] = 26.3803
halfWidths0p1[2500] = 30.8038
halfWidths0p1[2750] = 35.9283
halfWidths0p1[3000] = 38.7399
halfWidths0p1[3250] = 41.8178
halfWidths0p1[3500] = 40.2991
# mapping of couplings/masses to crossSections --> taken as input from theoretical calc.
totalXSecs0p01 = dict()
totalXSecs0p01[750] = 1.023e-02
totalXSecs0p01[1000] = 2.072e-03
totalXSecs0p01[1250] = 5.21e-04
totalXSecs0p01[1500] = 1.604e-04
totalXSecs0p01[1750] = 5.408e-05
totalXSecs0p01[2000] = 1.853e-05
totalXSecs0p01[3000] = 4.703e-07
# 0.05
totalXSecs0p05 = dict()
totalXSecs0p05[1750] = 1.331e-03
totalXSecs0p05[2000] = 4.665e-04
totalXSecs0p05[2500] = 7.226e-05
totalXSecs0p05[2750] = 2.803e-05
totalXSecs0p05[3000] = 1.169e-05
# 0.1
totalXSecs0p1 = dict()
totalXSecs0p1[2250] = 7.04e-04
totalXSecs0p1[2500] = 2.79e-04
totalXSecs0p1[2750] = 1.14e-04
totalXSecs0p1[3000] = 4.68e-05
totalXSecs0p1[3250] = 1.19e-05
totalXSecs0p1[3500] = 7.7e-06


def GetHalfWidth(coupling,mass):
  if coupling==0.01:
    dict = halfWidths0p01
  elif coupling==0.05:
    dict = halfWidths0p05
  elif coupling==0.1:
    dict = halfWidths0p1
  else:
    dict = dict()
  if mass in dict:
    return dict[mass]
  else:
    print 'GetHalfWidth: Coupling',coupling,'mass',mass,'not recognized; quitting'
    exit()


def GetXSec(coupling,mass):
  if coupling==0.01:
    dict = totalXSecs0p01
  elif coupling==0.05:
    dict = totalXSecs0p05
  elif coupling==0.1:
    dict = totalXSecs0p1
  else:
    dict = dict()
  if mass in dict:
    return dict[mass]
  else:
    print 'GetXSec: Coupling',coupling,'mass',mass,'not recognized; quitting'
    exit()


def GetSignalFileName(template,coupling,mass):
  if coupling==0.01:
    couplingStr = '001'
  elif coupling==0.05:
    couplingStr = '005'
  elif coupling==0.1:
    couplingStr = '01'
  else:
    print 'GetSignalFileName: Coupling',coupling,'not recognized; quitting'
    exit()
  massStr = str(mass)
  return template.format(coupling=couplingStr,mass=massStr)


# TODO: remove hardcoding; instead use list of couplings to run over
def DoLimitsAllPoints(cl95MacroPathAndName,lumi,lumiErr,limitsFileNameBase):
  print 'Run for coupling 0.01'
  with open(limitsFileNameBase+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p01,limitsFileNameBase+'0p01.txt')
  print 'Run for coupling 0.05'
  with open(limitsFileNameBase+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p05,limitsFileNameBase+'0p05.txt')
  print 'Run for coupling 0.1'
  with open(limitsFileNameBase+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p1,limitsFileNameBase+'0p1.txt')


## TODO
# function to read in the coupling we ran over
#def ReadResultsAllPoints():
#  resultsByCoupling = {}
#  for c in couplingsToRunOver:
#    # make coupling string
#    couplingString = f(c)
#    with open(limitsFileNameBase+couplingString+'.txt', 'r') as file:
#      readModelPoints = ReadFromFile(file)
#    resultsByCoupling[c] = readModelPoints
#  return resultsByCoupling


def DoPlotsAllPoints(lumi,rootFile,pathToTDRStyle):
  # set up tdrStyle
  gROOT.ProcessLine('.L '+pathToTDRStyle)
  gROOT.ProcessLine('setTDRStyle()')
  print 'Run for coupling 0.01'
  with open(limitsFileNameBase+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  for modelPoint in readModelPoints0p01:
    modelPoint.Print()
  PlotBands(readModelPoints0p01,lumi,rootFile)
  m0p01,xs,mExp0p01,xsE = GetMassLimit(readModelPoints0p01)
  print string.ljust('Coupling: '+str(readModelPoints0p01[0].coupling),14),
  print ' Observed limit mass: %0.2f'%m0p01
  print '                Expected limit mass: %0.2f'%mExp0p01
  print '                Observed XSec limit: %0.6f'%xs
  print '                Expected XSec limit: %0.6f'%xsE
  # 0.05
  print 'Run for coupling 0.05'
  with open(limitsFileNameBase+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  PlotBands(readModelPoints0p05,lumi,rootFile)
  m0p05,xs,mExp0p05,xsE = GetMassLimit(readModelPoints0p05)
  print string.ljust('Coupling: '+str(readModelPoints0p05[0].coupling),14),
  print ' Observed limit mass: %0.2f'%m0p05
  print '                Expected limit mass: %0.2f'%mExp0p05
  print '                Observed XSec limit: %0.6f'%xs
  print '                Expected XSec limit: %0.6f'%xsE
  # 0.1
  print 'Run for coupling 0.1'
  with open(limitsFileNameBase+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  PlotBands(readModelPoints0p1,lumi,rootFile)
  m0p1,xs,mExp0p1,xsE = GetMassLimit(readModelPoints0p1)
  print string.ljust('Coupling: '+str(readModelPoints0p1[0].coupling),14),
  print ' Observed limit mass: %0.2f'%m0p1
  print '                Expected limit mass: %0.2f'%mExp0p1
  print '                Observed XSec limit: %0.6f'%xs
  print '                Expected XSec limit: %0.6f'%xsE
  # 2-D results plot
  mLimObs = [m0p01,m0p05,m0p1]
  mLimExp = [mExp0p01,mExp0p05,mExp0p1]
  couplings = [0.01,0.05,0.1]
  CouplingVsMassPlot(couplings,mLimExp,mLimObs,rootFile,lumi)
  # efficiencies for all couplings/masses on same axes
  PlotAllEfficiencies([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile)
  # half widths
  PlotAllHalfWidths([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile)
  # exp BG
  PlotAllExpBGs([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile)
  # limit plot for all couplings on same axes
  PlotAllBands([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile)
  # make table
  DoTablesAllPoints(lumi)


def DoTablesAllPoints(lumi):
  with open(limitsFileNameBase+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  # 0.05
  with open(limitsFileNameBase+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  # 0.1
  with open(limitsFileNameBase+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  # now make table output
  tableTitleString=string.ljust('Coupling',9)+string.ljust('Mass',6)+string.ljust('Eff.',8)
  tableTitleString+=string.center('Exp. Sig.',9)
  tableTitleString+=string.center('Exp. Bkg.',17)+string.center('Obs.',10)
  tableTitleString+=string.center('Th. XSec',10)+string.center('Exp. Lim.',12)
  tableTitleString+=string.center('Obs. Lim.',10)
  print
  print tableTitleString
  for modelPoint in readModelPoints0p01+readModelPoints0p05+readModelPoints0p1:
    print modelPoint.StringTableLine(lumi)
  # latex table
  # k, mass, windowRange, sigEff, expSigEvts, expBgEvts, obs
  print
  print
  print '\\begin{table}[htpb]\n\t\\begin{center}'
  print '\t\t\\begin{tabular}{ccccccc}\n\t\t\\hline'
  print '\t\t$\\tilde{k}$ & $M_1$ & Window & Sig. Eff. & Exp. Sig. Evts. & Exp. Bg. Evts. & Obs. \\\\'
  print '\t\t\\hline'
  for modelPoint in readModelPoints0p01+readModelPoints0p05+readModelPoints0p1:
    print modelPoint.LatexTableLine(lumi)
  print '\t\t\\hline'
  print '\t\t\\end{tabular}'
  captionLine="\t\t\\caption[Event Yields]{Event yields of signal and data after selection.  "
  captionLine+="The columns show: coupling ``$\\tilde{k}$'', graviton mass ``$M_1$'' in GeV, mass window range ``Window'' in GeV, "
  captionLine+="signal efficiency ``Sig. Eff.'', expected number of signal events ``Exp. Sig. Evts.'', "
  captionLine+="expected number of background events and error ``Exp. Bg. Evts.'', and observed data events ``Obs.''.}"
  print captionLine
  print '\t\\label{table:eventYields}'
  print '\t\\end{center}'
  print '\\end{table}'
  # next table -- limits
  print
  print
  print '\\begin{table}[htpb]\n\t\\begin{center}'
  print '\t\t\\begin{tabular}{ccc}\n\t\t\\hline'
  print '\t\t\t$\\tilde{k}$ & $\sigma$ & $M_1$ \\\\'
  print '\t\t\t\\hline'
  # FIXME remove hardcoding of couplings to run over
  # 0.01
  latexLine='\t\t\t'+'0.01'+' & '
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p01)
  latexLine+='%.5f'%xsLimObs+' & '
  latexLine+=str(int(massLimObs))
  latexLine+=' \\\\'
  print latexLine
  # 0.05
  latexLine='\t\t\t'+'0.05'+' & '
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p05)
  latexLine+='%.5f'%xsLimObs+' & '
  latexLine+=str(int(massLimObs))
  latexLine+=' \\\\'
  print latexLine
  # 0.1
  latexLine='\t\t\t'+'0.1'+' & '
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p1)
  latexLine+='%.5f'%xsLimObs+' & '
  latexLine+=str(int(massLimObs))
  latexLine+=' \\\\'
  print latexLine
  print '\t\t\t\\hline'
  print '\t\t\\end{tabular}'
  captionLine="\t\t\\caption[Limit results]{ Limit results from the observed data. "
  captionLine+="The columns show: coupling ``$\\tilde{k}$'', observed 95\% confidence level upper limit on the graviton cross section in pb, "
  captionLine+="and the observed 95\% confidence level lower limit on the graviton mass ``$M_1$'' in GeV. }"
  print captionLine
  print '\t\\label{table:limitResults}'
  print '\t\\end{center}'
  print '\\end{table}'
  

def DoOptimizeAllPoints():
  maxWindowRange = 600 # bins/GeV
  useAsymmWindow = True
  print 'Run for coupling 0.01'
  colorIndex = 2 #TODO add this into modelpoint itself?
  with open(limitsFileNameBase+'0p01.txt', 'w') as file:
    graphOptHalfWindowsVsMass0p01,graphOptMinMaxWindowsVsMass0p01 = OptimizeSignalMassWindows(
        rootFileLocation,modelPointsC0p01,lumi,useAsymmWindow,maxWindowRange,file,rootFile,colorIndex)
  print 'Run for coupling 0.05'
  colorIndex = 4
  with open(limitsFileNameBase+'0p05.txt', 'w') as file:
    graphOptHalfWindowsVsMass0p05,graphOptMinMaxWindowsVsMass0p05 = OptimizeSignalMassWindows(
        rootFileLocation,modelPointsC0p05,lumi,useAsymmWindow,maxWindowRange,file,rootFile,colorIndex)
  print 'Run for coupling 0.1'
  colorIndex = 8
  with open(limitsFileNameBase+'0p1.txt', 'w') as file:
    graphOptHalfWindowsVsMass0p1,graphOptMinMaxWindowsVsMass0p1 = OptimizeSignalMassWindows(
        rootFileLocation,modelPointsC0p1,lumi,useAsymmWindow,maxWindowRange,file,rootFile,colorIndex)
  # make multigraphs for all masses/couplings
  MakeOptHalfWindowVsMassMultigraph(graphOptHalfWindowsVsMass0p01,graphOptHalfWindowsVsMass0p05,graphOptHalfWindowsVsMass0p1,rootFile)
  MakeOptMassWindowsVsMassMultiGraph(graphOptMinMaxWindowsVsMass0p01,graphOptMinMaxWindowsVsMass0p05,graphOptMinMaxWindowsVsMass0p1,rootFile)


def Usage():
  print 'Usage: python RunLimitsAndPlots.py [arg] where arg can be:'
  print '    all      --> run yields, limits, and plots (see below)'
  print '    limits   --> Compute limits for all model points & write out results'
  print '    plots    --> Read limit results from text files and make limit plots'
  print '    tables   --> Read limit results from text files and make results tables (text/latex)'
  print '    yields   --> Calculate event yields from histograms in root files (from CreateHistogramFiles)'
  print '    optimize --> Calculate optimal inv. mass windows using histograms in root files (from CreateHistogramFiles)'




#
#
# RUN
#
gROOT.Reset()
# get path to tdrStyle script (should always be in LimitScripts, but we don't hardcode it anyway)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
pathToTDRStyle = dname+'/tdrStyle.C'
# Define RooStats macro path and name
cl95MacroPath = os.environ['CMSSW_BASE']+'/src/StatisticalTools/RooStatsRoutines/root/'
cl95MacroName = 'roostats_cl95.C'


# Configurable stuff here
# limit results file base name
limitsFileNameBase = 'results_limits_k_'
# location of backgroundMC/data root files from CreateHistogramFiles code
rootFileLocation = '/afs/cern.ch/user/s/scooper/work/private/results/diPhotonHistogramsPF_deltaPhi2p8_19p6invFb/'
# Declarations of Lumi and model points to consider
lumi = 19620.
lumiErr = lumi*0.044
masses0p01 = [750,1000,1250,1500,1750,2000,3000]
masses0p05 = [1750,2000,2500,2750,3000]
masses0p1 = [2250,2500,2750,3000,3250,3500]
# root file for plots
plotFileName = 'plots.root'
rootFile = TFile(plotFileName,'recreate')


# List signal histogram file template; rest of quantities are filled from functions (xsec, width, etc.) or histograms in the files
signalHistogramFilesPathTemplate = rootFileLocation+'diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia_merged/histograms_diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia_merged.root'
# for now, we use the directory structure that the CreateHistogramFiles.C code uses
# setup k=0.01
modelPointsC0p01 = []
for mass in masses0p01:
  k = 0.01
  modelPointsC0p01.append(ModelPoint(coupling=k,mass=mass,totalXSec=GetXSec(k,mass),halfWidth=GetHalfWidth(k,mass),
    fileName=GetSignalFileName(signalHistogramFilesPathTemplate,k,mass))) 
# setup k=0.05
modelPointsC0p05 = []
for mass in masses0p05:
  k = 0.05
  modelPointsC0p05.append(ModelPoint(coupling=k,mass=mass,totalXSec=GetXSec(k,mass),halfWidth=GetHalfWidth(k,mass),
    fileName=GetSignalFileName(signalHistogramFilesPathTemplate,k,mass))) 
# setup k=0.1
modelPointsC0p1 = []
for mass in masses0p1:
  k = 0.1
  modelPointsC0p1.append(ModelPoint(coupling=k,mass=mass,totalXSec=GetXSec(k,mass),halfWidth=GetHalfWidth(k,mass),
    fileName=GetSignalFileName(signalHistogramFilesPathTemplate,k,mass))) 

# Parse arguments
if len(sys.argv)==1:
  Usage()
  sys.exit()
elif sys.argv[1]=='optimize':
  print 'optimize: OptimizeSignalMassWindows'
  DoOptimizeAllPoints()
  PlotAllEfficiencies([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile)
elif sys.argv[1]=='yields':
  print 'yields: CalculateYieldsForMassRanges'
  CalculateYieldsForMassRanges(rootFileLocation, modelPointsC0p01+modelPointsC0p05+modelPointsC0p1, lumi, 3)
elif sys.argv[1]=='limits':
  print 'limits: DoLimitsAllPoints'
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileNameBase)
elif sys.argv[1]=='plots':
  print 'plots: DoPlotsAllPoints'
  DoPlotsAllPoints(lumi,rootFile,pathToTDRStyle)
elif sys.argv[1]=='tables':
  print 'tables: DoTablesAllPoints'
  DoTablesAllPoints(lumi)
elif sys.argv[1]=='all':
  print 'all: yields+limits+plots'
  print 'yields: CalculateYieldsForMassRanges'
  CalculateYieldsForMassRanges(rootFileLocation, modelPointsC0p01+modelPointsC0p05+modelPointsC0p1, lumi, 3)
  print 'limits: DoLimitsAllPoints'
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileNameBase)
  print 'plots: DoPlotsAllPoints'
  DoPlotsAllPoints(lumi,rootFile,pathToTDRStyle)
  print 'tables: DoTablesAllPoints'
  DoTablesAllPoints(lumi)
else:
  print 'Did not understand input.'
  Usage()
  sys.exit()

rootFile.Close()


### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#  rep = ''
#  while not rep in [ 'q', 'Q' ]:
#    rep = raw_input( 'enter "q" to quit: ' )
#    if 1 < len(rep):
#      rep = rep[0]






