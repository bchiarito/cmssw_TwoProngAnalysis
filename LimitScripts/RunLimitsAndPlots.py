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




def DoLimitsAllPoints(cl95MacroPathAndName,lumi,lumiErr,limitsFileNameBase):
  print 'Run for coupling 0.01'
  with open(limitsFileNameBase+'0p01.txt', 'w') as file:
    ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,modelPointsC0p01,file)
  print 'Run for coupling 0.05'
  with open(limitsFileNameBase+'0p05.txt', 'w') as file:
    ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,modelPointsC0p05,file)
  print 'Run for coupling 0.1'
  with open(limitsFileNameBase+'0p1.txt', 'w') as file:
    ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,modelPointsC0p1,file)


def DoPlotsAllPoints(lumi):
  print 'Run for coupling 0.01'
  with open(limitsFileNameBase+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  PlotBands(readModelPoints0p01)
  m0p01,mExp0p01 = GetMassLimit(readModelPoints0p01)
  # 0.05
  print 'Run for coupling 0.05'
  with open(limitsFileNameBase+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  PlotBands(readModelPoints0p05)
  m0p05,mExp0p05 = GetMassLimit(readModelPoints0p05)
  # 0.1
  print 'Run for coupling 0.1'
  with open(limitsFileNameBase+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  PlotBands(readModelPoints0p1)
  m0p1,mExp0p1 = GetMassLimit(readModelPoints0p1)
  # 2-D results plot
  mLimObs = [m0p01,m0p05,m0p1]
  mLimExp = [mExp0p01,mExp0p05,mExp0p1]
  couplings = [0.01,0.05,0.1]
  #for i in range(0,len(mLimExp)):
  #  print 'coupling[',i,']:',couplings[i]
  #  print 'mLimExp[',i,']:',mLimExp[i]
  #  print 'mLimObs[',i,']:',mLimObs[i]
  CouplingVsMassPlot(couplings,mLimExp,mLimObs)
  # make table
  tableTitleString=string.ljust('Coupling',9)+string.ljust('Mass',6)+string.ljust('Eff.',8)
  tableTitleString+=string.center('Exp. Sig.',9)
  tableTitleString+=string.center('Exp. Bkg.',17)+string.center('Obs.',10)
  tableTitleString+=string.center('Th. XSec',10)+string.center('Exp. Lim.',12)
  tableTitleString+=string.center('Obs. Lim.',10)
  print
  print tableTitleString
  for modelPoint in readModelPoints0p01+readModelPoints0p05+readModelPoints0p1:
    print modelPoint.StringTableLine(lumi)
  


def Usage():
  print 'Usage: python RunLimitsAndPlots.py [arg] where arg can be:'
  print '    all --> run limits and plots (see below)'
  print '    limits --> Compute limits for all model points & write out results'
  print '    plots  --> Read limit results from text files and make limit plots'
  print '    yields --> Calculate event yields from histograms in root file'

# set up tdrStyle
gROOT.Reset()
gROOT.ProcessLine('.L tdrStyle.C')
gROOT.ProcessLine('setTDRStyle()')
# Define RooStats macro path and name
cl95MacroPath = os.environ['CMSSW_BASE']+'/src/StatisticalTools/RooStatsRoutines/root/'
cl95MacroName = 'roostats_cl95.C'
# location of root files from CreateHistogramFiles code
rootFileLocation = '/afs/cern.ch/user/s/scooper/work/private/results/diPhotonHistograms/'
# limit results file base name
limitsFileNameBase = 'results_limits_k_'
# Declarations of Lumi and model points
lumi = 10252.
lumiErr = lumi*0.044
# No Cuts
# coupling, mass, totalXSec, totalEff, width, nDataObs, nBG, nBGerr
##c = 0.01
#modelPointsC0p01 = []
#modelPointsC0p01.append(ModelPoint(0.01, 750,  1.023e-02, 30.1475/100, 5.2041,  4, 3.47787, 1.86490))
#modelPointsC0p01.append(ModelPoint(0.01, 1000, 2.072e-03, 34.7347/100, 6.50326, 0, 0.46004, 0.67826))
#modelPointsC0p01.append(ModelPoint(0.01, 1250, 5.21e-04,  38.7874/100, 8.53237, 0, 0.21650, 0.46530))
#modelPointsC0p01.append(ModelPoint(0.01, 1500, 1.604e-04, 42.5412/100, 10.0917, 0, 0.11362, 0.33708))
#modelPointsC0p01.append(ModelPoint(0.01, 1750, 5.408e-05, 46.6982/100, 11.5976, 0, 0.04732, 0.21752))
#modelPointsC0p01.append(ModelPoint(0.01, 2000, 1.853e-05, 49.1978/100, 13.6845, 0, 0.02484, 0.15760))
#modelPointsC0p01.append(ModelPoint(0.01, 3000, 4.703e-07, 54.5872/100, 22.4852, 0, 0.00193, 0.04389))
##c = 0.05
#modelPointsC0p05 = []
#modelPointsC0p05.append(ModelPoint(0.05, 1750, 1.331e-03, 46.1352/100, 13.7476, 0, 0.04732, 0.21752))
#modelPointsC0p05.append(ModelPoint(0.05, 2000, 4.665e-04, 49.1199/100, 16.795,  0, 0.02484, 0.15760))
#modelPointsC0p05.append(ModelPoint(0.05, 2500, 7.226e-05, 52.9777/100, 20.4539, 0, 0.00456, 0.06750))
#modelPointsC0p05.append(ModelPoint(0.05, 2750, 2.803e-05, 53.6758/100, 22.7374, 0, 0.00182, 0.04265))
#modelPointsC0p05.append(ModelPoint(0.05, 3000, 1.169e-05, 53.6002/100, 25.2982, 0, 0.00193, 0.04389))
##c = 0.1
#modelPointsC0p1 = []
#modelPointsC0p1.append(ModelPoint(0.1, 2250, 7.04e-04, 51.4743/100, 26.3803, 0, 0.01617, 0.12718))
#modelPointsC0p1.append(ModelPoint(0.1, 2500, 2.79e-04, 52.9381/100, 30.8038, 0, 0.00674, 0.08211))
#modelPointsC0p1.append(ModelPoint(0.1, 2750, 1.14e-04, 53.6588/100, 35.9283, 0, 0.00321, 0.05667))
#modelPointsC0p1.append(ModelPoint(0.1, 3000, 4.68e-05, 53.512/100,  38.7399, 0, 0.00193, 0.04390))
#modelPointsC0p1.append(ModelPoint(0.1, 3250, 1.19e-05, 53.5149/100, 41.8178, 0, 0.00326, 0.05711))
#modelPointsC0p1.append(ModelPoint(0.1, 3500, 7.7e-05,  52.9363/100, 40.2991, 0, 0.00326, 0.05711))

# With DeltaPhi 2.8 cut, but old efficiencies
# coupling, mass, totalXSec, totalEff, width, nDataObs, nBG, nBGerr
#c = 0.01
modelPointsC0p01 = []
modelPointsC0p01.append(ModelPoint(0.01, 750,  1.023e-02, 30.1475/100, 5.2041,  4, 3.29455, 1.81509))
modelPointsC0p01.append(ModelPoint(0.01, 1000, 2.072e-03, 34.7347/100, 6.50326, 0, 0.44277, 0.66541))
modelPointsC0p01.append(ModelPoint(0.01, 1250, 5.21e-04,  38.7874/100, 8.53237, 0, 0.21078, 0.45911))
modelPointsC0p01.append(ModelPoint(0.01, 1500, 1.604e-04, 42.5412/100, 10.0917, 0, 0.10940, 0.33076))
modelPointsC0p01.append(ModelPoint(0.01, 1750, 5.408e-05, 46.6982/100, 11.5976, 0, 0.04682, 0.21638))
modelPointsC0p01.append(ModelPoint(0.01, 2000, 1.853e-05, 49.1978/100, 13.6845, 0, 0.02463, 0.15695))
modelPointsC0p01.append(ModelPoint(0.01, 3000, 4.703e-07, 54.5872/100, 22.4852, 0, 0.00193, 0.04389))
#c = 0.05
modelPointsC0p05 = []
modelPointsC0p05.append(ModelPoint(0.05, 1750, 1.331e-03, 46.1352/100, 13.7476, 0, 0.04682, 0.21638))
modelPointsC0p05.append(ModelPoint(0.05, 2000, 4.665e-04, 49.1199/100, 16.795,  0, 0.02463, 0.15695))
modelPointsC0p05.append(ModelPoint(0.05, 2500, 7.226e-05, 52.9777/100, 20.4539, 0, 0.00445, 0.06668))
modelPointsC0p05.append(ModelPoint(0.05, 2750, 2.803e-05, 53.6758/100, 22.7374, 0, 0.00182, 0.04265))
modelPointsC0p05.append(ModelPoint(0.05, 3000, 1.169e-05, 53.6002/100, 25.2982, 0, 0.00193, 0.04389))
#c = 0.1
modelPointsC0p1 = []
modelPointsC0p1.append(ModelPoint(0.1, 2250, 7.04e-04, 51.4743/100, 26.3803, 0, 0.01589, 0.12605))
modelPointsC0p1.append(ModelPoint(0.1, 2500, 2.79e-04, 52.9381/100, 30.8038, 0, 0.00663, 0.08143))
modelPointsC0p1.append(ModelPoint(0.1, 2750, 1.14e-04, 53.6588/100, 35.9283, 0, 0.00321, 0.05667))
modelPointsC0p1.append(ModelPoint(0.1, 3000, 4.68e-05, 53.512/100,  38.7399, 0, 0.00193, 0.04390))
modelPointsC0p1.append(ModelPoint(0.1, 3250, 1.19e-05, 53.5149/100, 41.8178, 0, 0.00326, 0.05711))
modelPointsC0p1.append(ModelPoint(0.1, 3500, 7.7e-05,  52.9363/100, 40.2991, 0, 0.00326, 0.05711))
#
#
# RUN
#
# Parse arguments
if len(sys.argv)==1:
  Usage()
  sys.exit()
elif sys.argv[1]=='yields':
  print 'yields: MakeYieldsTableForMassRanges'
  MakeYieldsTableForMassRanges(rootFileLocation, modelPointsC0p01, lumi, 3)
elif sys.argv[1]=='limits':
  print 'limits: DoLimitsAllPoints'
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileNameBase)
elif sys.argv[1]=='plots':
  print 'plots: DoPlotsAllPoints'
  DoPlotsAllPoints(lumi)
# TODO: add yield calculation into all
# TODO: calculate signal efficiency from here
elif sys.argv[1]=='all':
  print 'all: limits+plots'
  print 'limits: DoLimitsAllPoints'
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileNameBase)
  print 'plots: DoPlotsAllPoints'
  DoPlotsAllPoints(lumi)
else:
  print 'Did not understand input.'
  Usage()
  sys.exit()


### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#  rep = ''
#  while not rep in [ 'q', 'Q' ]:
#    rep = raw_input( 'enter "q" to quit: ' )
#    if 1 < len(rep):
#      rep = rep[0]






