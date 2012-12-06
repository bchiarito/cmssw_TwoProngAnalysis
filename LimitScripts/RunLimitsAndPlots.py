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




# TODO: remove hardcoding; instead use list of couplings to run over
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


def DoPlotsAllPoints(lumi):
  print 'Run for coupling 0.01'
  with open(limitsFileNameBase+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  PlotBands(readModelPoints0p01,lumi)
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
  PlotBands(readModelPoints0p05,lumi)
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
  PlotBands(readModelPoints0p1,lumi)
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
  CouplingVsMassPlot(couplings,mLimExp,mLimObs)
  # limit plot for all couplings on same axes
  PlotAllBands([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi)
  # make table
  #DoTablesAllPoints(lumi,3)


def DoTablesAllPoints(lumi,numsigmas):
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
    print modelPoint.LatexTableLine(lumi,numsigmas) # numsigmas for mass windows
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
# No Cuts -- updated signal eff.
# coupling, mass, totalXSec, totalEff, width, nDataObs, nBG, nBGerr
##c = 0.01
#modelPointsC0p01 = []
#modelPointsC0p01.append(ModelPoint(0.01, 750,  1.023e-02, 0.300598, 5.2041,  4, 3.47787, 1.86490))
#modelPointsC0p01.append(ModelPoint(0.01, 1000, 2.072e-03, 0.346988, 6.50326, 0, 0.46004, 0.67826))
#modelPointsC0p01.append(ModelPoint(0.01, 1250, 5.21e-04,  0.387196, 8.53237, 0, 0.21650, 0.46530))
#modelPointsC0p01.append(ModelPoint(0.01, 1500, 1.604e-04, 0.424459, 10.0917, 0, 0.11362, 0.33708))
#modelPointsC0p01.append(ModelPoint(0.01, 1750, 5.408e-05, 0.466624, 11.5976, 0, 0.04732, 0.21752))
#modelPointsC0p01.append(ModelPoint(0.01, 2000, 1.853e-05, 0.491459, 13.6845, 0, 0.02484, 0.15760))
#modelPointsC0p01.append(ModelPoint(0.01, 3000, 4.703e-07, 0.545592, 22.4852, 0, 0.00193, 0.04389))
##c = 0.05
#modelPointsC0p05 = []
#modelPointsC0p05.append(ModelPoint(0.05, 1750, 1.331e-03, 0.46056,  13.7476, 0, 0.04732, 0.21752))
#modelPointsC0p05.append(ModelPoint(0.05, 2000, 4.665e-04, 0.490243, 16.795,  0, 0.02484, 0.15760))
#modelPointsC0p05.append(ModelPoint(0.05, 2500, 7.226e-05, 0.528822, 20.4539, 0, 0.00456, 0.06750))
#modelPointsC0p05.append(ModelPoint(0.05, 2750, 2.803e-05, 0.5364,   22.7374, 0, 0.00182, 0.04265))
#modelPointsC0p05.append(ModelPoint(0.05, 3000, 1.169e-05, 0.535603, 25.2982, 0, 0.00193, 0.04389))
##c = 0.1
#modelPointsC0p1 = []
#modelPointsC0p1.append(ModelPoint(0.1, 2250, 7.04e-04, 0.514227, 26.3803, 0, 0.01617, 0.12718))
#modelPointsC0p1.append(ModelPoint(0.1, 2500, 2.79e-04, 0.528942, 30.8038, 0, 0.00674, 0.08211))
#modelPointsC0p1.append(ModelPoint(0.1, 2750, 1.14e-04, 0.536109, 35.9283, 0, 0.00321, 0.05667))
#modelPointsC0p1.append(ModelPoint(0.1, 3000, 4.68e-05, 0.53504,  38.7399, 0, 0.00193, 0.04390))
#modelPointsC0p1.append(ModelPoint(0.1, 3250, 1.19e-05, 0.534551, 41.8178, 0, 0.00326, 0.05711))
#modelPointsC0p1.append(ModelPoint(0.1, 3500, 7.7e-05,  0.529164, 40.2991, 0, 0.00326, 0.05711))

# With DeltaPhi 2.8 cut, updated efficiencies
# coupling, mass, totalXSec, totalEff, width, nDataObs, nBG, nBGerr
#c = 0.01
modelPointsC0p01 = []
modelPointsC0p01.append(ModelPoint(0.01, 750,  1.023e-02, 0.245853, 5.2041,  4, 3.29455, 1.81509))
modelPointsC0p01.append(ModelPoint(0.01, 1000, 2.072e-03, 0.302633, 6.50326, 0, 0.44277, 0.66541))
modelPointsC0p01.append(ModelPoint(0.01, 1250, 5.21e-04,  0.349581, 8.53237, 0, 0.21078, 0.45911))
modelPointsC0p01.append(ModelPoint(0.01, 1500, 1.604e-04, 0.395791, 10.0917, 0, 0.10940, 0.33076))
modelPointsC0p01.append(ModelPoint(0.01, 1750, 5.408e-05, 0.444453, 11.5976, 0, 0.04682, 0.21638))
modelPointsC0p01.append(ModelPoint(0.01, 2000, 1.853e-05, 0.47302,  13.6845, 0, 0.02463, 0.15695))
modelPointsC0p01.append(ModelPoint(0.01, 3000, 4.703e-07, 0.53921,  22.4852, 0, 0.00193, 0.04389))
#c = 0.05
modelPointsC0p05 = []
modelPointsC0p05.append(ModelPoint(0.05, 1750, 1.331e-03, 0.437371, 13.7476, 0, 0.04682, 0.21638))
modelPointsC0p05.append(ModelPoint(0.05, 2000, 4.665e-04, 0.471764, 16.795,  0, 0.02463, 0.15695))
modelPointsC0p05.append(ModelPoint(0.05, 2500, 7.226e-05, 0.517763, 20.4539, 0, 0.00445, 0.06668))
modelPointsC0p05.append(ModelPoint(0.05, 2750, 2.803e-05, 0.528236, 22.7374, 0, 0.00182, 0.04265))
modelPointsC0p05.append(ModelPoint(0.05, 3000, 1.169e-05, 0.529311, 25.2982, 0, 0.00193, 0.04389))
#c = 0.1
modelPointsC0p1 = []
modelPointsC0p1.append(ModelPoint(0.1, 2250, 7.04e-04, 0.499742, 26.3803, 0, 0.01589, 0.12605))
modelPointsC0p1.append(ModelPoint(0.1, 2500, 2.79e-04, 0.517166, 30.8038, 0, 0.00663, 0.08143))
modelPointsC0p1.append(ModelPoint(0.1, 2750, 1.14e-04, 0.527167, 35.9283, 0, 0.00321, 0.05667))
modelPointsC0p1.append(ModelPoint(0.1, 3000, 4.68e-05, 0.527422, 38.7399, 0, 0.00193, 0.04390))
modelPointsC0p1.append(ModelPoint(0.1, 3250, 1.19e-05, 0.529048, 41.8178, 0, 0.00326, 0.05711))
modelPointsC0p1.append(ModelPoint(0.1, 3500, 7.7e-05,  0.524787, 40.2991, 0, 0.00326, 0.05711))

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
elif sys.argv[1]=='tables':
  print 'tables: DoTablesAllPoints'
  DoTablesAllPoints(lumi,3)
# TODO: add yield calculation into all
# TODO: calculate signal efficiency from here
elif sys.argv[1]=='all':
  print 'all: limits+plots'
  print 'limits: DoLimitsAllPoints'
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileNameBase)
  print 'plots: DoPlotsAllPoints'
  DoPlotsAllPoints(lumi)
  print 'tables: DoTablesAllPoints'
  DoTablesAllPoints(lumi,3)
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






