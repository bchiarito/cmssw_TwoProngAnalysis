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
import datetime
import itertools

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
# mapping of couplings/masses to crossSections --> taken from AN-12-305
totalXSecs0p01 = dict()
totalXSecs0p01[750] =  1.023e-02
totalXSecs0p01[1000] = 2.072e-03
totalXSecs0p01[1250] = 5.21e-04
totalXSecs0p01[1500] = 1.604e-04
totalXSecs0p01[1750] = 5.408e-05
totalXSecs0p01[2000] = 1.853e-05
totalXSecs0p01[2250] = 7.207e-06
totalXSecs0p01[2500] = 2.945e-06
totalXSecs0p01[3000] = 4.703e-07
# 0.05
totalXSecs0p05 = dict()
totalXSecs0p05[1250] = 1.28e-02
totalXSecs0p05[1500] = 3.866e-03
totalXSecs0p05[1750] = 1.331e-03
totalXSecs0p05[2000] = 4.665e-04
totalXSecs0p05[2250] = 1.774e-04
totalXSecs0p05[2500] = 7.226e-05
totalXSecs0p05[2750] = 2.803e-05
totalXSecs0p05[3000] = 1.169e-05
# 0.1
totalXSecs0p1 = dict()
totalXSecs0p1[1500] = 1.5e-02
totalXSecs0p1[1750] = 5.2e-03
totalXSecs0p1[2000] = 1.8e-03
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
    return -1
    #print 'GetHalfWidth: Coupling',coupling,'mass',mass,'not recognized; quitting'
    #exit()


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


def FillKFactors(modelPointArray):
  myCoupling = modelPointArray[0].coupling
  massToKFactorDict = dict()
  with open(kFactorFile, 'r') as file:
    lines = file.readlines()
  i = 0
  while i < len(lines):
    line = lines[i]
    line.strip('\n')
    if "Coupling" in line:
      coupling = float(line.split('= ')[1])
      if coupling == myCoupling:
        break
    i+=1
  if i>=len(lines):
    print 'Could not find k-factors in file:',kFactorFile,'for coupling=',myCoupling
    return
  i+=2 # get past ----- line
  line = lines[i]
  while i < len(lines):
    if '--' in line:
      break
    lineSplit = line.split()
    massToKFactorDict[int(float(lineSplit[0]))] = float(lineSplit[1])
    i+=1
    line = lines[i]
  for mp in modelPointArray:
    mp.kFactor = massToKFactorDict[mp.mass]


# TODO: remove hardcoding; instead use list of couplings to run over
def DoLimitsAllPoints(cl95MacroPathAndName,lumi,lumiErr,limitsFileName):
  print 'Run for coupling 0.01'
  with open(optimizationFileName+'0p01.txt', 'r') as file:
    lines = file.readlines()
    readModelPointsC0p01 = ReadFromLines(lines)
    if UseKFactor:
      FillKFactors(readModelPointsC0p01)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p01,limitsFileName+'0p01.txt',UseKFactor)
  print 'Run for coupling 0.05'
  with open(optimizationFileName+'0p05.txt', 'r') as file:
    lines = file.readlines()
    readModelPointsC0p05 = ReadFromLines(lines)
    if UseKFactor:
      FillKFactors(readModelPointsC0p05)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p05,limitsFileName+'0p05.txt',UseKFactor)
  print 'Run for coupling 0.1'
  with open(optimizationFileName+'0p1.txt', 'r') as file:
    lines = file.readlines()
    readModelPointsC0p1 = ReadFromLines(lines)
    if UseKFactor:
      FillKFactors(readModelPointsC0p1)
  ComputeLimits(cl95MacroPathAndName,lumi,lumiErr,readModelPointsC0p1,limitsFileName+'0p1.txt',UseKFactor)
  subprocess.call(['mv','cls_plots',limitsOutputDir])


## TODO
# function to read in the coupling we ran over
#def ReadResultsAllPoints():
#  resultsByCoupling = {}
#  for c in couplingsToRunOver:
#    # make coupling string
#    couplingString = f(c)
#    with open(limitsFileName+couplingString+'.txt', 'r') as file:
#      readModelPoints = ReadFromFile(file)
#    resultsByCoupling[c] = readModelPoints
#  return resultsByCoupling


def DoPlotsAllPoints(lumi,rootFile,pathToTDRStyle):
  # set up tdrStyle
  gROOT.ProcessLine('.L '+pathToTDRStyle)
  gROOT.ProcessLine('setTDRStyle()')
  print 'Run for coupling 0.01'
  with open(limitsFileName+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  for plotObsLim in [True,False]:
    PlotBands(readModelPoints0p01,plotObsLim,lumi,rootFile,plotsOutputDir)
  m0p01,xs,mExp0p01,xsE = GetMassLimit(readModelPoints0p01)
  print string.ljust('Coupling: '+str(readModelPoints0p01[0].coupling),14),
  print ' Observed limit mass: %0.2f'%m0p01
  print '                Expected limit mass: %0.2f'%mExp0p01
  print '                Observed XSec limit: %0.6f'%xs
  print '                Expected XSec limit: %0.6f'%xsE
  # 0.05
  print 'Run for coupling 0.05'
  with open(limitsFileName+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  for plotObsLim in [True,False]:
    PlotBands(readModelPoints0p05,plotObsLim,lumi,rootFile,plotsOutputDir)
  m0p05,xs,mExp0p05,xsE = GetMassLimit(readModelPoints0p05)
  print string.ljust('Coupling: '+str(readModelPoints0p05[0].coupling),14),
  print ' Observed limit mass: %0.2f'%m0p05
  print '                Expected limit mass: %0.2f'%mExp0p05
  print '                Observed XSec limit: %0.6f'%xs
  print '                Expected XSec limit: %0.6f'%xsE
  # 0.1
  print 'Run for coupling 0.1'
  with open(limitsFileName+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  for plotObsLim in [True,False]:
    PlotBands(readModelPoints0p1,plotObsLim,lumi,rootFile,plotsOutputDir)
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
  CouplingVsMassPlot(couplings,mLimExp,mLimObs,rootFile,lumi,plotsOutputDir)
  # efficiencies for all couplings/masses on same axes
  print 'PlotAllEfficiencies'
  PlotAllEfficiencies([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesMScale([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1], lumi, rootFile, plotsOutputDir)
  PlotAllEfficienciesMRes([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesPileup([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  # half widths
  print 'PlotAllHalfWidths'
  PlotAllHalfWidths([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  # exp BG
  print 'PlotAllExpBGs'
  PlotAllExpBGs([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  # limit plot for all couplings on same axes
  print 'PlotAllBands'
  PlotAllBands([readModelPoints0p01,readModelPoints0p05,readModelPoints0p1],lumi,rootFile,plotsOutputDir)
  # make table
  print 'DoTablesAllPoints'
  DoTablesAllPoints(lumi)


def DoTablesAllPoints(lumi):
  with open(limitsFileName+'0p01.txt', 'r') as file:
    readModelPoints0p01 = ReadFromFile(file)
  # 0.05
  with open(limitsFileName+'0p05.txt', 'r') as file:
    readModelPoints0p05 = ReadFromFile(file)
  # 0.1
  with open(limitsFileName+'0p1.txt', 'r') as file:
    readModelPoints0p1 = ReadFromFile(file)
  # now make table output
  tableTitleString=string.ljust('Kmpl',6)+string.ljust('Mass',6)+string.ljust('Window',11)
  tableTitleString+=string.ljust('Eff.',8)
  tableTitleString+=string.ljust('ExpSig.',9)
  tableTitleString+=string.center('ExpBkg.',41)+string.center('Obs.',10)
  tableTitleString+=string.center('ThXSec',10)+string.center('ExpLim.',12)
  tableTitleString+=string.center('ObsLim.',12)
  print
  print tableTitleString
  for modelPoint in readModelPoints0p01+readModelPoints0p05+readModelPoints0p1:
    print modelPoint.StringTableLine(lumi)
  # twiki table
  twikiTitleString='|* '+string.ljust('Coupling',9)+' *|* '+string.ljust('Mass [GeV]',10)+' *|* '
  twikiTitleString+='Mass Window *|* '
  twikiTitleString+=string.ljust('Efficiency',8)+' *|* '
  twikiTitleString+=string.center('Exp. Sig. Evts.',10)+' *|* '
  twikiTitleString+=string.center('Exp. Bkg. Evts.',17)+' *|* '+string.center('Obs. Data',10)+' *|* '
  twikiTitleString+=string.center('Exp. Lim. [pb]',13)+' *|* '
  twikiTitleString+=string.center('Obs. Lim. [pb]',13)+' *|'
  print
  print twikiTitleString
  for modelPoint in readModelPoints0p01+readModelPoints0p05+readModelPoints0p1:
    print modelPoint.TwikiTableLine(lumi)
  # twiki mass limits
  twikiMassTitleString = '|* Coupling *|* Expected XSec Limit [pb] *|* Observed XSec Limit [pb] *|* Expected Mass Limit [GeV] *|* Observed Mass Limit [GeV] *|'
  twikiLine = ''
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p01)
  twikiLine+='| 0.01 | '+str(round(xsLimExp,5))+' | '+str(round(xsLimObs,5))+' | '+str(int(massLimExp))+' | '+str(int(massLimObs))+' |\n'
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p05)
  twikiLine+='| 0.05 | '+str(round(xsLimExp,5))+' | '+str(round(xsLimObs,5))+' | '+str(int(massLimExp))+' | '+str(int(massLimObs))+' |\n'
  massLimObs, xsLimObs, massLimExp, xsLimExp = GetMassLimit(readModelPoints0p1)
  twikiLine+='| 0.1  | '+str(round(xsLimExp,5))+' | '+str(round(xsLimObs,5))+' | '+str(int(massLimExp))+' | '+str(int(massLimObs))+' |\n'
  print
  print twikiMassTitleString
  print twikiLine
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
  # TODO remove hardcoding of couplings to run over
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
  

def DoOptimizeAllPoints(optimizationRootFile):
  # configurable options for optimization
  maxWindowRange = 600 # bins/GeV
  useAsymmWindow = False
  useSSB = True
  print 'Run for coupling 0.01'
  colorIndex = 2 #TODO add this into modelpoint itself?
  with open(optimizationFileName+'0p01.txt', 'w') as file:
    OptimizeSignalMassWindows(
        rootFileLocationDataFake,rootFileBackgroundMC,modelPointsC0p01,lumi,useAsymmWindow,useSSB,maxWindowRange,file,optimizationRootFile,colorIndex,optimizationOutputDir,extraWindowMargin,DataSample)
  print 'Run for coupling 0.05'
  colorIndex = 4
  with open(optimizationFileName+'0p05.txt', 'w') as file:
    OptimizeSignalMassWindows(
        rootFileLocationDataFake,rootFileBackgroundMC,modelPointsC0p05,lumi,useAsymmWindow,useSSB,maxWindowRange,file,optimizationRootFile,colorIndex,optimizationOutputDir,extraWindowMargin,DataSample)
  print 'Run for coupling 0.1'
  colorIndex = 8
  with open(optimizationFileName+'0p1.txt', 'w') as file:
    OptimizeSignalMassWindows(
        rootFileLocationDataFake,rootFileBackgroundMC,modelPointsC0p1,lumi,useAsymmWindow,useSSB,maxWindowRange,file,optimizationRootFile,colorIndex,optimizationOutputDir,extraWindowMargin,DataSample)
  # make multigraphs for all masses/couplings
  MakeOptHalfWindowVsMassMultigraph(optimizationRootFile)
  MakeOptMassWindowsVsMassMultiGraph(optimizationRootFile)
  MakeOptSSBValueVsMassMultigraph(optimizationRootFile)
  # print wiki-style table with links to plots (made/copied later by plotting script)
  print 'Twiki-style table of mass windows'
  print
  print '| *Coupling* | *Mass (GeV)* | *Mass Window Low* | *Mass Window High* | *Mass Window Half-Width* | *Optimization Plot* | *PDF* |'
  for mp in modelPointsC0p01:
    print '|',mp.coupling,'|',mp.mass,'|',mp.optMassWindowLow,'|',mp.optMassWindowHigh,'|',(mp.optMassWindowHigh-mp.optMassWindowLow)/2+0.5,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p01_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p01_m'+str(mp.mass)+'.png" alt="optimization_K0p01_m'+str(mp.mass)+'" width="400" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p01_m'+str(mp.mass)+'.pdf>PDF Version</a>|'
  print '| |||||'
  for mp in modelPointsC0p05:
    print '|',mp.coupling,'|',mp.mass,'|',mp.optMassWindowLow,'|',mp.optMassWindowHigh,'|',(mp.optMassWindowHigh-mp.optMassWindowLow)/2+0.5,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p05_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p05_m'+str(mp.mass)+'.png" alt="optimization_K0p05_m'+str(mp.mass)+'" width="400" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p05_m'+str(mp.mass)+'.pdf>PDF Version</a>|'
  print '| |||||'
  for mp in modelPointsC0p1:
    print '|',mp.coupling,'|',mp.mass,'|',mp.optMassWindowLow,'|',mp.optMassWindowHigh,'|',(mp.optMassWindowHigh-mp.optMassWindowLow)/2+0.5,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p1_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p1_m'+str(mp.mass)+'.png" alt="optimization_K0p1_m'+str(mp.mass)+'" width="400" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/ssbOpt_k0p1_m'+str(mp.mass)+'.pdf>PDF Version</a>|'
  # print out signal/background mass window plot code for twiki
  print 'Twiki-style table of signal/background mass plot'
  print '| *Coupling* | *Mass (!GeV)* | *Plot* | *PDF* |'
  for mp in modelPointsC0p01:
    print '|0.01|',mp.mass,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p01_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p01_m'+str(mp.mass)+'.png" alt="optimization_K0p01_m'+str(mp.mass)+'" width="600" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p01_m'+str(mp.mass)+'.pdf>PDF Version</a>|'
  print '| ||||'
  for mp in modelPointsC0p05:
    print '|0.05|',mp.mass,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p05_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p05_m'+str(mp.mass)+'.png" alt="optimization_K0p05_m'+str(mp.mass)+'" width="600" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p05_m'+str(mp.mass)+'.pdf>PDF Version</a>|'
  print '| ||||'
  for mp in modelPointsC0p1:
    print '|0.1|',mp.mass,'|',
    print '<a href="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p1_m'+str(mp.mass)+'.png"><img src="http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p1_m'+str(mp.mass)+'.png" alt="optimization_K0p1_m'+str(mp.mass)+'" width="600" /></a>|<a href=http://scooper.web.cern.ch/scooper/exoDiPhotons/twikiPlots/signalBackgroundOptWindows_k0p1_m'+str(mp.mass)+'.pdf>PDF Version</a>|'


def DoOptimizationPlots(optimizationRootFile):
  # toggle plot s/sqrt(b) off
  MakeOptGraphImages(optimizationRootFile,optimizationOutputDir,modelPointsC0p01+modelPointsC0p05+modelPointsC0p1,True)
  # above takes a lot of time--many points
  # do images
  MakeOptMassWindowVsMassImages(optimizationRootFile,optimizationOutputDir)
  MakeSmoothedMassWindowVsMassImages(optimizationRootFile,optimizationOutputDir)
  MakeOptSSBVsMassImages(optimizationRootFile, optimizationOutputDir)
  MakeOptMassWindowSignalBackgroundImages(optimizationRootFile,optimizationOutputDir,modelPointsC0p01+modelPointsC0p05+modelPointsC0p1)


def DoCalculateYieldsAllPoints():
  print 'Yields: Run for coupling 0.01'
  with open(optimizationFileName+'0p01.txt', 'w') as file:
    CalculateYieldsForMassRanges(rootFileLocationDataFake, rootFileBackgroundMC, modelPointsC0p01, lumi, 3, file,DataSample)
  print 'Yields: Run for coupling 0.05'
  with open(optimizationFileName+'0p05.txt', 'w') as file:
    CalculateYieldsForMassRanges(rootFileLocationDataFake, rootFileBackgroundMC, modelPointsC0p05, lumi, 3, file,DataSample)
  print 'Yields: Run for coupling 0.1'
  with open(optimizationFileName+'0p1.txt', 'w') as file:
    CalculateYieldsForMassRanges(rootFileLocationDataFake, rootFileBackgroundMC, modelPointsC0p1, lumi, 3, file,DataSample)


def GetConfigurationString():
  configString='-----------------------------------------------\n'+'Running with configuration:'+'\n'
  configString+='lumi='+str(lumi)+'\n'
  configString+='lumiErr='+str(lumiErr)+'\n'+'-----------------------------------------------\n'
  configString+='signalRootFileLocation='+signalRootFileLocation+'\n'
  configString+='signalPoints K0p01='+str(masses0p01)+'\n'
  configString+='signalPoints K0p05='+str(masses0p05)+'\n'
  configString+='signalPoints K0p1='+str(masses0p1)+'\n'
  configString+='Data/Fake rootFileLocation='+rootFileLocationDataFake+'\n'
  configString+='DataSample='+DataSample+'\n'
  configString+='BackgroundMC root file='+rootFileBackgroundMC+'\n'+'-----------------------------------------------\n'
  configString+='extraWindowMargin='+str(extraWindowMargin)+'\n'
  configString+='UseKFactor='+str(UseKFactor)+', k-factor file='+kFactorFile+'\n'
  # overall systematics
  configString+='SigPUSyst='+str(SigPUSyst)+'\n'
  configString+='SigPDFSyst='+str(SigPDFSyst)+'\n'
  configString+='SigScaleFactorSystOneGamma='+str(SigScaleFactorSystOneGamma)+'\n'
  configString+='SigPtSFSystOneGamma='+str(SigPtSFSystOneGamma)+'\n'
  configString+='BGSystsFromHists'+'\n'+'-----------------------------------------------\n'
  configString+='optimizationOutputDir='+optimizationOutputDir+'\n'
  configString+='limitsOutputDir='+limitsOutputDir+'\n'
  configString+='limitsFileNameBase='+limitsFileNameBase+'\n'
  configString+='optimizationFileNameBase='+optimizationFileNameBase+'\n'+'-----------------------------------------------\n'
  return configString


def Usage():
  print 'Usage: python RunLimitsAndPlots.py [arg] where arg can be:'
  print '    all           --> run yields, limits, and plots (see below)'
  print '    limits        --> Compute limits for all model points & write out results'
  print '    plots         --> Read limit results from text files and make limit plots'
  print '    tables        --> Read limit results from text files and make results tables (text/latex)'
  print '    yields        --> Calculate event yields from histograms in root files (from CreateHistogramFiles)'
  print '    optimize      --> Calculate optimal inv. mass windows using histograms in root files (from CreateHistogramFiles)'
  print '    optimizePlots --> Print optimization graphs as images (use after optimize has been run)'




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
now = datetime.datetime.now()
Date = now.strftime("%b%d")
#outputDirBase = Date.lower()+'_results_lumiErr2p6_newFakeAndBGSysts_sherpaBG_noDPhi_symmWindowSSBOpt_BGSystsFromHists_0pctOptMWindowMarginSmoothed_loosePFID'
outputDirBase = 'mar07_results_lumiErr2p6_newFakeAndBGSysts_sherpaBG_noDPhi_symmWindowSSBOpt_BGSystsFromHists_0pctOptMWindowMarginSmoothed_mediumPFID'
#outputDirBase = 'mar08_results_lumiErr2p6_newFakeAndBGSysts_sherpaBG_noDPhi_symmWindowSSBOpt_BGSystsFromHists_0pctOptMWindowMarginSmoothed_loosePFID'
#outputDirBase = 'oct17_TightPFID_optLimitsPlots/oct17_results_lumiErr2p6_newFakeAndBGSysts_sherpaBG_noDPhi_symmWindowSSBOpt_BGSystsFromHists_0pctOptMWindowMarginSmoothed'
optimizationOutputDir = outputDirBase+'_optimization'
limitsOutputDir = outputDirBase+'_limits'
plotsOutputDir = outputDirBase+'_plots3'
# limit results file base name
limitsFileNameBase = 'limits_k_'
limitsFileName = limitsOutputDir+'/'+limitsFileNameBase
limitsPlotFileName = limitsOutputDir+'/plots.root'
# optimization results file base name
optimizationFileNameBase = 'optimization_k_'
optimizationFileName = optimizationOutputDir+'/'+optimizationFileNameBase
optimizationPlotFileName = optimizationOutputDir+'/plots.root'
# location of signal root files from CreateHistogramFiles code
#signalRootFileLocation = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_deltaPhi2p8_19p6invFb/'
#signalRootFileLocation = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb/'
#signalRootFileLocation = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_mediumID/'
signalRootFileLocation = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_looseID/'
# location of data/fake root files from CreateHistogramFiles code
## OLD #rootFileLocationDataFake = '/afs/cern.ch/work/c/charaf/public/DiPhotonTrees/Histograms/ExoDiPhotonAnalyzer_PFDec14th_DataABCD/'
#rootFileLocationDataFake = '/afs/cern.ch/work/c/charaf/public/ForDiPhoton/'
#DataSample = 'ExoDiPhotonAnalyzer_PFDec14th_DataABCD'
#rootFileLocationDataFake = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_mediumID/'
#DataSample = 'ExoDiPhotonAnalyzer_MediumDataABCD2012'
rootFileLocationDataFake = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_looseID/'
DataSample = 'ExoDiPhotonAnalyzer_LooseDataABCD2012'
# location of backgroundMC root files from CreateHistogramFiles code
#rootFileBackgroundMC = '/afs/cern.ch/work/c/charaf/public/ForDiPhoton/histograms_diphoton_tree_MC_all.root' # SHERPA, no dphi
## OLD #rootFileBackgroundMC= '/afs/cern.ch/work/c/charaf/public/DiPhotonTrees/Histograms/diphoton_tree_MC_all/histograms_diphoton_tree_MC_all.root' # PYTHIA
#rootFileBackgroundMC = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_mediumID/histograms_diphoton_tree_MC_all.root'
rootFileBackgroundMC = '/afs/cern.ch/user/s/scooper/work/public/DiPhotonHistograms/PFID_19p6invFb_looseID/histograms_diphoton_tree_MC_all.root'
kFactorFile = 'RS-KF-LHC-8TeV-y1.4442-ptcut80.dat'
# use k-factors to compute limit (optimization always done with k-factor=1)
UseKFactor = False
# use extra window margin to minimize mScale/mRes systematic effects on signal
extraWindowMargin = 0.00 # fraction to expand window by, e.g., 0.01 --> expand edges by 1%
# Declarations of Lumi and model points to consider -- must have xsec, etc. defined above
lumi = 19620.
lumiErr = lumi*0.026
masses0p01 = [750,1000,1250,1500,1750,2000,2250,2500,3000]
masses0p05 = [1250,1500,1750,2000,2250,2500,2750,3000]
masses0p1 = [1500,1750,2000,2250,2500,2750,3000,3250,3500]
#

configString = GetConfigurationString()
print configString

# List signal histogram file template; rest of quantities are filled from functions (xsec, width, etc.) or histograms in the files
# initialize signal points
#signalHistogramFilesPathTemplate = signalRootFileLocation+'diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia_merged/histograms_diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia_merged.root'
#FIXME --> added pythia6 for med/loose signal samples
signalHistogramFilesPathTemplate = signalRootFileLocation+'diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia6_merged/histograms_diphoton_tree_RSGravToGG_kMpl-{coupling}_M-{mass}_TuneZ2star_8TeV-pythia6_merged.root'
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
elif sys.argv[1]=='optimizePlots':
  print 'optimizePlots: print plot images from optimization'
  rootFile = TFile(optimizationPlotFileName,'update')
  DoOptimizationPlots(rootFile)
elif sys.argv[1]=='optimize':
  print 'optimize: OptimizeSignalMassWindows'
  print 'warning: will overwrite file',optimizationPlotFileName,'if it already exists.'
  if not os.path.isdir(optimizationOutputDir):
    os.mkdir(optimizationOutputDir)
  with open(optimizationOutputDir+'/config.txt', 'w') as file:
    file.write(configString)
  # recreate tfile for optimize step
  rootFile = TFile(optimizationPlotFileName,'recreate')
  DoOptimizeAllPoints(rootFile)
  DoOptimizationPlots(rootFile)
  if not os.path.isdir(plotsOutputDir):
    os.mkdir(plotsOutputDir)
  PlotAllEfficiencies([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesMScale([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesMRes([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesPileup([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
elif sys.argv[1]=='yields':
  print 'yields: CalculateYieldsForMassRanges'
  if not os.path.isdir(optimizationOutputDir):
    os.mkdir(optimizationOutputDir)
  with open(optimizationOutputDir+'/config.txt', 'w') as file:
    file.write(configString)
  rootFile = TFile(optimizationPlotFileName,'recreate')
  DoCalculateYieldsAllPoints()
  if not os.path.isdir(plotsOutputDir):
    os.mkdir(plotsOutputDir)
  PlotAllEfficiencies([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
elif sys.argv[1]=='limits':
  print 'limits: DoLimitsAllPoints'
  print 'warning: will overwrite file',limitsPlotFileName,'if it already exists.'
  if not os.path.isdir(limitsOutputDir):
    os.mkdir(limitsOutputDir)
  with open(limitsOutputDir+'/config.txt', 'w') as file:
    file.write(configString)
  rootFile = TFile(limitsPlotFileName,'recreate')
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileName)
elif sys.argv[1]=='plots':
  print 'plots: DoPlotsAllPoints'
  if not os.path.isdir(plotsOutputDir):
    os.mkdir(plotsOutputDir)
  rootFile = TFile(limitsPlotFileName,'update')
  DoPlotsAllPoints(lumi,rootFile,pathToTDRStyle)
elif sys.argv[1]=='tables':
  print 'tables: DoTablesAllPoints'
  rootFile = TFile(limitsPlotFileName,'update')
  DoTablesAllPoints(lumi)
elif sys.argv[1]=='all':
  if not os.path.isdir(optimizationOutputDir):
    os.mkdir(optimizationOutputDir)
  print 'all: optimize+limits+plots'
  print 'optimize: OptimizeSignalMassWindows'
  print 'warning: will overwrite file',optimizationPlotFileName,'if it already exists.'
  # recreate tfile for optimize step
  rootFile = TFile(optimizationPlotFileName,'recreate')
  with open(optimizationOutputDir+'/config.txt', 'w') as file:
    file.write(configString)
  DoOptimizeAllPoints(rootFile)
  DoOptimizationPlots(rootFile)
  #
  #print 'all: stdYields+limits+plots'
  #rootFile = TFile(optimizationPlotFileName,'recreate')
  #DoCalculateYieldsAllPoints() # std/old mass windows
  if not os.path.isdir(plotsOutputDir):
    os.mkdir(plotsOutputDir)
  PlotAllEfficiencies([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesMScale([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesMRes([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  PlotAllEfficienciesPileup([modelPointsC0p01,modelPointsC0p05,modelPointsC0p1],lumi,rootFile,plotsOutputDir)
  rootFile.Close()
  print 'limits: DoLimitsAllPoints'
  print 'warning: will overwrite file',limitsPlotFileName,'if it already exists.'
  if not os.path.isdir(limitsOutputDir):
    os.mkdir(limitsOutputDir)
  with open(limitsOutputDir+'/config.txt', 'w') as file:
    file.write(configString)
  rootFile = TFile(limitsPlotFileName,'recreate')
  DoLimitsAllPoints(cl95MacroPath+cl95MacroName,lumi,lumiErr,limitsFileName)
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






