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


