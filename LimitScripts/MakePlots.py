#!/usr/bin/env python

#
# Support functions to make limits plots from input txt files
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


def ConvertToPng(baseName):
  subprocess.call(['gs','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+baseName+'.png',baseName+'.eps'])
  print 'file',baseName+'.png','created.'


def makeExpBGArrays(modelPointArray):
  masses = []
  expBGs = []
  for modelPoint in modelPointArray:
    masses.append(modelPoint.mass)
    expBGs.append(modelPoint.nBackground)
  return masses, expBGs


def makeWidthArrays(modelPointArray):
  masses = []
  halfWidths = []
  for modelPoint in modelPointArray:
    masses.append(modelPoint.mass)
    halfWidths.append(modelPoint.halfWidth)
  return masses, halfWidths


def makeEffArrays(modelPointArray):
  masses = []
  efficiency = []
  for modelPoint in modelPointArray:
    masses.append(modelPoint.mass)
    efficiency.append(modelPoint.totalEff)
  return masses, efficiency


def makeArrays(modelPointArray):
  errMass = [0]*len(modelPointArray)
  errUp = []
  errDn = []
  err2Up = []
  err2Dn = []
  masses = []
  expUp = []
  expDn = []
  exp2Up = []
  exp2Dn = []
  limitExp = []
  limitObs = []
  totalXSec = []
  # test last model point to see if there are limits for all points
  # if not, return empty arrays
  try:
    value = float(modelPointArray[len(modelPointArray)-1].expLimit)
  except ValueError:
    return errMass, errUp, errDn, err2Up, err2Dn, masses, expUp, expDn, exp2Up, exp2Dn, limitExp, limitObs, totalXSec
  for modelPoint in modelPointArray:
    errUp.append(modelPoint.expLimitOneSigmaHigh-modelPoint.expLimit)
    errDn.append(modelPoint.expLimit-modelPoint.expLimitOneSigmaLow)
    err2Up.append(modelPoint.expLimitTwoSigmaHigh-modelPoint.expLimit)
    err2Dn.append(modelPoint.expLimit-modelPoint.expLimitTwoSigmaLow)
    masses.append(modelPoint.mass)
    expUp.append(modelPoint.expLimitOneSigmaHigh)
    expDn.append(modelPoint.expLimitOneSigmaLow)
    exp2Up.append(modelPoint.expLimitTwoSigmaHigh)
    exp2Dn.append(modelPoint.expLimitTwoSigmaLow)
    limitExp.append(modelPoint.expLimit)
    limitObs.append(modelPoint.obsLimit)
    totalXSec.append(modelPoint.totalXSec)
  return errMass, errUp, errDn, err2Up, err2Dn, masses, expUp, expDn, exp2Up, exp2Dn, limitExp, limitObs, totalXSec


def SetCustomGStyle():
  gStyle.SetOptStat(0)
  gStyle.SetTextSize(18)
  gStyle.SetTitleBorderSize(0)
  gStyle.SetTitleFillColor(0)
  gStyle.SetTitleFontSize(0.04)
  gStyle.SetTitleStyle(00)
  gStyle.SetStatBorderSize(1)
  gStyle.SetStatColor(0)
  gStyle.SetOptStat(110)
  gStyle.SetCanvasColor(0)
  gStyle.SetCanvasBorderMode(0)
  gStyle.SetPadColor(0)
  gStyle.SetPadBorderMode(0)
  gStyle.SetFrameLineWidth(2)
  #  gStyle.SetFillColor(0)
  gStyle.SetHistLineWidth(2)
  gStyle.SetLineWidth(2)
  gStyle.SetFuncWidth(2)
  gStyle.SetFuncColor(2)
  gStyle.SetLabelSize(0.073)
  gStyle.SetOptStat(0)
  gStyle.SetOptFit(111111)


def MakeOptimizationGraph(peakMass,modelPoint,minMassTried,maxMassTried,massRangesTried,ssbTried,optSSBIndex,useAsymmWindow,rootFile):
  rootFile.cd()
  # make optimization graph (1-D for symm window; 2-D for asymm window)
  if useAsymmWindow:
    graph = TH2F('test','test',int(peakMass-minMassTried)+1,minMassTried,peakMass+1,int(maxMassTried-peakMass),peakMass,maxMassTried)
    for massRange,ssb in itertools.izip(massRangesTried,ssbTried):
      graph.Fill(massRange[0],massRange[1]-1,ssb)
    graph.SetTitle('Optimization for RS Graviton mass='+str(modelPoint.mass)+'GeV, coupling='+str(modelPoint.coupling))
    graph.SetName('ssbOpt_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graph.GetXaxis().SetTitle('Mass window low [GeV]')
    graph.GetYaxis().SetTitle('Mass window high [GeV]')
    #graph.GetZaxis().SetTitle('S/#sqrt{S+B}')
  else:
    halfWindowSizesTried = [(massRangeTried[1]-massRangeTried[0])/2.0 for massRangeTried in massRangesTried]
    graph = TGraph(len(halfWindowSizesTried), array.array("f",halfWindowSizesTried),array.array("f",ssbTried))
    graph.SetTitle('Optimization for RS Graviton mass='+str(modelPoint.mass)+'GeV, coupling='+str(modelPoint.coupling))
    graph.SetName('ssbOpt_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graph.GetXaxis().SetTitle('Mass window half size [GeV]')
    graph.GetYaxis().SetTitle('S/#sqrt{S+B}')
    graphOpt = TGraph(1, array.array("f",[halfWindowSizesTried[optSSBIndex]]),array.array("f",[ssbTried[optSSBIndex]]))
    graphOpt.SetName('ssbOptPoint_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graphOpt.Write()
  graph.Write()


def MakeOptimizationGraphSRootB(peakMass,modelPoint,minMassTried,maxMassTried,massRangesTried,sRootBTried,optSRootBIndex,useAsymmWindow,rootFile):
  rootFile.cd()
  # make optimization graph (1-D for symm window; 2-D for asymm window)
  if useAsymmWindow:
    graph = TH2F('test','test',int(peakMass-minMassTried)+1,minMassTried,peakMass+1,int(maxMassTried-peakMass),peakMass,maxMassTried)
    for massRange,ssb in itertools.izip(massRangesTried,sRootBTried):
      graph.Fill(massRange[0],massRange[1]-1,ssb)
    graph.SetTitle('Optimization for RS Graviton mass='+str(modelPoint.mass)+'GeV, coupling='+str(modelPoint.coupling))
    graph.SetName('srootbOpt_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graph.GetXaxis().SetTitle('Mass window low [GeV]')
    graph.GetYaxis().SetTitle('Mass window high [GeV]')
    #graph.GetZaxis().SetTitle('S/#sqrt{S+B}')
  else:
    halfWindowSizesTried = [(massRangeTried[1]-massRangeTried[0])/2.0 for massRangeTried in massRangesTried]
    graph = TGraph(len(halfWindowSizesTried), array.array("f",halfWindowSizesTried),array.array("f",sRootBTried))
    graph.SetTitle('Optimization for RS Graviton mass='+str(modelPoint.mass)+'GeV, coupling='+str(modelPoint.coupling))
    graph.SetName('srootbOpt_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graph.GetXaxis().SetTitle('Mass window half size [GeV]')
    graph.GetYaxis().SetTitle('S/#sqrt{B}')
    graphOpt = TGraph(1, array.array("f",[halfWindowSizesTried[optSRootBIndex]]),array.array("f",[sRootBTried[optSRootBIndex]]))
    graphOpt.SetName('srootbOptPoint_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
    graphOpt.Write()
  graph.Write()


def MakeSmoothedWindowVsMassPlot(coupling,masses,optMinMasses,optMaxMasses,colorIndex,rootFile):
  # make plot of opt. halfWindowSize vs. mass/coupling
  if len(masses) < 1:
    return
  optHalfWindowSizes = [(optMassHigh-optMassLow-1) for optMassHigh,optMassLow in itertools.izip(optMaxMasses,optMinMasses)]
  #                   # -1, since we get the top edge of the maxBin as the upper mass limit
  graph = TGraph(len(masses), array.array("f",masses),array.array("f",optHalfWindowSizes))
  plotname = "optSmoothedWindowVsMassK"+str(coupling)
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  graph.SetName(savename.Data())
  graph.SetMarkerColor(colorIndex)
  graph.SetLineColor(colorIndex)
  graph.Write()


def MakeOptHalfWindowVsMassMultigraph(rootFile):
  rootFile.cd()
  graph0p01 = MakeNullPointer(TGraph)
  graph0p05 = MakeNullPointer(TGraph)
  graph0p1 = MakeNullPointer(TGraph)
  try:
    rootFile.GetObject('optSmoothedWindowVsMassK0p01',graph0p01)
    rootFile.GetObject('optSmoothedWindowVsMassK0p05',graph0p05)
    rootFile.GetObject('optSmoothedWindowVsMassK0p1',graph0p1)
  except LookupError:
    pass
  c = TCanvas()
  c.SetName('optSmoothedWindowsVsMassAllCanvas')
  c.SetTitle('')
  c.cd()
  mg = TMultiGraph()
  if graph0p01:
    mg.Add(graph0p01)
  if graph0p05:
    mg.Add(graph0p05)
  if graph0p1:
    mg.Add(graph0p1)
  mg.Draw('ap')
  mg.GetXaxis().SetTitle("Mass [GeV]")
  mg.GetYaxis().SetTitle("Opt. HalfWindow size [GeV]")
  mg.SetName('optSmoothedWindowsVsMassAll')
  legend = TLegend(0.42,0.71,0.73,0.92)
  if graph0p01:
    legend.AddEntry(graph0p01," #tilde{k} = "+str(0.01),"l")
  if graph0p05:
    legend.AddEntry(graph0p05," #tilde{k} = "+str(0.05),"l")
  if graph0p1:
    legend.AddEntry(graph0p1," #tilde{k} = "+str(0.1),"l")
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  c.Write()
  mg.Write()


#def MakeOptMassLowHighVsMassMultigraph():
#def MakeOptMassLowGraph():
  ## make plot of low/high opt. mass window values vs. mass/coupling
  #graph0p01min = TGraph(len(massesC0p01), array.array("f",massesC0p01),array.array("f",optMinMassC0p01))
  #graph0p01min.SetName('optMinMassVsMassK0p01')
  #graph0p01min.SetMarkerColor(2)
  #graph0p01min.SetLineColor(2)
  #graph0p01min.Write()
  #graph0p01max = TGraph(len(massesC0p01), array.array("f",massesC0p01),array.array("f",optMaxMassC0p01))
  #graph0p01max.SetName('optMaxMassVsMassK0p01')
  #graph0p01max.SetMarkerColor(2)
  #graph0p01max.SetLineColor(2)
  #graph0p01max.Write()
  #graph0p05min = TGraph(len(massesC0p05), array.array("f",massesC0p05),array.array("f",optMinMassC0p05))
  #graph0p05min.SetName('optMinMassVsMassK0p05')
  #graph0p05min.SetMarkerColor(4)
  #graph0p05min.SetLineColor(4)
  #graph0p05min.Write()
  #graph0p05max = TGraph(len(massesC0p05), array.array("f",massesC0p05),array.array("f",optMaxMassC0p05))
  #graph0p05max.SetName('optMaxMassVsMassK0p05')
  #graph0p05max.SetMarkerColor(4)
  #graph0p05max.SetLineColor(4)
  #graph0p05max.Write()
  #graph0p1min = TGraph(len(massesC0p1), array.array("f",massesC0p1),array.array("f",optMinMassC0p1))
  #graph0p1min.SetName('optMinMassVsMassK0p1')
  #graph0p1min.SetMarkerColor(8)
  #graph0p1min.SetLineColor(8)
  #graph0p1min.Write()
  #graph0p1max = TGraph(len(massesC0p1), array.array("f",massesC0p1),array.array("f",optMaxMassC0p1))
  #graph0p1max.SetName('optMaxMassVsMassK0p1')
  #graph0p1max.SetMarkerColor(8)
  #graph0p1max.SetLineColor(8)
  #graph0p1max.Write()
  #c = TCanvas()
  #c.SetName('optMinMaxMassWindowsVsMassAllCanvas')
  #c.SetTitle('')
  #c.cd()
  #mg = TMultiGraph()
  #mg.Add(graph0p01min)
  #mg.Add(graph0p01max)
  #mg.Add(graph0p05min)
  #mg.Add(graph0p05max)
  #mg.Add(graph0p1min)
  #mg.Add(graph0p1max)
  #mg.Draw('ap')
  #mg.GetXaxis().SetTitle("Mass [GeV]")
  #mg.GetYaxis().SetTitle("Opt. mass window [GeV]")
  #mg.SetName('optMinMaxMassWindowsVsMassAll')
  #legend = TLegend(0.42,0.71,0.73,0.92)
  #legend.AddEntry(graph0p01min," #tilde{k} = "+str(0.01),"l")
  #legend.AddEntry(graph0p05min," #tilde{k} = "+str(0.05),"l")
  #legend.AddEntry(graph0p1min," #tilde{k} = "+str(0.1),"l")
  #legend.Draw()
  #c.Write()
  #mg.Write()


def MakeOptMassWindowGraphs(coupling,masses,optMinMasses,optMaxMasses,peakMasses,colorIndex,rootFile):
  # filled graphs
  if len(masses) < 1:
    return
  rootFile.cd()
  massErrs = [0]*len(masses)
  massWindErrsUp= [optMaxMass-peakMass for optMaxMass,peakMass in itertools.izip(optMaxMasses,peakMasses)]
  massWindErrsDown = [peakMass-optMinMass for optMinMass,peakMass in itertools.izip(optMinMasses,peakMasses)]
  graphErrs = TGraphAsymmErrors(len(masses),array.array("f",masses),array.array("f",peakMasses),array.array("f",massErrs),array.array("f",massErrs),array.array("f",massWindErrsDown),array.array("f",massWindErrsUp))
  plotname = "optMassWindowsVsMassK"+str(coupling)
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  graphErrs.SetName(savename.Data())
  graphErrs.SetMarkerColor(colorIndex)
  graphErrs.SetLineColor(colorIndex)
  graphErrs.SetFillColor(colorIndex)
  graphErrs.SetFillStyle(3003)
  graphErrs.Write()
  massWindowWidths = [optMaxMass-optMinMass for optMaxMass,optMinMass in itertools.izip(optMaxMasses,optMinMasses)]
  graph = TGraph(len(masses),array.array("f",masses),array.array("f",massWindowWidths))
  plotname = "optMassWindowWidthVsMassK"+str(coupling)
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  graph.SetName(savename.Data())
  graph.SetMarkerColor(colorIndex)
  graph.SetLineColor(colorIndex)
  graph.SetFillColor(colorIndex)
  graph.SetFillStyle(3003)
  graph.Fit('pol1')
  graph.Write()
  myFunc = graph.GetFunction('pol1')
  return myFunc.GetParameter(0),myFunc.GetParameter(1)


def MakeOptMassWindowsVsMassMultiGraph(rootFile):
  rootFile.cd()
  graph0p01Errs = MakeNullPointer(TGraphAsymmErrors)
  graph0p05Errs = MakeNullPointer(TGraphAsymmErrors)
  graph0p1Errs = MakeNullPointer(TGraphAsymmErrors)
  try:
    rootFile.GetObject('optMassWindowsVsMassK0p01',graph0p01Errs)
    rootFile.GetObject('optMassWindowsVsMassK0p05',graph0p05Errs)
    rootFile.GetObject('optMassWindowsVsMassK0p1',graph0p1Errs)
  except LookupError:
    pass
  c = TCanvas()
  c.SetName('optMassWindowsVsMassAllCanvas')
  c.SetTitle('')
  c.cd()
  mg = TMultiGraph()
  if graph0p01Errs:
    mg.Add(graph0p01Errs,'3')
  if graph0p05Errs:
    mg.Add(graph0p05Errs,'3')
  if graph0p1Errs:
    mg.Add(graph0p1Errs,'3')
  mg.Draw('ap')
  mg.GetXaxis().SetTitle("Mass [GeV]")
  mg.GetYaxis().SetTitle("Opt. mass window [GeV]")
  mg.SetName('optMassWindowsVsMassAll')
  legend = TLegend(0.42,0.71,0.73,0.92)
  if graph0p01Errs:
    legend.AddEntry(graph0p01Errs," #tilde{k} = "+str(0.01),"l")
  if graph0p05Errs:
    legend.AddEntry(graph0p05Errs," #tilde{k} = "+str(0.05),"l")
  if graph0p1Errs:
    legend.AddEntry(graph0p1Errs," #tilde{k} = "+str(0.1),"l")
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  c.Write()
  mg.Write()


def MakeOptSSBValuesGraph(coupling,masses,optSSBValues,colorIndex,rootFile):
  # make plot of opt. SSB vs. mass/coupling
  if len(masses) < 1:
    return
  graph = TGraph(len(masses), array.array("f",masses),array.array("f",optSSBValues))
  plotname = "optSSBVsMassK"+str(coupling)
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  graph.SetName(savename.Data())
  graph.SetMarkerColor(colorIndex)
  graph.SetLineColor(colorIndex)
  graph.Write()


def MakeOptSSBValueVsMassMultigraph(rootFile):
  rootFile.cd()
  graph0p01 = MakeNullPointer(TGraph)
  graph0p05 = MakeNullPointer(TGraph)
  graph0p1 = MakeNullPointer(TGraph)
  try:
    rootFile.GetObject('optSSBVsMassK0p01',graph0p01)
    rootFile.GetObject('optSSBVsMassK0p05',graph0p05)
    rootFile.GetObject('optSSBVsMassK0p1',graph0p1)
  except LookupError:
    pass
  c = TCanvas()
  c.SetName('optSSBVsMassAllCanvas')
  c.SetTitle('')
  c.cd()
  mg = TMultiGraph()
  if graph0p01:
    mg.Add(graph0p01,'lp')
  if graph0p05:
    mg.Add(graph0p05,'lp')
  if graph0p1:
    mg.Add(graph0p1,'lp')
  mg.Draw('ap')
  mg.GetXaxis().SetTitle("Mass [GeV]")
  mg.GetYaxis().SetTitle("Opt. S/#sqrt{S+B}")
  mg.SetName('optSSBVsMassAll')
  legend = TLegend(0.42,0.71,0.73,0.92)
  if graph0p01:
    legend.AddEntry(graph0p01," #tilde{k} = "+str(0.01),"l")
  if graph0p05:
    legend.AddEntry(graph0p05," #tilde{k} = "+str(0.05),"l")
  if graph0p1:
    legend.AddEntry(graph0p1," #tilde{k} = "+str(0.1),"l")
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  c.Write()
  mg.Write()


def MakeOptMassWindowSignalBackgroundPlot(rootFile,signalHistogram,backgroundHist,optMassLow,optMassHigh,modelPoint,outputDir):
  # mkdir if not there already
  if not os.path.isdir(outputDir):
    os.mkdir(outputDir)
  rootFile.cd()
  c = TCanvas()
  c.SetName('signalBackgroundOptWindowsCanvas_k'+str(modelPoint.coupling).replace('.','p')+'_m'+str(modelPoint.mass))
  c.SetTitle('')
  c.cd()
  c.SetLogy()
  signalHistogram.SetStats(False)
  signalHistogram.SetLineColor(2)
  signalHistogram.SetMarkerColor(2)
  signalHistogram.Draw()
  signalHistogram.GetXaxis().SetRangeUser(optMassLow-50,optMassHigh+50)
  signalHistogram.Draw()
  backgroundHist.Draw('same')
  lineLow = TLine(optMassLow,0,optMassLow,signalHistogram.GetMaximum())
  lineLow.SetLineColor(4)
  lineLow.Draw()
  lineHigh = TLine(optMassHigh,0,optMassHigh,signalHistogram.GetMaximum())
  lineHigh.SetLineColor(4)
  lineHigh.Draw()
  c.Write()


def MakeOptMassWindowSignalBackgroundImages(rootFile,outputDir,modelPointArray):
  # mkdir if not there already
  if not os.path.isdir(outputDir):
    os.mkdir(outputDir)
  for mp in modelPointArray:
    canvasName = 'signalBackgroundOptWindowsCanvas_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
    c = rootFile.Get(canvasName)
    fileName = 'signalBackgroundOptWindows_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
    savePath = outputDir+'/'+fileName
    c.Print(savePath+'.pdf')
    c.Print(savePath+'.eps')
    ConvertToPng(savePath)


def MakeOptGraphImages(rootFile,outputDir,modelPointArray, plotSRootBCurve):
  # mkdir if not there already
  if not os.path.isdir(outputDir):
    os.mkdir(outputDir)
  gStyle.SetOptFit(1)
  gStyle.SetOptStat(0)
  gStyle.SetOptTitle(0)
  c = TCanvas("c", "c",0,0,600,600)
  c.SetHighLightColor(2)
  c.Range(200,-3.536669,3637.5,-2.469521)
  c.SetFillColor(0)
  c.SetBorderMode(0)
  c.SetBorderSize(2)
  c.SetTickx(1)
  c.SetTicky(1)
  c.SetLeftMargin(0.18)
  c.SetRightMargin(0.04)
  c.SetTopMargin(0.05)
  c.SetBottomMargin(0.15)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.cd()
  for mp in modelPointArray:
    graphName = 'ssbOpt_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
    #print 'Get graph',graphName,'from',rootFile.GetName()
    optGraph = rootFile.Get(graphName)
    optPtGraphName = 'ssbOptPoint_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
    optPtGraph = rootFile.Get(optPtGraphName)
    optPtGraph.SetMarkerColor(2)
    optPtGraph.SetLineColor(2)
    optPtGraph.SetMarkerStyle(3)
    optPtGraph.SetMarkerSize(1.4)
    mg = TMultiGraph()
    drawOpt = 'p'
    mg.Add(optGraph,drawOpt)
    mg.Add(optPtGraph,drawOpt)
    leg = TLegend(0.75,0.65,0.95,0.86)
    leg.SetBorderSize(0)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(2)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    legDrawOpt = 'p'
    leg.AddEntry(optGraph,"s/#sqrt{s+b}",legDrawOpt)
    leg.AddEntry(optPtGraph,"Opt. s/#sqrt{s+b}",legDrawOpt)
    if plotSRootBCurve:
      srbGraphName = 'srootbOpt_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
      srbOptGraph = rootFile.Get(srbGraphName)
      srbOptGraph.SetLineColor(4)
      srbOptGraph.SetMarkerColor(4)
      srbOptGraph.SetMarkerStyle(3)
      srbOptGraph.SetMarkerSize(0.5)
      srbOptPtGraphName = 'srootbOptPoint_k'+str(mp.coupling).replace('.','p')+'_m'+str(mp.mass)
      srbOptPtGraph = rootFile.Get(srbOptPtGraphName)
      srbOptPtGraph.SetLineColor(8)
      srbOptPtGraph.SetMarkerColor(8)
      srbOptPtGraph.SetMarkerStyle(3)
      srbOptPtGraph.SetMarkerSize(1.4)
      mg.Add(srbOptGraph,drawOpt)
      mg.Add(srbOptPtGraph,drawOpt)
      leg.AddEntry(srbOptGraph,"s/#sqrt{b}",legDrawOpt)
      leg.AddEntry(srbOptPtGraph,"Opt. s/#sqrt{b}",legDrawOpt)
      lineOptSRB = TLine(srbOptPtGraph.GetX()[0],0,srbOptPtGraph.GetX()[0],srbOptPtGraph.GetY()[0])
      lineOptSRB.SetVertical()
      lineOptSRB.SetLineColor(8)
      lineOptSRB.Draw()
      lineOptSSB = TLine(optPtGraph.GetX()[0],0,optPtGraph.GetX()[0],optPtGraph.GetY()[0])
      lineOptSSB.SetVertical()
      lineOptSSB.SetLineColor(2)
      lineOptSSB.Draw()
    graphNameTemplate = 'ssbOpt_k{coupling}_m{mass}'
    graphName = graphNameTemplate.format(coupling=str(mp.coupling).replace('.','p'),mass=str(mp.mass))
    mg.Draw('a')
    #mg.GetXaxis().SetTitle('Window Size (GeV)')
    #mg.GetYaxis().SetTitle('Opt. Var.')
    leg.Draw()
    #
    savePath = outputDir+'/'+graphName
    c.Print(savePath+'.pdf')
    c.Print(savePath+'.eps')
    #c.Print(savePath+'.C')
    # png output looks strange, so convert from pdf instead
    ConvertToPng(savePath)
    # eps too big
    os.unlink(savePath+'.eps')


def MakeSmoothedMassWindowVsMassImages(rootFile,imageDir):
  #TODO Make this take the model array like the above...
  # mkdir if not there already
  if not os.path.isdir(imageDir):
    os.mkdir(imageDir)
  c = TCanvas("c", "c",0,0,600,600)
  gStyle.SetOptFit(1)
  gStyle.SetOptStat(0)
  gStyle.SetOptTitle(0)
  c.SetHighLightColor(2)
  c.Range(200,-3.536669,3637.5,-2.469521)
  c.SetFillColor(0)
  c.SetBorderMode(0)
  c.SetBorderSize(2)
  c.SetTickx(1)
  c.SetTicky(1)
  c.SetLeftMargin(0.18)
  c.SetRightMargin(0.04)
  c.SetTopMargin(0.05)
  c.SetBottomMargin(0.15)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  #
  symmWindowOptGraphk0p01 = MakeNullPointer(TGraph)
  symmWindowOptGraphk0p05 = MakeNullPointer(TGraph)
  symmWindowOptGraphk0p1 = MakeNullPointer(TGraph)
  try:
    rootFile.GetObject('optSmoothedWindowVsMassK0p01',symmWindowOptGraphk0p01)
    rootFile.GetObject('optSmoothedWindowVsMassK0p05',symmWindowOptGraphk0p05)
    rootFile.GetObject('optSmoothedWindowVsMassK0p1',symmWindowOptGraphk0p1)
  except:
    pass
  test = TH1F("test","test",10,750,3500)
  test.SetMinimum(600)
  test.SetMaximum(4100)
  test.SetStats(0)
  test.SetLineStyle(0)
  test.SetLineWidth(2)
  test.SetMarkerStyle(20)
  test.SetMarkerSize(0.8)
  test.GetXaxis().SetTitle("M_{1} [GeV]")
  test.GetXaxis().SetLabelFont(42)
  test.GetXaxis().SetLabelOffset(0.007)
  test.GetXaxis().SetTitleOffset(1.2)
  test.GetXaxis().SetTitleFont(42)
  test.GetYaxis().SetTitle("Mass window width [GeV]")
  test.GetYaxis().SetLabelFont(42)
  test.GetYaxis().SetLabelOffset(0.007)
  test.GetYaxis().SetTitleOffset(1.5)
  test.GetYaxis().SetTitleFont(42)
  test.Draw()
  #
  mg = TMultiGraph()
  if symmWindowOptGraphk0p01:
    mg.Add(symmWindowOptGraphk0p01,"p")
  if symmWindowOptGraphk0p05:
    mg.Add(symmWindowOptGraphk0p05,"p")
  if symmWindowOptGraphk0p1:
    mg.Add(symmWindowOptGraphk0p1,"p")
  mg.Draw("l")
  #
  leg = TLegend(0.18,0.65,0.5,0.86)
  leg.SetBorderSize(0)
  leg.SetLineColor(1)
  leg.SetLineStyle(1)
  leg.SetLineWidth(2)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  legDrawOpt = 'lp'
  entry=leg.AddEntry(symmWindowSmoothedraphk0p01,"Smoothed Symm. Windows #tilde{k} = 0.01",legDrawSmoothed
  entry.SetLineColor(2)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(2)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  entry=leg.AddEntry(symmWindowSmoothedraphk0p05,"Smoothed Symm. Windows #tilde{k} = 0.05",legDrawSmoothed
  entry.SetLineColor(4)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(4)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  entry=leg.AddEntry(symmWindowSmoothedraphk0p1,"Smoothed Symm. Windows #tilde{k} = 0.1",legDrawSmoothed
  entry.SetLineColor(8)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(8)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  leg.Draw()
  #
  test__6 = TH1F("test__6","test",10,750,3500)
  test__6.SetMinimum(600)
  test__6.SetMaximum(4100)
  test__6.SetDirectory(0)
  test__6.SetStats(0)
  test__6.SetLineStyle(0)
  test__6.SetLineWidth(2)
  test__6.SetMarkerStyle(20)
  test__6.SetMarkerSize(0.8)
  test__6.GetXaxis().SetTitle("M_{1} [GeV]")
  test__6.GetXaxis().SetLabelFont(42)
  test__6.GetXaxis().SetLabelOffset(0.007)
  test__6.GetXaxis().SetTitleOffset(1.2)
  test__6.GetXaxis().SetTitleFont(42)
  test__6.GetYaxis().SetTitle("Mass Window width [GeV]")
  test__6.GetYaxis().SetLabelFont(42)
  test__6.GetYaxis().SetLabelOffset(0.007)
  test__6.GetYaxis().SetTitleOffset(1.5)
  test__6.GetYaxis().SetTitleFont(42)
  test__6.Draw("sameaxis")
  fileName = 'smoothedMassWindowsVsMass'
  savePath = imageDir+'/'+fileName
  c.Print(savePath+'.pdf')
  c.Print(savePath+'.eps')
  c.Print(savePath+'.C')
  ConvertToPng(savePath)


def MakeOptMassWindowVsMassImages(rootFile,imageDir):
  #TODO Make this take the model array like the above...
  # mkdir if not there already
  if not os.path.isdir(imageDir):
    os.mkdir(imageDir)
  c = TCanvas("c", "c",0,0,600,600)
  gStyle.SetOptFit(1)
  gStyle.SetOptStat(0)
  gStyle.SetOptTitle(0)
  c.SetHighLightColor(2)
  c.Range(200,-3.536669,3637.5,-2.469521)
  c.SetFillColor(0)
  c.SetBorderMode(0)
  c.SetBorderSize(2)
  c.SetTickx(1)
  c.SetTicky(1)
  c.SetLeftMargin(0.18)
  c.SetRightMargin(0.04)
  c.SetTopMargin(0.05)
  c.SetBottomMargin(0.15)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  #
  symmWindowOptGraphk0p01 = MakeNullPointer(TGraphAsymmErrors)
  symmWindowOptGraphk0p05 = MakeNullPointer(TGraphAsymmErrors)
  symmWindowOptGraphk0p1 = MakeNullPointer(TGraphAsymmErrors)
  try:
    rootFile.GetObject('optMassWindowsVsMassK0p01',symmWindowOptGraphk0p01)
    rootFile.GetObject('optMassWindowsVsMassK0p05',symmWindowOptGraphk0p05)
    rootFile.GetObject('optMassWindowsVsMassK0p1',symmWindowOptGraphk0p1)
  except:
    pass
  if symmWindowOptGraphk0p01:
    symmWindowOptGraphk0p01.SetFillStyle(3002)
  if symmWindowOptGraphk0p05:
    symmWindowOptGraphk0p05.SetFillStyle(3007)
  if symmWindowOptGraphk0p1:
    symmWindowOptGraphk0p1.SetFillStyle(3006)
  test = TH1F("test","test",10,750,3500)
  test.SetMinimum(600)
  test.SetMaximum(4100)
  test.SetStats(0)
  test.SetLineStyle(0)
  test.SetLineWidth(2)
  test.SetMarkerStyle(20)
  test.SetMarkerSize(0.8)
  test.GetXaxis().SetTitle("M_{1} [GeV]")
  test.GetXaxis().SetLabelFont(42)
  test.GetXaxis().SetLabelOffset(0.007)
  test.GetXaxis().SetTitleOffset(1.2)
  test.GetXaxis().SetTitleFont(42)
  test.GetYaxis().SetTitle("Mass window [GeV]")
  test.GetYaxis().SetLabelFont(42)
  test.GetYaxis().SetLabelOffset(0.007)
  test.GetYaxis().SetTitleOffset(1.5)
  test.GetYaxis().SetTitleFont(42)
  test.GetZaxis().SetLabelFont(42)
  test.GetZaxis().SetLabelOffset(0.007)
  test.GetZaxis().SetLabelSize(0.05)
  test.GetZaxis().SetTitleSize(0.06)
  test.GetZaxis().SetTitleOffset(1.1)
  test.GetZaxis().SetTitleFont(42)
  test.Draw()
  #
  mg = TMultiGraph()
  if symmWindowOptGraphk0p01:
    mg.Add(symmWindowOptGraphk0p01,"3")
  if symmWindowOptGraphk0p05:
    mg.Add(symmWindowOptGraphk0p05,"3")
  if symmWindowOptGraphk0p1:
    mg.Add(symmWindowOptGraphk0p1,"3")
  mg.Draw("l")
  #
  leg = TLegend(0.18,0.65,0.5,0.86)
  leg.SetBorderSize(0)
  leg.SetLineColor(1)
  leg.SetLineStyle(1)
  leg.SetLineWidth(2)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  legDrawOpt = 'f'
  entry=leg.AddEntry(symmWindowOptGraphk0p01,"Opt. Symm. Windows #tilde{k} = 0.01",legDrawOpt)
  entry.SetLineColor(2)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(2)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  entry=leg.AddEntry(symmWindowOptGraphk0p05,"Opt. Symm. Windows #tilde{k} = 0.05",legDrawOpt)
  entry.SetLineColor(4)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(4)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  entry=leg.AddEntry(symmWindowOptGraphk0p1,"Opt. Symm. Windows #tilde{k} = 0.1",legDrawOpt)
  entry.SetLineColor(8)
  entry.SetLineStyle(2)
  entry.SetLineWidth(4)
  entry.SetMarkerColor(8)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  leg.Draw()
  #
  test__6 = TH1F("test__6","test",10,750,3500)
  test__6.SetMinimum(600)
  test__6.SetMaximum(4100)
  test__6.SetDirectory(0)
  test__6.SetStats(0)
  test__6.SetLineStyle(0)
  test__6.SetLineWidth(2)
  test__6.SetMarkerStyle(20)
  test__6.SetMarkerSize(0.8)
  test__6.GetXaxis().SetTitle("M_{1} [GeV]")
  test__6.GetXaxis().SetLabelFont(42)
  test__6.GetXaxis().SetLabelOffset(0.007)
  test__6.GetXaxis().SetTitleOffset(1.2)
  test__6.GetXaxis().SetTitleFont(42)
  test__6.GetYaxis().SetTitle("Mass Window [GeV]")
  test__6.GetYaxis().SetLabelFont(42)
  test__6.GetYaxis().SetLabelOffset(0.007)
  test__6.GetYaxis().SetTitleOffset(1.5)
  test__6.GetYaxis().SetTitleFont(42)
  test__6.GetZaxis().SetLabelFont(42)
  test__6.GetZaxis().SetLabelOffset(0.007)
  test__6.GetZaxis().SetLabelSize(0.05)
  test__6.GetZaxis().SetTitleSize(0.06)
  test__6.GetZaxis().SetTitleOffset(1.1)
  test__6.GetZaxis().SetTitleFont(42)
  test__6.Draw("sameaxis")
  fileName = 'optimizedMassWindowsVsMass'
  savePath = imageDir+'/'+fileName
  c.Print(savePath+'.pdf')
  c.Print(savePath+'.eps')
  c.Print(savePath+'.C')
  ConvertToPng(savePath)


def MakeOptSSBVsMassImages(rootFile, outputDir):
  symmWindowOptGraphk0p01 = MakeNullPointer(TGraph)
  symmWindowOptGraphk0p05 = MakeNullPointer(TGraph)
  symmWindowOptGraphk0p1 = MakeNullPointer(TGraph)
  try:
    rootFile.GetObject('optSSBVsMassK0p01',symmWindowOptGraphk0p01)
    rootFile.GetObject('optSSBVsMassK0p05',symmWindowOptGraphk0p05)
    rootFile.GetObject('optSSBVsMassK0p1',symmWindowOptGraphk0p1)
  except LookupError:
    pass
  #
  c = TCanvas("c", "c",0,0,600,600)
  gStyle.SetOptFit(1)
  gStyle.SetOptStat(0)
  gStyle.SetOptTitle(0)
  c.SetHighLightColor(2)
  c.Range(200,-3.536669,3637.5,-2.469521)
  c.SetFillColor(0)
  c.SetBorderMode(0)
  c.SetBorderSize(2)
  c.SetTickx(1)
  c.SetTicky(1)
  c.SetLeftMargin(0.18)
  c.SetRightMargin(0.04)
  c.SetTopMargin(0.05)
  c.SetBottomMargin(0.15)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.SetFrameFillStyle(0)
  c.SetFrameLineWidth(2)
  c.SetFrameBorderMode(0)
  c.SetLogy()
  #
  minYrange = 1e-1
  maxYrange = 8
  test = TH1F("test","test",10,750,3500)
  test.SetMinimum(minYrange)
  test.SetMaximum(maxYrange)
  test.SetStats(0)
  test.SetLineStyle(0)
  test.SetLineWidth(2)
  test.SetMarkerStyle(20)
  test.SetMarkerSize(0.8)
  test.GetXaxis().SetTitle("M_{1} [GeV]")
  test.GetXaxis().SetLabelFont(42)
  test.GetXaxis().SetLabelOffset(0.007)
  test.GetXaxis().SetTitleOffset(1.2)
  test.GetXaxis().SetTitleFont(42)
  test.GetYaxis().SetTitle("Opt. S/#sqrt{S+B}")
  test.GetYaxis().SetLabelFont(42)
  test.GetYaxis().SetLabelOffset(0.007)
  test.GetYaxis().SetTitleOffset(1.5)
  test.GetYaxis().SetTitleFont(42)
  test.GetZaxis().SetLabelFont(42)
  test.GetZaxis().SetLabelOffset(0.007)
  test.GetZaxis().SetLabelSize(0.05)
  test.GetZaxis().SetTitleSize(0.06)
  test.GetZaxis().SetTitleOffset(1.1)
  test.GetZaxis().SetTitleFont(42)
  test.Draw()
  #
  mg = TMultiGraph()
  drawOpt = 'pl'
  if symmWindowOptGraphk0p01:
    mg.Add(symmWindowOptGraphk0p01,drawOpt)
  if symmWindowOptGraphk0p05:
    mg.Add(symmWindowOptGraphk0p05,drawOpt)
  if symmWindowOptGraphk0p1:
    mg.Add(symmWindowOptGraphk0p1,drawOpt)
  mg.Draw("l")
  #
  leg = TLegend(0.6,0.7,0.8,0.9)
  leg.SetBorderSize(0)
  leg.SetLineColor(1)
  leg.SetLineStyle(1)
  leg.SetLineWidth(2)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  legDrawOpt = 'lp'
  if symmWindowOptGraphk0p01:
    entry=leg.AddEntry(symmWindowOptGraphk0p01,"Symm. Window #tilde{k} = 0.01",legDrawOpt)
    entry.SetLineColor(2)
    entry.SetLineStyle(2)
    entry.SetLineWidth(4)
    entry.SetMarkerColor(2)
    entry.SetMarkerStyle(21)
    entry.SetMarkerSize(1)
  if symmWindowOptGraphk0p05:
    entry=leg.AddEntry(symmWindowOptGraphk0p05,"Symm. Window #tilde{k} = 0.05",legDrawOpt)
    entry.SetLineColor(4)
    entry.SetLineStyle(2)
    entry.SetLineWidth(4)
    entry.SetMarkerColor(4)
    entry.SetMarkerStyle(21)
    entry.SetMarkerSize(1)
  if symmWindowOptGraphk0p1:
    entry=leg.AddEntry(symmWindowOptGraphk0p1,"Symm. Window #tilde{k} = 0.1",legDrawOpt)
    entry.SetLineColor(8)
    entry.SetLineStyle(2)
    entry.SetLineWidth(4)
    entry.SetMarkerColor(8)
    entry.SetMarkerStyle(21)
    entry.SetMarkerSize(1)
  leg.Draw()
  #
  test__6 = TH1F("test__6","test",10,750,3500)
  test__6.SetMinimum(minYrange)
  test__6.SetMaximum(maxYrange)
  test__6.SetDirectory(0)
  test__6.SetStats(0)
  test__6.SetLineStyle(0)
  test__6.SetLineWidth(2)
  test__6.SetMarkerStyle(20)
  test__6.SetMarkerSize(0.8)
  test__6.GetXaxis().SetTitle("M_{1} [GeV]")
  test__6.GetXaxis().SetLabelFont(42)
  test__6.GetXaxis().SetLabelOffset(0.007)
  test__6.GetXaxis().SetTitleOffset(1.2)
  test__6.GetXaxis().SetTitleFont(42)
  test__6.GetYaxis().SetTitle("Opt. S/#sqrt{S+B}")
  test__6.GetYaxis().SetLabelFont(42)
  test__6.GetYaxis().SetLabelOffset(0.007)
  test__6.GetYaxis().SetTitleOffset(1.5)
  test__6.GetYaxis().SetTitleFont(42)
  test__6.GetZaxis().SetLabelFont(42)
  test__6.GetZaxis().SetLabelOffset(0.007)
  test__6.GetZaxis().SetLabelSize(0.05)
  test__6.GetZaxis().SetTitleSize(0.06)
  test__6.GetZaxis().SetTitleOffset(1.1)
  test__6.GetZaxis().SetTitleFont(42)
  test__6.Draw("sameaxis")
  #
  fileName = 'optSSBVsMassAll'
  savePath = outputDir+'/'+fileName
  c.Print(savePath+'.pdf')
  c.Print(savePath+'.eps')
  c.Print(savePath+'.C')
  ConvertToPng(savePath)


def PlotAllBands(modelPointArrays, lumi, rootFile, imageDir):
  print 'PlotAllBands'
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  limitExpArrs = []
  limitObsArrs = []
  totalXSecArrs = []
  couplings = []
  for mpArray in modelPointArrays:
    errMass, errUp, errDn, err2Up, err2Dn, masses, expUp, expDn, exp2Up, exp2Dn, limitExp, limitObs, totalXSec = makeArrays(mpArray)
    if len(errUp)==0:
      continue
    massesArrs.append(masses)
    limitExpArrs.append(limitExp)
    limitObsArrs.append(limitObs)
    totalXSecArrs.append(totalXSec)
    couplings.append(mpArray[0].coupling)

  if len(massesArrs)==0:
    return

  rootFile.cd()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  c.SetLogy()
  c.SetRightMargin(0.04)
  # turn the arrays into graphs
  limitObsGraphs = []
  limitExpGraphs = []
  thGraphs = []
  for massesArr, limitExpArr, limitObsArr, thArr in itertools.izip(massesArrs,limitExpArrs,limitObsArrs,totalXSecArrs):
    limitObsGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",limitObsArr)))
    limitExpGraphs.append(TGraph(len(massesArr), array.array("f",massesArr), array.array("f",limitExpArr)))
    thGraphs.append(TGraph(len(massesArr), array.array("f",massesArr), array.array("f",thArr)))

  SetCustomGStyle()

  # Multigraph
  mg = TMultiGraph()
  #legend = TLegend(0.42,0.71,0.73,0.92)
  legend = TLegend(0.6,0.71,0.9,0.92)
  colorIndex = 2 
  for limitObsGraph, limitExpGraph, thGraph, coupling in itertools.izip(limitObsGraphs,limitExpGraphs,thGraphs,couplings):
    limitExpGraph.SetMarkerSize(2)
    limitExpGraph.SetMarkerColor(colorIndex)
    limitExpGraph.SetLineColor(colorIndex)
    limitExpGraph.SetLineWidth(4)
    limitExpGraph.SetLineStyle(2)
    limitExpGraph.SetName("limitExpGraph_k"+str(coupling))
    limitExpGraph.Write()
    #
    limitObsGraph.SetMarkerStyle(21)
    limitObsGraph.SetMarkerSize(1.0)
    limitObsGraph.SetMarkerColor(colorIndex)
    limitObsGraph.SetLineColor(colorIndex)
    limitObsGraph.SetLineWidth(4)
    limitObsGraph.SetName("limitObsGraph_k"+str(coupling))
    limitObsGraph.Write()
    #
    thGraph.SetLineStyle(7)
    thGraph.SetLineColor(colorIndex)
    thGraph.SetLineWidth(4)
    thGraph.SetName("thGraph_k"+str(coupling))
    thGraph.Write()
    #
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(limitObsGraph,"LP")
    legend.AddEntry(limitObsGraph,"95%"+" CL limit #tilde{k} = "+str(coupling),"LP")
    mg.Add(thGraph,"L")
    legend.AddEntry(thGraph,"G_{KK} #tilde{k} = "+str(coupling),"l")
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1


  # To make the separate graphs share the same axes
  h = TH1F("test","test",10,750,3500)
  h.SetStats(False)
  #h.GetYaxis().SetRangeUser(3e-4,3e-3)
  h.GetYaxis().SetRangeUser(2.4e-4,3e-3)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("RS graviton #sigma [pb]     ")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw("L")

  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()

  # CMS
  pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  pt.SetName("CMS Preliminary")
  pt.SetBorderSize(1)
  pt.SetLineColor(0)
  pt.SetFillColor(0)
  pt.SetTextSize(0.0354545)
  text = pt.AddText("CMS Preliminary")
  pt.Draw()
  # lumi
  pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  pt2.SetFillColor(0)
  pt2.SetBorderSize(1)
  pt2.SetLineColor(0)
  pt2.SetTextSize(0.0354545)
  text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  pt2.Draw()

  gPad.RedrawAxis()

  plotname = "allLimits.pdf"
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = "allLimits.C"
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allLimits.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  # write
  mg.SetName("limitsAllCouplingsMultiGraph")
  mg.Write()


def PlotBands(modelPointArray, lumi, rootFile, imageDir):
  #print '----- PlotBands -----'
  #print 'Coupling:',modelPointArray[0].coupling
  # fill arrays
  errMass, errUp, errDn, err2Up, err2Dn, masses, expUp, expDn, exp2Up, exp2Dn, limitExp, limitObs, totalXSec = makeArrays(modelPointArray)

  if len(errUp)==0:
    return

  #turn the arrays into graphs
  #canvas = TCanvas("canvas","canvas",100,200,500,500)
  #canvas.cd()
  #canvas.SetLeftMargin(0.139262)
  #canvas.SetRightMargin(0.0604027)
  #canvas.SetTopMargin(0.0804196)
  #canvas.SetBottomMargin(0.14)
  #g_up = TGraph(len(modelPointArray), array.array("f",masses),array.array("f",expUp));
  #g_dn = TGraph(len(modelPointArray), array.array("f",masses), array.array("f",expDn));
  #g_2up = TGraph(len(modelPointArray), array.array("f",masses), array.array("f",exp2Up));
  #g_2dn = TGraph(len(modelPointArray), array.array("f",masses), array.array("f",exp2Dn));
  #g_exp = TGraph(len(modelPointArray), array.array("f",masses), array.array("f",limitExp));
  #g_up.SetMarkerStyle(20)
  #g_up.Draw("AP")
  #g_dn.Draw("P")
  #g_2up.Draw("P")
  #g_2dn.Draw("P")
  #g_exp.Draw("P")
  #smooth them (didn't work)
  #gs = TGraphSmooth("normnal");
  #gs_up = gs.SmoothSuper(g_up,"",3,0);
  #gs_dn = gs.SmoothSuper(g_dn,"",3,0);
  #gs_2up = gs.SmoothSuper(g_2up,"",3,0);
  #gs_2dn = gs.SmoothSuper(g_2dn,"",3,0);
  #gs_exp =  gs.SmoothSuper(g_exp,"",3,0);
  ##gs_up = gs.SmoothKern(g_up);
  ##gs_dn = gs.SmoothKern(g_dn);
  ##gs_2up = gs.SmoothKern(g_2up);
  ##gs_2dn = gs.SmoothKern(g_2dn);
  ##gs_exp =  gs.SmoothKern(g_exp);
  ##gs_up = gs.SmoothLowess(g_up,"",1);
  ##gs_dn = gs.SmoothLowess(g_dn,"",1);
  ##gs_2up = gs.SmoothLowess(g_2up,"",1);
  ##gs_2dn = gs.SmoothLowess(g_2dn,"",1);
  ##gs_exp =  gs.SmoothLowess(g_exp,"",1);
  # gs_up.SetLineColor(kRed);
  # gs_up.Draw("L");
  # gs_dn.Draw("L");
  # gs_2up.Draw("L");
  # gs_2dn.Draw("L");
  # gs_exp.Draw("L");
  #   //fit them instead (pol8 seems to work), not for CLs, try pol2
  #   gs_up = new TF1("gs_up","pol2",500,2000);
  #   gs_dn = new TF1("gs_dn","pol2",500,2000);
  #   gs_2up = new TF1("gs_2up","pol2",500,2000);
  #   gs_dn = new TF1("gs_2dn","pol2",500,2000);
  #   //gs_exp = new TF1("gs_exp","pol8",500,2000);
  #   g_up->Fit(gs_up,"QS");
  #   g_dn->Fit(gs_dn,"QS");
  #   g_2up->Fit(gs_2up,"QS");
  #   g_2dn->Fit(gs_2dn,"QS");
  #   //g_exp->Fit(gs_exp,"QS");
  #   gs_up->Draw("same");
  #   gs_dn->Draw("same");
  #   gs_2up->Draw("same");
  #   gs_2dn->Draw("same");
  #   //gs_exp->Draw("same");

  #problem with CLs smoothing, turn it off for now.
  #  for (int i = 0; i< 5; i++) {
  #   err_up[i] = gs_up->Eval(MASSES[i]) - gs_exp->Eval(MASSES[i]);
  #   err_2up[i] = gs_2up->Eval(MASSES[i]) - gs_exp->Eval(MASSES[i]);
  #   if (i>10) {
  #   err_dn[i] = -( gs_dn->Eval(MASSES[i]) - gs_exp->Eval(MASSES[i]) );
  #   err_2dn[i] = -( gs_2dn->Eval(MASSES[i]) - gs_exp->Eval(MASSES[i]) );
  #   } //bottom part got smoothed away otherwise...
  #   if (MASSES[i]>1300 && (err_up[i] < 0.001) ) {
  #     err_up[i] = 0;
  #     err_2up[i] = 0;
  #     err_dn[i] = 0;
  #     err_2dn[i] = 0;
  #  }
  #  }
  #cout<<gs_up->Eval(700)<<endl;;

  rootFile.cd()

  g01_obs = TGraph(len(masses), array.array("f",masses),array.array("f",limitObs));
  g01_exp = TGraph(len(masses), array.array("f",masses), array.array("f",limitExp));

  g01_exp_2s = TGraphAsymmErrors(len(masses), array.array("f",masses), array.array("f",limitExp),
                                 array.array("f",errMass),array.array("f",errMass),
                                 array.array("f",err2Dn),array.array("f",err2Up))
  g01_exp_1s = TGraphAsymmErrors(len(masses), array.array("f",masses), array.array("f",limitExp),
                                 array.array("f",errMass),array.array("f",errMass),
                                 array.array("f",errDn),array.array("f",errUp))


  g01_theory = TGraph(len(masses), array.array("f",masses), array.array("f",totalXSec))

  SetCustomGStyle()

  g01_exp_2s.SetTitle("")
  g01_exp_2s.GetYaxis().SetLabelSize(0.03)
  g01_exp_2s.GetXaxis().SetLabelSize(0.03)
  g01_exp_2s.GetYaxis().SetTitle("RS graviton #sigma [pb]     ")
  g01_exp_2s.GetYaxis().SetTitleOffset(1.3)
  g01_exp_2s.GetXaxis().SetTitleOffset(1.0)
  g01_exp_2s.GetXaxis().SetTitle("M_{1} [GeV]")

  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  c.SetLogy()
  c.SetRightMargin(0.04)

  g01_exp_2s.SetMarkerStyle(20)
  g01_exp_2s.SetMarkerSize(2)
  g01_exp_2s.SetMarkerColor(2)
  g01_exp_2s.SetLineColor(5)
  g01_exp_2s.SetLineWidth(4)
  #  g01_exp_2s.SetLineStyle(2)
  g01_exp_2s.SetFillColor(5)
  #g01_exp_2s.GetXaxis().SetRangeUser(MASSES[0],3000)
  #g01_exp_2s.GetXaxis().SetRangeUser(750,3000)
  g01_exp_2s.GetXaxis().SetNdivisions(510)
  #g01_exp_2s.GetYaxis().SetRangeUser(1.e-4,totalXSec[0]*10.)
  #g01_exp_2s.GetYaxis().SetRangeUser(2.e-4,8.e-3)

  #g01_exp_2s.Draw("3AC")

  #1s bands
  g01_exp_1s.SetMarkerSize(2)
  g01_exp_1s.SetMarkerColor(2)
  g01_exp_1s.SetLineColor(3)
  g01_exp_1s.SetLineWidth(4)
  #  g01_exp_1s.SetLineStyle(2)
  g01_exp_1s.SetFillColor(3)
  #g01_exp_1s.GetXaxis().SetRangeUser(MASSES[0],3000)
  #g01_exp_1s.GetXaxis().SetRangeUser(750,3000)

  #g01_exp_1s.Draw("3C")

  #expected dashed line
  #g01_exp_smooth = TGraphSmooth("normal")
  #   gr_exp_out = g01_exp_smooth.SmoothSuper(g01_exp,"",0,0)
  #   gr_exp_out.SetMarkerSize(2)
  #   gr_exp_out.SetMarkerColor(2)
  #   gr_exp_out.SetLineColor(2)
  #   gr_exp_out.SetLineWidth(4)
  #   gr_exp_out.SetLineStyle(2)
  #gr_exp_out.Draw("C")

  g01_exp.SetMarkerSize(2)
  g01_exp.SetMarkerColor(2)
  g01_exp.SetLineColor(2)
  g01_exp.SetLineWidth(4)
  g01_exp.SetLineStyle(2)
  #g01_exp.Draw("C")


  g01_obs.SetMarkerStyle(21)
  g01_obs.SetMarkerSize(1)
  g01_obs.SetMarkerColor(1)
  g01_obs.SetLineColor(1)
  g01_obs.SetLineWidth(4)

  #g01_obs.Draw("sameL")
  ###XXX SIC
  #g01_theory.SetMarkerSize(2)
  #g01_theory.SetMarkerColor(3)
  g01_theory.SetLineStyle(9) #big dashes, small spaces
  g01_theory.SetLineColor(4)
  g01_theory.SetLineWidth(4)
  ###XXX SIC
  #  g01_theory.SetMarkerStyle(20)
  #g01_theory.Draw("sameL")
  #g01_theory.GetXaxis().SetRangeUser(MASSES[0],3000)
  #g01_theory.GetXaxis().SetRangeUser(750,3000)
  #g01_theory.GetYaxis().SetRangeUser(1.e-5,0.6)


  # Multigraph
  mg = TMultiGraph()
  #mg.Add(g_up,"P")
  #mg.Add(g_dn,"P")
  #mg.Add(g_2up,"P")
  #mg.Add(g_2dn,"P")
  #mg.Add(g_exp,"P")
  mg.Add(g01_exp_2s,"3C")
  mg.Add(g01_exp_1s,"3C")
  mg.Add(g01_exp,"C")
  mg.Add(g01_obs,"LP")
  mg.Add(g01_theory,"L")
  mg.Draw("AL")
  # To have diff. x-axes for each limit graph (TMultiGraph doesn't exactly respect zoomed x-axis)
  #mg.GetYaxis().SetRangeUser(1.e-4,totalXSec[0]*10.)
  #mg.GetYaxis().SetRangeUser(2e-4,2e-2)
  mg.GetYaxis().SetRangeUser(limitExp[len(limitExp)-1]/2,totalXSec[0]*2)
  #mg.GetXaxis().SetRangeUser(totalXSec[0],totalXSec[len(totalXSec)-1])
  mg.GetXaxis().SetRangeUser(masses[0],masses[len(masses)-1])
  #mg.GetXaxis().SetRangeUser(750,3500)
  mg.GetXaxis().SetTitle("M_{1} [GeV]")
  mg.GetYaxis().SetTitle("RS graviton #sigma [pb]     ")
  mg.GetXaxis().SetLabelFont(42)
  mg.GetYaxis().SetLabelFont(42)
  mg.GetXaxis().SetLabelSize(0.04)
  mg.GetYaxis().SetLabelSize(0.04)
  mg.GetYaxis().SetTitleOffset(1.5)
  mg.GetXaxis().SetTitleOffset(1.2)
  mg.GetXaxis().SetTitleSize(0.04)
  mg.GetYaxis().SetTitleSize(0.04)
  #
  ## To make the separate graphs share the same axes
  #h = TH1F("test","test",10,750,3500)
  #h.SetStats(False)
  #h.GetYaxis().SetRangeUser(2e-4,2e-2)
  #h.GetXaxis().SetTitle("M_{1} [GeV]")
  #h.GetYaxis().SetTitle("RS graviton #sigma [pb]     ")
  #h.GetXaxis().SetLabelFont(42)
  #h.GetYaxis().SetLabelFont(42)
  #h.GetXaxis().SetLabelSize(0.04)
  #h.GetYaxis().SetLabelSize(0.04)
  #h.GetYaxis().SetTitleOffset(1.5)
  #h.GetXaxis().SetTitleOffset(1.2)
  #h.GetXaxis().SetTitleSize(0.04)
  #h.GetYaxis().SetTitleSize(0.04)
  #h.Draw()
  #mg.Draw("L")
  #

  titlename = "G_{KK} #tilde{k} = "+str(modelPointArray[0].coupling)
  legend = TLegend(.42,0.71,.73,.88)
  #  legend.AddEntry((TObject*)0,title,"")
  #legend.AddEntry(gr_exp_out ,"median expected","l")
  legend.AddEntry(g01_exp,"median expected","l")
  legend.AddEntry(g01_exp_1s ,"68% expected","lf")
  legend.AddEntry(g01_exp_2s ,"95% expected","lf")
  legend.AddEntry(g01_theory ,titlename,"l")
  legend.AddEntry(g01_obs ,"95% CL limit","l")
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()

  # CMS
  pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  pt.SetName("CMS Preliminary")
  pt.SetBorderSize(1)
  pt.SetLineColor(0)
  pt.SetFillColor(0)
  pt.SetTextSize(0.0354545)
  text = pt.AddText("CMS Preliminary")
  pt.Draw()
  # lumi
  pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  pt2.SetFillColor(0)
  pt2.SetBorderSize(1)
  pt2.SetLineColor(0)
  pt2.SetTextSize(0.0354545)
  text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  pt2.Draw()

  gPad.RedrawAxis()

  plotname = "limit_k_%.2f.pdf" % modelPointArray[0].coupling
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = "limit_k_%.2f.C" % modelPointArray[0].coupling
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "limit_k_%.2f.png" % modelPointArray[0].coupling
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  # write
  mg.SetName("limit_k_%.2f_MultiGraph" % modelPointArray[0].coupling)
  mg.Write()


def GetMassLimit(modelPointArray):
  #print '----- GetMassLimit -----'
  #print 'Coupling:',modelPointArray[0].coupling
  # fill arrays
  masses = []
  limitExp = []
  limitObs = []
  totalXSec = []
  # test last model point to see if there are limits for all points
  # if not, return -1
  try:
    value = float(modelPointArray[len(modelPointArray)-1].expLimit)
  except ValueError:
    return -1, -1, -1, -1
  for modelPoint in modelPointArray:
    masses.append(modelPoint.mass)
    limitExp.append(modelPoint.expLimit)
    limitObs.append(modelPoint.obsLimit)
    totalXSec.append(modelPoint.totalXSec)
 
  # scan for r = 1
  for i in range(len(modelPointArray)):
    if modelPointArray[i].obsLimit/modelPointArray[i].totalXSec < 1 and \
       modelPointArray[i+1].obsLimit/modelPointArray[i+1].totalXSec > 1:
      rLow = modelPointArray[i].obsLimit/modelPointArray[i].totalXSec
      rHigh = modelPointArray[i+1].obsLimit/modelPointArray[i+1].totalXSec
      xsLow = modelPointArray[i].obsLimit
      xsHigh = modelPointArray[i+1].obsLimit
      mLow = modelPointArray[i].mass
      mHigh = modelPointArray[i+1].mass
      break

  slope =  (rLow - rHigh)/(mLow - mHigh);
  m = (1.0 - rLow)/slope + mLow;
  xsSlope = (xsLow - xsHigh)/(mLow - mHigh);
  xs = xsSlope * (m - mHigh) + xsHigh;

  #print 'obs: rLow:',rLow
  #print 'obs: rHigh:',rHigh
  #print 'obs: mLow:',mLow
  #print 'obs: mHigh:',mHigh

  #   EXPECTED LIMITS SECTION
  #     do it again for expected limits
  for i in range(len(modelPointArray)):
    if modelPointArray[i].expLimit/modelPointArray[i].totalXSec < 1 and \
       modelPointArray[i+1].expLimit/modelPointArray[i+1].totalXSec > 1:
      rLowExp = modelPointArray[i].expLimit/modelPointArray[i].totalXSec
      rHighExp = modelPointArray[i+1].expLimit/modelPointArray[i+1].totalXSec
      xsLowExp = modelPointArray[i].expLimit
      xsHighExp = modelPointArray[i+1].expLimit
      mLowExp = modelPointArray[i].mass
      mHighExp = modelPointArray[i+1].mass
      break

  slopeExp =  (rLowExp - rHighExp)/(mLowExp - mHighExp);
  mExp= (1.0 - rLowExp)/slopeExp + mLowExp;
  xsSlopeExp = (xsLowExp - xsHighExp)/(mLowExp - mHighExp);
  xsExp= xsSlopeExp * (mExp - mHighExp) + xsHighExp;

  #print 'exp: rLow:',rLowExp
  #print 'exp: rHigh:',rHighExp
  #print 'exp: mLow:',mLowExp
  #print 'exp: mHigh:',mHighExp
  #print string.ljust('Coupling: '+str(modelPointArray[0].coupling),14),
  #print ' Observed limit mass: %0.2f'%m
  #print '                Expected limit mass: %0.2f'%mExp
  ##
  #print '                Observed XSec limit: %0.6f'%xs
  #print '                Expected XSec limit: %0.6f'%xsExp
  #   R = 1, M=?
  #     R-r_low/(M-m_low) = (r_low - r_high)/(m_low - m_high)
  return m,xs,mExp,xsExp


def CouplingVsMassPlot(couplingList, expMassLimList, obsMassLimList, rootFile, lumi, imageDir):
  rootFile.cd()
  cLimit = TCanvas("cLimit","cLimit",800,600)
  gStyle.SetOptStat(0)
  #cLimit.Range(595.5973,-0.03483694,1345.283,0.2539571)
  cLimit.SetFillColor(0)
  cLimit.SetBorderMode(0)
  cLimit.SetBorderSize(2)
  cLimit.SetLeftMargin(0.139262)
  cLimit.SetRightMargin(0.0604027)
  cLimit.SetTopMargin(0.0804196)
  cLimit.SetBottomMargin(0.14)
  cLimit.SetFrameBorderMode(0)
  cLimit.SetFrameBorderMode(0)

  # Observed limits
  graph = TGraph(len(couplingList))
  graph.SetName("Graph")
  graph.SetTitle("")
  obsLimFillColor = kRed-4
  graph.SetFillColor(obsLimFillColor)
  graph.SetFillStyle(3001)
  graph.SetLineColor(obsLimFillColor)
  graph.SetLineWidth(10000)
  graph.SetMarkerColor(obsLimFillColor)
  graph.SetMarkerStyle(20)
  graph.SetMarkerSize(1.0)
  for i in range(0,len(couplingList)):
    graph.SetPoint(i,obsMassLimList[i],couplingList[i])

  # Expected limits
  graph2 = TGraph(len(couplingList))
  graph2.SetMarkerColor(kRed)
  graph2.SetMarkerStyle(20)
  graph2.SetMarkerSize(0.0)
  for i in range(0,len(couplingList)):
    graph2.SetPoint(i,expMassLimList[i],couplingList[i])
  graph2.SetLineStyle(7)
  graph2.SetLineColor(kRed+2)
  graph2.SetLineWidth(2)
  graph2.GetXaxis().SetLabelFont(42)
  graph2.GetYaxis().SetLabelFont(42)

  # Previous CMS limits at 7 TeV, 2.2/fb; see arXiv:1112.0688v2
  graphPrevCMS = TGraph(10)
  graphPrevCMS.SetName("graphPrevCMS")
  graphPrevCMS.SetLineStyle(2)
  graphPrevCMS.SetPoint(0,860,0.01)
  graphPrevCMS.SetPoint(1,1130,0.02)
  graphPrevCMS.SetPoint(2,1270,0.03)
  graphPrevCMS.SetPoint(3,1390,0.04)
  graphPrevCMS.SetPoint(4,1500,0.05)
  graphPrevCMS.SetPoint(5,1590,0.06)
  graphPrevCMS.SetPoint(6,1670,0.07)
  graphPrevCMS.SetPoint(7,1740,0.08)
  graphPrevCMS.SetPoint(8,1800,0.09)
  graphPrevCMS.SetPoint(9,1840,0.1)
  #graphPrevCMS.SetLineColor(kRed+2)
  #graphPrevCMS.SetLineColor(kRed+3)
  graphPrevCMS.SetLineColor(kRed+4)

  # ATLAS limits at 7 TeV, 2.12/fb; see Physics Letters B 710 (2012) 538-556
  # See table 3, use result with k-factor=1.75, all channels
  graphAtlas = TGraph(4)
  graphAtlas.SetName("graphAtlas")
  graphAtlas.SetLineStyle(3)
  graphAtlas.SetPoint(0,800,0.01)
  graphAtlas.SetPoint(1,1370,0.03)
  graphAtlas.SetPoint(2,1550,0.05)
  graphAtlas.SetPoint(3,1950,0.1)
  #graphAtlas.SetLineColor(kGreen+2)
  graphAtlas.SetLineColor(kGreen+3)


  # M_D > 10 TeV
  LambdaPi = TF1("LambdaPi","pol1",250,2500)
  LambdaPi.SetParameter(0,0)
  LambdaPi.SetParError(0,0)
  LambdaPi.SetParLimits(0,0,0)
  LambdaPi.SetParameter(1,2.61097e-05)
  LambdaPi.SetParError(1,0)
  LambdaPi.SetParLimits(1,0,0)
  LambdaPiGraph = TGraph(LambdaPi)
  LambdaPiGraph.SetFillColor(kOrange-8)
  LambdaPiGraph.SetLineColor(kBlack)
  LambdaPiGraph.SetFillStyle(1001)
  LambdaPiGraph.SetLineWidth(-10000)

  # Electroweak limits
  graphEW = TGraph(27)
  graphEW.SetName("graphEW")
  graphEW.SetTitle("graphEW")
  graphEW.SetLineStyle(5)
  graphEW.SetPoint(0,180,0.1071)
  graphEW.SetPoint(1,183,0.1062)
  graphEW.SetPoint(2,190,0.1043)
  graphEW.SetPoint(3,200,0.1016)
  graphEW.SetPoint(4,210,0.0989)
  graphEW.SetPoint(5,220,0.0963)
  graphEW.SetPoint(6,230,0.0938)
  graphEW.SetPoint(7,240,0.0913)
  graphEW.SetPoint(8,250,0.0889)
  graphEW.SetPoint(9,260,0.0866)
  graphEW.SetPoint(10,270,0.0843)
  graphEW.SetPoint(11,280,0.0821)
  graphEW.SetPoint(12,290,0.0799)
  graphEW.SetPoint(13,300,0.0778)
  graphEW.SetPoint(14,400,0.0603)
  graphEW.SetPoint(15,500,0.0481)
  graphEW.SetPoint(16,600,0.0397)
  graphEW.SetPoint(17,700,0.0337)
  graphEW.SetPoint(18,800,0.0292)
  graphEW.SetPoint(19,900,0.0258)
  graphEW.SetPoint(20,1000,0.0231)
  graphEW.SetPoint(21,1100,0.0209)
  graphEW.SetPoint(22,1200,0.0191)
  graphEW.SetPoint(23,1300,0.0176)
  graphEW.SetPoint(24,1400,0.0163)
  graphEW.SetPoint(25,1500,0.0152)
  graphEW.SetPoint(26,2200,0.01)
  graphEW.GetXaxis().SetLabelFont(42)
  graphEW.GetYaxis().SetLabelFont(42)
  graphEW.SetLineWidth(-10000)
  graphEW.SetLineColor(kBlack)
  graphEW.SetFillColor(kAzure-4)
  graphEW.SetFillStyle(1001)

  # Draw limExp, limObs, EW limits -- same axes
  mg = TMultiGraph()
  mg.Add(LambdaPiGraph,'C')
  mg.Add(graphEW,'C')
  mg.Add(graph,'L')
  mg.Add(graph2,'L')
  mg.Add(graphPrevCMS,'L')
  mg.Add(graphAtlas,'L')
  xArray = graph.GetX()
  xArray.SetSize(graph.GetN())
  xList = list(xArray)
  cLimit.Clear()
  mg.Draw('a')
  #mg.GetXaxis().SetRangeUser(xList[0],xList[len(xList)-1])
  mg.GetXaxis().SetRangeUser(xList[0],xList[len(xList)-1]-200)
  mg.SetMinimum(0)
  mg.SetMaximum(0.12)
  mg.GetXaxis().SetTitle("M_{1} [GeV]")
  mg.GetYaxis().SetTitle("Coupling k/#bar{M}_{Pl}")
  mg.GetXaxis().SetLabelFont(42)
  mg.GetYaxis().SetLabelFont(42)
  mg.GetXaxis().SetLabelSize(0.04)
  mg.GetYaxis().SetLabelSize(0.04)
  mg.GetYaxis().SetTitleOffset(1.19)
  mg.GetXaxis().SetTitleOffset(1.0)
  mg.GetXaxis().SetTickLength()
  mg.GetYaxis().SetTickLength()
  mg.GetXaxis().Draw()
  mg.GetYaxis().Draw()
  gPad.RedrawAxis()
  gPad.Update()

  # Legend
  #leg = TLegend(0.207,0.467,0.401,0.752,"","brNDC")
  leg = TLegend(0.16,0.59,0.43,0.89,"","brNDC")
  leg.SetBorderSize(1)
  leg.SetTextFont(62)
  leg.SetTextSize(0.0275)
  leg.SetLineColor(0)
  leg.SetLineStyle(0)
  leg.SetLineWidth(0)
  leg.SetFillColor(kWhite)
  leg.SetFillStyle(1001)
  entry=leg.AddEntry(graphEW,"Electroweak Limits","lf")
  entry.SetLineColor(1)
  entry.SetLineStyle(5)
  entry.SetLineWidth(3)
  #entry.SetMarkerColor(1)
  #entry.SetMarkerStyle(21)
  #entry.SetMarkerSize(1)
  #entry.SetMarkerStyle(20)
  #entry.SetMarkerSize(1.3)
  entry=leg.AddEntry(graph,"CMS 95% CL Limit","f")
  leg.AddEntry(graph2,"CMS Expected Limit","l")
  #entry.SetMarkerStyle(21)
  #entry.SetMarkerSize(1.3)
  leg.AddEntry(graphPrevCMS,"CMS 2.2 fb^{-1} at 7 TeV","l")
  leg.AddEntry(graphAtlas,"ATLAS 2.2 fb^{-1} at 7 TeV","l")
  entry=leg.AddEntry(LambdaPiGraph,"M_{D} > 10TeV","lf")
  leg.SetHeader("CMS %.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  leg.Draw()
  ## lumi
  #pt = TPaveText(0.629397,0.798951,0.8002513,0.8653846,"blNDC")
  #pt.SetFillColor(0)
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetTextSize(0.06)
  #text = pt.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt.Draw()
  ## CMS
  #ptCMS = TPaveText(0.2236181,0.7884615,0.4736181,0.8583916,"blNDC")
  #ptCMS.SetName("CMS")
  #ptCMS.SetBorderSize(1)
  #ptCMS.SetLineColor(0)
  #ptCMS.SetFillColor(0)
  #ptCMS.SetTextSize(0.06)
  #text = ptCMS.AddText("CMS")
  #ptCMS.Draw()

  cLimit.SetLogy(0)
  cLimit.Modified()
  cLimit.cd()
  cLimit.SetSelected(cLimit)
  cLimit.Update()
  cLimit.Print(imageDir+'/RSMassVsCouplingLimits.C')
  cLimit.Print(imageDir+'/RSMassVsCouplingLimits.pdf')
  cLimit.Print(imageDir+'/RSMassVsCouplingLimits.eps')
  #cLimit.Print(imageDir+'/RSMassVsCouplingLimits.png')
  baseName=imageDir+'/RSMassVsCouplingLimits'
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+baseName+'.png',baseName+'.eps'])
  cLimit.SetName("RSMassVsCouplingLimits")
  cLimit.Write()


def PlotAllEfficienciesMScale(modelPointArrays, lumi, rootFile, imageDir):
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  effArrs = []
  effMScaleUpArrs = []
  effMScaleDownArrs = []
  effMScaleUpNonAbsArrs = []
  effMScaleDownNonAbsArrs = []
  for mpArray in modelPointArrays:
    if len(mpArray) < 1:
      continue
    masses, effs = makeEffArrays(mpArray)
    massesArrs.append(masses)
    effArrs.append(effs)
    couplings.append(mpArray[0].coupling)
    effMScaleUp = []
    effMScaleDown = []
    effMScaleUpNonAbs = []
    effMScaleDownNonAbs = []
    for mp in mpArray:
      effMScaleUp.append(math.fabs(mp.totalEffMScaleSystUp-mp.totalEff)/mp.totalEff)
      effMScaleDown.append(math.fabs(mp.totalEffMScaleSystDown-mp.totalEff)/mp.totalEff)
      effMScaleUpNonAbs.append((mp.totalEffMScaleSystUp-mp.totalEff)/mp.totalEff)
      effMScaleDownNonAbs.append((mp.totalEffMScaleSystDown-mp.totalEff)/mp.totalEff)
    effMScaleUpArrs.append(effMScaleUp)
    effMScaleDownArrs.append(effMScaleDown)
    effMScaleUpNonAbsArrs.append(effMScaleUpNonAbs)
    effMScaleDownNonAbsArrs.append(effMScaleDownNonAbs)

  rootFile.cd()
  # turn the arrays into graphs
  effGraphs = []
  effGraphsNonAbsUp = []
  effGraphsNonAbsDown = []
  for massesArr, effsArr, effMScaleDownArr, effMScaleUpArr, effMScaleUpNonAbsArr, effMScaleDownNonAbsArr in itertools.izip(massesArrs,effArrs,effMScaleDownArrs,effMScaleUpArrs,effMScaleUpNonAbsArrs,effMScaleDownNonAbsArrs):
    massErrs = [0]*len(massesArr)
    nominals = [1]*len(massesArr)
    effGraphs.append(TGraphAsymmErrors(len(massesArr), array.array("f",massesArr),array.array("f",nominals),array.array("f",massErrs),array.array("f",massErrs),array.array("f",effMScaleDownArr),array.array("f",effMScaleUpArr)))
    effGraphsNonAbsUp.append(TGraph(len(massesArr), array.array("f",massesArr), array.array("f",effMScaleUpNonAbsArr)))
    effGraphsNonAbsDown.append(TGraph(len(massesArr), array.array("f",massesArr), array.array("f",effMScaleDownNonAbsArr)))
    #print 'effMScaleDownArr:',effMScaleDownArr
    #print 'effMScaleDownArr:',effMScaleDownArr

  SetCustomGStyle()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  c.SetGridy()
  # Multigraph
  mg = TMultiGraph()
  mgNonAbs = TMultiGraph()
  legend = TLegend(0.2,0.75,0.5,0.9)
  legendNonAbs = TLegend(0.2,0.75,0.5,0.9)
  colorIndex = 2 
  drawOpt = 'l3'
  legDrawOpt = 'f'
  drawOptNonAbs = 'lp'
  legDrawOptNonAbs = 'lp'
  for effGraph, effGraphNAUp, effGraphNADown, coupling in itertools.izip(effGraphs,effGraphsNonAbsUp,effGraphsNonAbsDown,couplings):
    effGraph.SetMarkerSize(0.8)
    effGraph.SetMarkerColor(colorIndex)
    effGraph.SetLineColor(colorIndex)
    effGraph.SetLineWidth(4)
    effGraph.SetLineStyle(2)
    effGraph.SetFillColor(colorIndex)
    effGraph.SetFillStyle(3003)
    mg.Add(effGraph,drawOpt)
    legend.AddEntry(effGraph," #tilde{k} = "+str(coupling),legDrawOpt)
    effGraphNAUp.SetMarkerSize(0.8)
    effGraphNAUp.SetMarkerColor(colorIndex)
    effGraphNAUp.SetLineColor(colorIndex)
    effGraphNAUp.SetLineWidth(4)
    effGraphNAUp.SetName('mscaleShiftUpK='+str(coupling).replace('.','p'))
    effGraphNAUp.Write()
    effGraphNADown.SetMarkerSize(0.8)
    effGraphNADown.SetMarkerColor(colorIndex)
    effGraphNADown.SetLineColor(colorIndex)
    effGraphNADown.SetLineWidth(4)
    effGraphNADown.SetLineStyle(2)
    effGraphNADown.SetName('mscaleShiftDownK'+str(coupling).replace('.','p'))
    effGraphNADown.Write()
    mgNonAbs.Add(effGraphNAUp,drawOptNonAbs)
    mgNonAbs.Add(effGraphNADown,drawOptNonAbs)
    legendNonAbs.AddEntry(effGraphNAUp," scaleShiftUp: #tilde{k} = "+str(coupling),legDrawOptNonAbs)
    legendNonAbs.AddEntry(effGraphNADown," scaleShiftDown: #tilde{k} = "+str(coupling),legDrawOptNonAbs)
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1
  # To make the separate graphs share the same axes
  h = TH1F("test","",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0.8,1.2)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("rel. efficiency*acc change")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw()
  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()
  gPad.RedrawAxis()
  basename = 'allEfficienciesMScaleSyst'
  plotname = basename+'.pdf'
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.C'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.eps'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = basename+'.png'
  savename = TString(plotname)
  fullBasename = imageDir+'/'+basename
  #subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+fullBasename+'.png',fullBasename+'.eps'])
  # write
  name = TString(basename)
  mg.SetName(name.Data())
  mg.Write()
  #
  # for non-abs
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  c.SetRightMargin(0.04)
  c.SetGridy()
  # To make the separate graphs share the same axes
  h = TH1F("test","",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(-0.2,0.2)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("rel. efficiency*acc change")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mgNonAbs.Draw()
  # draw legend
  legendNonAbs.SetBorderSize(0)
  legendNonAbs.SetFillColor(0)
  legendNonAbs.Draw()
  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()
  gPad.RedrawAxis()
  basename = 'allEfficienciesMScaleSystNonAbs'
  plotname = basename+'.pdf'
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.C'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.eps'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = basename+'.png'
  savename = TString(plotname)
  fullBasename = imageDir+'/'+basename
  #subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+fullBasename+'.png',fullBasename+'.eps'])
  # write
  name = TString(basename)
  mgNonAbs.SetName(name.Data())
  mgNonAbs.Write()


def PlotAllEfficienciesMRes(modelPointArrays, lumi, rootFile, imageDir):
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  effArrs = []
  effMResUpArrs = []
  effMResDownArrs = []
  for mpArray in modelPointArrays:
    if len(mpArray) < 1:
      continue
    masses, effs = makeEffArrays(mpArray)
    massesArrs.append(masses)
    effArrs.append(effs)
    couplings.append(mpArray[0].coupling)
    effMResUp = []
    effMResDown = []
    for mp in mpArray:
      #print 'mp.totalEffMResSystUp=',mp.totalEffMResSystUp,'mp.totalEff=',mp.totalEff,'relDiff=',(mp.totalEffMResSystUp-mp.totalEff)/mp.totalEff
      #print 'mp.totalEffMResSystDown=',mp.totalEffMResSystDown,'mp.totalEff=',mp.totalEff,'relDiff=',(mp.totalEffMResSystDown-mp.totalEff)/mp.totalEff
      #effMResUp.append(math.fabs(mp.totalEffMResSystUp-mp.totalEff)/mp.totalEff)
      effMResUp.append((mp.totalEffMResSystUp-mp.totalEff)/mp.totalEff)
      effMResDown.append(math.fabs(mp.totalEffMResSystDown-mp.totalEff)/mp.totalEff) # not used
    effMResUpArrs.append(effMResUp)
    effMResDownArrs.append(effMResDown)

  rootFile.cd()
  # turn the arrays into graphs
  effGraphs = []
  for massesArr, effsArr, effMResDownArr, effMResUpArr in itertools.izip(massesArrs,effArrs,effMResDownArrs,effMResUpArrs):
    massErrs = [0]*len(massesArr)
    nominals = [1]*len(massesArr)
    #effGraphs.append(TGraphAsymmErrors(len(massesArr), array.array("f",massesArr),array.array("f",nominals),array.array("f",massErrs),array.array("f",massErrs),array.array("f",effMResDownArr),array.array("f",effMResUpArr)))
    effGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",effMResUpArr)))

  SetCustomGStyle()
  gStyle.SetFuncColor(1)
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  c.SetGridy()
  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.2,0.75,0.5,0.9)
  colorIndex = 2 
  #drawOpt = 'l3'
  #legDrawOpt = 'f'
  drawOpt = 'lp'
  legDrawOpt = 'lp'
  for effGraph, coupling in itertools.izip(effGraphs,couplings):
    effGraph.SetMarkerSize(0.8)
    effGraph.SetMarkerColor(colorIndex)
    effGraph.SetLineColor(colorIndex)
    effGraph.SetLineWidth(4)
    #effGraph.SetLineStyle(2)
    #effGraph.SetFillColor(colorIndex)
    #effGraph.SetFillStyle(3003)
    effGraph.SetName('mResShiftK'+str(coupling).replace('.','p'))
    effGraph.Fit('pol1')
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(effGraph,drawOpt)
    legend.AddEntry(effGraph," #tilde{k} = "+str(coupling),legDrawOpt)
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1
  # To make the separate graphs share the same axes
  h = TH1F("test","",10,750,3500)
  h.SetStats(False)
  #h.GetYaxis().SetRangeUser(0.8,1.2)
  h.GetYaxis().SetRangeUser(-0.25,0.25)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("rel. efficiency*acc change")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw()
  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()
  gPad.RedrawAxis()
  basename = 'allEfficienciesMResSyst'
  plotname = basename+'.pdf'
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.C'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.eps'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = basename+'.png'
  savename = TString(plotname)
  fullBasename = imageDir+'/'+basename
  #subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+fullBasename+'.png',fullBasename+'.eps'])
  # write
  name = TString(basename)
  mg.SetName(name.Data())
  mg.Write()


def PlotAllEfficienciesPileup(modelPointArrays, lumi, rootFile, imageDir):
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  effArrs = []
  effPileupUpArrs = []
  effPileupDownArrs = []
  for mpArray in modelPointArrays:
    if len(mpArray) < 1:
      continue
    masses, effs = makeEffArrays(mpArray)
    massesArrs.append(masses)
    effArrs.append(effs)
    couplings.append(mpArray[0].coupling)
    effPileupUp = []
    effPileupDown = []
    for mp in mpArray:
      effPileupUp.append(math.fabs(mp.totalEffPileupSystUp-mp.totalEff)/mp.totalEff)
      effPileupDown.append(math.fabs(mp.totalEffPileupSystDown-mp.totalEff)/mp.totalEff)
    effPileupUpArrs.append(effPileupUp)
    effPileupDownArrs.append(effPileupDown)

  rootFile.cd()
  # turn the arrays into graphs
  effGraphs = []
  for massesArr, effsArr, effPileupDownArr, effPileupUpArr in itertools.izip(massesArrs,effArrs,effPileupDownArrs,effPileupUpArrs):
    massErrs = [0]*len(massesArr)
    nominals = [1]*len(massesArr)
    effGraphs.append(TGraphAsymmErrors(len(massesArr), array.array("f",massesArr),array.array("f",nominals),array.array("f",massErrs),array.array("f",massErrs),array.array("f",effPileupDownArr),array.array("f",effPileupUpArr)))

  SetCustomGStyle()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  c.SetGridy()
  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.2,0.75,0.5,0.9)
  colorIndex = 2 
  drawOpt = 'l3'
  legDrawOpt = 'f'
  for effGraph, coupling in itertools.izip(effGraphs,couplings):
    effGraph.SetMarkerSize(0.8)
    effGraph.SetMarkerColor(colorIndex)
    effGraph.SetLineColor(colorIndex)
    effGraph.SetLineWidth(4)
    effGraph.SetLineStyle(2)
    effGraph.SetFillColor(colorIndex)
    effGraph.SetFillStyle(3003)
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(effGraph,drawOpt)
    legend.AddEntry(effGraph," #tilde{k} = "+str(coupling),legDrawOpt)
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1
  # To make the separate graphs share the same axes
  h = TH1F("test","",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0.99,1.01)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("rel. efficiency*acc change")
  h.GetYaxis().SetNdivisions(515)
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.7)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw()
  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()
  gPad.RedrawAxis()
  basename = 'allEfficienciesPileupSyst'
  plotname = basename+'.pdf'
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.C'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.eps'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = basename+'.png'
  savename = TString(plotname)
  fullBasename = imageDir+'/'+basename
  #subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+fullBasename+'.png',fullBasename+'.eps'])
  # write
  name = TString(basename)
  mg.SetName(name.Data())
  mg.Write()


def PlotAllEfficiencies(modelPointArrays, lumi, rootFile, imageDir):
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  effArrs = []
  for mpArray in modelPointArrays:
    if len(mpArray) < 1:
      continue
    masses, effs = makeEffArrays(mpArray)
    massesArrs.append(masses)
    effArrs.append(effs)
    couplings.append(mpArray[0].coupling)

  rootFile.cd()
  # turn the arrays into graphs
  effGraphs = []
  for massesArr, effsArr in itertools.izip(massesArrs,effArrs):
    effGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",effsArr)))

  SetCustomGStyle()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.42,0.71,0.73,0.92)
  colorIndex = 2 
  drawOpt = 'lx'
  for effGraph, coupling in itertools.izip(effGraphs,couplings):
    effGraph.SetMarkerSize(0.8)
    effGraph.SetMarkerColor(colorIndex)
    effGraph.SetLineColor(colorIndex)
    effGraph.SetLineWidth(4)
    effGraph.SetLineStyle(2)
    effGraph.SetFillColor(colorIndex)
    effGraph.SetFillStyle(3003)
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(effGraph,drawOpt)
    legend.AddEntry(effGraph," #tilde{k} = "+str(coupling),drawOpt)
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1
  # To make the separate graphs share the same axes
  h = TH1F("test","",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0.2,0.7)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("efficiency*acc")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw("L")
  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()
  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()
  gPad.RedrawAxis()
  basename = 'allEfficiencies'
  plotname = basename+'.pdf'
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = basename+'.C'
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = basename+'.png'
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  # write
  name = TString(basename)
  mg.SetName(name.Data())
  mg.Write()


def PlotAllHalfWidths(modelPointArrays, lumi, rootFile, imageDir):
  print 'PlotAllHalfWidths'
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  halfWidthArrs = []
  for mpArray in modelPointArrays:
    masses, halfWidths = makeWidthArrays(mpArray)
    massesArrs.append(masses)
    couplings.append(mpArray[0].coupling)
    halfWidthArrs.append(halfWidths)

  rootFile.cd()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  # turn the arrays into graphs
  halfWidthGraphs = []
  for massesArr, halfWidthsArr in itertools.izip(massesArrs,halfWidthArrs):
    halfWidthGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",halfWidthsArr)))

  SetCustomGStyle()

  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.42,0.71,0.73,0.92)
  colorIndex = 2 
  for halfWidthGraph, coupling in itertools.izip(halfWidthGraphs,couplings):
    halfWidthGraph.SetMarkerSize(0.8)
    halfWidthGraph.SetMarkerColor(colorIndex)
    halfWidthGraph.SetLineColor(colorIndex)
    halfWidthGraph.SetLineWidth(4)
    halfWidthGraph.SetLineStyle(2)
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(halfWidthGraph,"l")
    legend.AddEntry(halfWidthGraph," halfWidth #tilde{k} = "+str(coupling),"l")
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1


  # To make the separate graphs share the same axes
  h = TH1F("test","test",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0,50)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("halfWidth [GeV]")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw("L")

  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()

  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  #pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  #pt2.SetFillColor(0)
  #pt2.SetBorderSize(1)
  #pt2.SetLineColor(0)
  #pt2.SetTextSize(0.0354545)
  #text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  #pt2.Draw()

  gPad.RedrawAxis()

  plotname = "allHalfWidths.pdf"
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = "allHalfWidths.C"
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allHalfWidths.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  # write
  mg.SetName("allHalfWidths")
  mg.Write()


def PlotAllExpBGs(modelPointArrays, lumi, rootFile, imageDir):
  print 'PlotAllExpBGs'
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  bgArrs = []
  for mpArray in modelPointArrays:
    masses, bgs = makeExpBGArrays(mpArray)
    massesArrs.append(masses)
    couplings.append(mpArray[0].coupling)
    bgArrs.append(bgs)

  rootFile.cd()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  c.SetLogy()
  c.SetRightMargin(0.04)
  # turn the arrays into graphs
  bgGraphs = []
  for massesArr, bgsArr in itertools.izip(massesArrs,bgArrs):
    bgGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",bgsArr)))

  SetCustomGStyle()

  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.42,0.71,0.73,0.92)
  colorIndex = 2 
  for bgGraph, coupling in itertools.izip(bgGraphs,couplings):
    bgGraph.SetMarkerSize(0.8)
    bgGraph.SetMarkerColor(colorIndex)
    bgGraph.SetLineColor(colorIndex)
    bgGraph.SetLineWidth(4)
    bgGraph.SetLineStyle(2)
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(bgGraph,"l")
    legend.AddEntry(bgGraph," Exp. BG #tilde{k} = "+str(coupling),"l")
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1


  # To make the separate graphs share the same axes
  h = TH1F("test","test",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0.001,4.0)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("Exp. BG")
  h.GetXaxis().SetLabelFont(42)
  h.GetYaxis().SetLabelFont(42)
  h.GetXaxis().SetLabelSize(0.04)
  h.GetYaxis().SetLabelSize(0.04)
  h.GetYaxis().SetTitleOffset(1.5)
  h.GetXaxis().SetTitleOffset(1.2)
  h.GetXaxis().SetTitleSize(0.04)
  h.GetYaxis().SetTitleSize(0.04)
  h.Draw()
  mg.Draw("L")

  # draw legend
  legend.SetBorderSize(0)
  legend.SetFillColor(0)
  legend.Draw()

  ## CMS
  #pt = TPaveText(0.645973,0.629371,0.845638,0.699301,"blNDC")
  #pt.SetName("CMS Preliminary")
  #pt.SetBorderSize(1)
  #pt.SetLineColor(0)
  #pt.SetFillColor(0)
  #pt.SetTextSize(0.0354545)
  #text = pt.AddText("CMS Preliminary")
  #pt.Draw()
  ## lumi
  pt2 = TPaveText(0.654362,0.585664,0.825503,0.652098,"blNDC")
  pt2.SetFillColor(0)
  pt2.SetBorderSize(1)
  pt2.SetLineColor(0)
  pt2.SetTextSize(0.0354545)
  text = pt2.AddText("%.1f" % (lumi/1000)+" fb^{-1} at 8 TeV")
  pt2.Draw()

  gPad.RedrawAxis()

  plotname = "allExpBGs.pdf"
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(imageDir+'/'+savename.Data())
  plotname = "allExpBGs.C"
  savename = TString(plotname)
  c.SaveAs(imageDir+'/'+savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allExpBGs.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',imageDir+'/'+pdfName,imageDir+'/'+savename.Data()])
  # write
  mg.SetName("allExpBGs")
  mg.Write()

