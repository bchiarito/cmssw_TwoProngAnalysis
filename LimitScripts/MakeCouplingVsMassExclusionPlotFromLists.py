#!/usr/bin/env python

import string
import os
import sys
import subprocess

from ROOT import *


# Macro to make RS Model exclusion plot in 2-D plane from lists of exp/obs limits + couplings
# Taken from MakePlots.py code (added SetNdivisions on y-axis to get proper tick labeling)
#
# May 20, 2014
# Seth I. Cooper (Alabama)
#
# Modify lists of limits/couplings and run below


def CouplingVsMassPlot(couplingList, expMassLimList, obsMassLimList, lumi):
  cLimit = TCanvas("cLimit","cLimit",1000,750)
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
  #for i in range(0,len(couplingList)):
  #  graph.SetPoint(i,obsMassLimList[i],couplingList[i])
  # slope+x-intercept at y=1
  slope = (couplingList[1]-couplingList[0])/(obsMassLimList[1]-obsMassLimList[0])
  xInt = (couplingList[0]-1)/slope + obsMassLimList[0]
  graph.SetPoint(0,xInt,1)
  for i in range(0,len(couplingList)):
    graph.SetPoint(i+1,obsMassLimList[i],couplingList[i])
  graph.SetPoint(len(couplingList)+1,obsMassLimList[len(couplingList)-1],1)

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
  #mg.Add(graph,'L')
  mg.Add(graph,'FL')
  mg.Add(graph2,'L')
  mg.Add(graphPrevCMS,'L')
  mg.Add(graphAtlas,'L')
  cLimit.Clear()
  gPad.DrawFrame(obsMassLimList[0]-20,0,obsMassLimList[len(obsMassLimList)-1],couplingList[len(couplingList)-1]+0.02)
  mg.Draw()
  mg.GetXaxis().SetTitle("M_{1} [GeV]")
  mg.GetYaxis().SetTitle("Coupling k/#bar{M}_{Pl}")
  mg.GetXaxis().SetLabelFont(42)
  mg.GetYaxis().SetLabelFont(42)
  mg.GetXaxis().SetLabelSize(0.04)
  mg.GetYaxis().SetLabelSize(0.04)
  mg.GetXaxis().SetLabelOffset(0.01)
  mg.GetYaxis().SetTitleOffset(1.18)
  mg.GetXaxis().SetTitleOffset(1.0)
  mg.GetXaxis().SetTickLength()
  mg.GetYaxis().SetTickLength()
  #
  mg.GetXaxis().SetNdivisions(510)
  #
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
  entry=leg.AddEntry(LambdaPiGraph,"#Lambda_{#pi} < 10TeV","lf")
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
  cLimit.Print('RSMassVsCouplingLimits.C')
  cLimit.Print('RSMassVsCouplingLimits.pdf')
  cLimit.Print('RSMassVsCouplingLimits.eps')
  #cLimit.Print('RSMassVsCouplingLimits.png')
  baseName='RSMassVsCouplingLimits'
  subprocess.call(['gs','-dTextAlphaBits=4','-dBATCH','-dNOPAUSE','-dQUIET','-dEPSCrop','-sDEVICE=png16m','-sOutputFile='+baseName+'.png',baseName+'.eps'])




###########
### RUN ###
###########

couplings = [ 0.01, 0.05, 0.1]
massLimExp = [1277,2153,2538]
massLimObs = [1233,2161,2544]
lumi = 19620

CouplingVsMassPlot(couplings,massLimExp,massLimObs,lumi)

