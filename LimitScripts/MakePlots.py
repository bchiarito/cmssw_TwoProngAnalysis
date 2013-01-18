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


def ReadFromFile(file):
  outputModelPoints = []
  for line in file:
    #print line,
    if "Coupling" in line:
      mp = ModelPoint()
      mp.coupling = float(line.split(': ')[1])
    elif "Mass" in line:
      mp.mass = float(line.split(': ')[1])
    elif "TotalXSection" in line:
      mp.totalXSec = float(line.split(': ')[1])
    elif "TotalEff" in line:
      mp.totalEff = float(line.split(': ')[1])
    elif "HalfWidth" in line:
      mp.halfWidth = float(line.split(': ')[1])
    elif "NDataObs" in line:
      mp.nDataObs = float(line.split(': ')[1])
    elif "NBackground" in line and not "Err" in line:
      mp.nBackground = float(line.split(': ')[1])
    elif "NBackgroundErr" in line:
      mp.nBackgroundErr = float(line.split(': ')[1])
    elif "ExpectedLimit" in line and not "Sigma" in line:
      mp.expLimit = float(line.split(': ')[1])
    elif "ExpectedLimitOneSigmaHigh" in line:
      mp.expLimitOneSigmaHigh = float(line.split(': ')[1])
    elif "ExpectedLimitOneSigmaLow" in line:
      mp.expLimitOneSigmaLow = float(line.split(': ')[1])
    elif "ExpectedLimitTwoSigmaHigh" in line:
      mp.expLimitTwoSigmaHigh = float(line.split(': ')[1])
    elif "ExpectedLimitTwoSigmaLow" in line:
      mp.expLimitTwoSigmaLow = float(line.split(': ')[1])
    elif "ObservedLimit" in line:
      mp.obsLimit = float(line.split(': ')[1])
      outputModelPoints.append(mp)
  return outputModelPoints


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


def PlotAllBands(modelPointArrays, lumi, rootFile):
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
    massesArrs.append(masses)
    limitExpArrs.append(limitExp)
    limitObsArrs.append(limitObs)
    totalXSecArrs.append(totalXSec)
    couplings.append(mpArray[0].coupling)

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
  legend = TLegend(0.42,0.71,0.73,0.92)
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
  h.GetYaxis().SetRangeUser(4e-4,3e-3)
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
  c.SaveAs(savename.Data())
  plotname = "allLimits.C"
  savename = TString(plotname)
  c.SaveAs(savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allLimits.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',pdfName,savename.Data()])
  # write
  mg.SetName("limitsAllCouplingsMultiGraph")
  mg.Write()


def PlotBands(modelPointArray, lumi, rootFile):
  #print '----- PlotBands -----'
  #print 'Coupling:',modelPointArray[0].coupling
  # fill arrays
  errMass, errUp, errDn, err2Up, err2Dn, masses, expUp, expDn, exp2Up, exp2Dn, limitExp, limitObs, totalXSec = makeArrays(modelPointArray)

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


  #  g01_obs.SetMarkerStyle(20)
  g01_obs.SetMarkerSize(0.4)
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
  mg.Add(g01_obs,"L")
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
  c.SaveAs(savename.Data())
  plotname = "limit_k_%.2f.C" % modelPointArray[0].coupling
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  c.SaveAs(savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "limit_k_%.2f.png" % modelPointArray[0].coupling
  savename = TString(plotname)
  indexstring = savename.Index(".")
  savename.Replace(indexstring,1,"p")
  subprocess.Popen(['convert','-trim',pdfName,savename.Data()])
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


def CouplingVsMassPlot(couplingList, expMassLimList, obsMassLimList, rootFile):
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
  #graph = TGraph(3)
  graph = TGraph(len(couplingList))
  graph.SetName("Graph")
  graph.SetTitle("")
  graph.SetFillColor(1)
  ci = TColor.GetColor("#ff0000")
  graph.SetLineColor(ci)
  graph.SetLineWidth(3)
  graph.SetMarkerColor(ci)
  graph.SetMarkerStyle(20)
  graph.SetMarkerSize(1.0)
  #graph.SetPoint(0,1164.15,0.01)
  #graph.SetPoint(1,1885.73,0.05)
  #graph.SetPoint(2,2293.09,0.1)
  for i in range(0,len(couplingList)):
    graph.SetPoint(i,obsMassLimList[i],couplingList[i])

  # Expected limits
  #graph2 = TGraph(3)
  graph2 = TGraph(len(couplingList))
  graph2.SetMarkerColor(ci)
  graph2.SetMarkerStyle(20)
  graph2.SetMarkerSize(0.0)
  #graph2.SetPoint(0,  1157.02,0.01)
  #graph2.SetPoint(1, 1922.33,0.05)
  #graph2.SetPoint(2, 2299.17,0.1)
  for i in range(0,len(couplingList)):
    graph2.SetPoint(i,expMassLimList[i],couplingList[i])
  graph2.SetLineStyle(kDashed)
  graph2.SetLineColor(ci)
  graph2.SetLineWidth(3)
  graph2.GetXaxis().SetLabelFont(42)
  graph2.GetYaxis().SetLabelFont(42)

  # M_D > 10 TeV
  LambdaPi = TF1("LambdaPi","pol1",250,2500)
  LambdaPi.SetFillColor(15)
  LambdaPi.SetFillStyle(3004)
  LambdaPi.SetLineColor(15)
  LambdaPi.SetLineWidth(1)
  LambdaPi.SetParameter(0,0)
  LambdaPi.SetParError(0,0)
  LambdaPi.SetParLimits(0,0,0)
  LambdaPi.SetParameter(1,2.61097e-05)
  LambdaPi.SetParError(1,0)
  LambdaPi.SetParLimits(1,0,0)
  LambdaPi.GetXaxis().SetLabelFont(42)
  LambdaPi.GetYaxis().SetLabelFont(42)

  # Electroweak limits
  graphEW = TGraph(27)
  graphEW.SetName("graphEW")
  graphEW.SetTitle("graphEW")
  graphEW.SetFillColor(1)
  graphEW.SetFillStyle(3004)    
  graphEW.SetLineStyle(5)
  graphEW.SetLineWidth(3)
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
  graphEW.SetFillColor(kBlack)
  graphEW.SetFillStyle(3004)
  graphEW.SetLineWidth(1)
  graphEW.SetFillStyle(3002)
  graphEW.SetFillColor(kBlack)
  graphEW.SetLineWidth(-10000)

  # Draw limExp, limObs, EW limits -- same axes
  mg = TMultiGraph()
  mg.Add(graph)
  mg.Add(graph2)
  mg.Add(graphEW)
  xArray = graph.GetX()
  xArray.SetSize(graph.GetN())
  xList = list(xArray)
  mg.Draw("al")
  mg.GetXaxis().SetRangeUser(xList[0],xList[len(xList)-1])
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
  # Draw LambdaPi
  LambdaPi.Draw("same")

  # Legend
  leg = TLegend(0.2072864,0.4667832,0.4007538,0.7517483,"","brNDC")
  leg.SetBorderSize(1)
  leg.SetTextFont(62)
  leg.SetTextSize(0.05)
  leg.SetLineColor(0)
  leg.SetLineStyle(0)
  leg.SetLineWidth(0)
  leg.SetFillColor(0)
  leg.SetFillStyle(0)
  entry=leg.AddEntry("graph","Electroweak Limits","lf")
  entry.SetLineColor(1)
  entry.SetLineStyle(5)
  entry.SetLineWidth(3)
  entry.SetMarkerColor(1)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1)
  ci = TColor.GetColor("#ff0000")
  entry.SetMarkerColor(ci)
  entry.SetMarkerStyle(20)
  entry.SetMarkerSize(1.3)
  entry=leg.AddEntry("Graph","95% CL Limit","l")
  leg.AddEntry(graph2,"Expected Limit","l")
  ci = TColor.GetColor("#00ff00")
  entry.SetMarkerColor(ci)
  entry.SetMarkerStyle(21)
  entry.SetMarkerSize(1.3)
  entry=leg.AddEntry("LambdaPi","M_{D} > 10TeV","lf")
  leg.Draw()
  # 10.3/fb
  pt = TPaveText(0.629397,0.798951,0.8002513,0.8653846,"blNDC")
  pt.SetFillColor(0)
  pt.SetBorderSize(1)
  pt.SetLineColor(0)
  pt.SetTextSize(0.06)
  text = pt.AddText("10.3 fb^{-1} at 8 TeV")
  pt.Draw()
  # CMS
  ptCMS = TPaveText(0.2236181,0.7884615,0.4736181,0.8583916,"blNDC")
  ptCMS.SetName("CMS")
  ptCMS.SetBorderSize(1)
  ptCMS.SetLineColor(0)
  ptCMS.SetFillColor(0)
  ptCMS.SetTextSize(0.06)
  text = ptCMS.AddText("CMS")
  ptCMS.Draw()

  cLimit.SetLogy(0)
  cLimit.Modified()
  cLimit.cd()
  cLimit.SetSelected(cLimit)
  cLimit.Update()
  cLimit.Print("RSMassVsCouplingLimits.C")
  cLimit.Print("RSMassVsCouplingLimits.png")
  cLimit.Print("RSMassVsCouplingLimits.pdf")
  cLimit.SetName("RSMassVsCouplingLimits")
  cLimit.Write()


def PlotAllEfficiencies(modelPointArrays, lumi, rootFile):
  # this function takes a list of modelPointArrays: [mpsCoupling0.01, mpsCoupling0.05, ...]
  # fill arrays that we care about
  # TODO reuse this code in the plotbands section
  massesArrs = []
  couplings = []
  effArrs = []
  for mpArray in modelPointArrays:
    masses, effs = makeEffArrays(mpArray)
    massesArrs.append(masses)
    couplings.append(mpArray[0].coupling)
    effArrs.append(effs)

  rootFile.cd()
  c = TCanvas("c","c",100,100,600,600)
  c.cd()
  #c.SetLogy()
  c.SetRightMargin(0.04)
  # turn the arrays into graphs
  effGraphs = []
  for massesArr, effsArr in itertools.izip(massesArrs,effArrs):
    effGraphs.append(TGraph(len(massesArr), array.array("f",massesArr),array.array("f",effsArr)))

  SetCustomGStyle()

  # Multigraph
  mg = TMultiGraph()
  legend = TLegend(0.42,0.71,0.73,0.92)
  colorIndex = 2 
  for effGraph, coupling in itertools.izip(effGraphs,couplings):
    effGraph.SetMarkerSize(0.8)
    effGraph.SetMarkerColor(colorIndex)
    effGraph.SetLineColor(colorIndex)
    effGraph.SetLineWidth(4)
    effGraph.SetLineStyle(2)
    # no exp graph for now
    #mg.Add(limitExpGraph,"C")
    #legend.AddEntry(limitExpGraph,"Exp. limit #tilde{k} = "+str(coupling),"l")
    mg.Add(effGraph,"l")
    legend.AddEntry(effGraph," efficiency #tilde{k} = "+str(coupling),"l")
    if colorIndex==2:
      colorIndex = 4
    elif colorIndex==4:
      colorIndex = 8
    elif colorIndex>=8:
      colorIndex+=1


  # To make the separate graphs share the same axes
  h = TH1F("test","test",10,750,3500)
  h.SetStats(False)
  h.GetYaxis().SetRangeUser(0.2,0.7)
  h.GetXaxis().SetTitle("M_{1} [GeV]")
  h.GetYaxis().SetTitle("efficiency")
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

  plotname = "allEfficiencies.pdf"
  savename = TString(plotname)
  pdfName = savename.Data()
  c.SaveAs(savename.Data())
  plotname = "allEfficiencies.C"
  savename = TString(plotname)
  c.SaveAs(savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allEfficiencies.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',pdfName,savename.Data()])
  # write
  mg.SetName("allEfficiencies")
  mg.Write()


def PlotAllHalfWidths(modelPointArrays, lumi, rootFile):
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
  c.SaveAs(savename.Data())
  plotname = "allHalfWidths.C"
  savename = TString(plotname)
  c.SaveAs(savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allHalfWidths.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',pdfName,savename.Data()])
  # write
  mg.SetName("allHalfWidths")
  mg.Write()


def PlotAllExpBGs(modelPointArrays, lumi, rootFile):
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
  c.SaveAs(savename.Data())
  plotname = "allExpBGs.C"
  savename = TString(plotname)
  c.SaveAs(savename.Data())
  # png output looks strange, so convert from pdf instead
  plotname = "allExpBGs.png"
  savename = TString(plotname)
  subprocess.Popen(['convert','-trim',pdfName,savename.Data()])
  # write
  mg.SetName("allExpBGs")
  mg.Write()


