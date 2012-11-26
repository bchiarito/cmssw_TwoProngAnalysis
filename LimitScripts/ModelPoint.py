#!/usr/bin/env python

#
# Define ModelPoint class
#
# Seth I. Cooper, U. Alabama
# November 21 2012

import string

class ModelPoint:
  def __init__(self,coupling=-1,mass=-1,totalXSec=-1,totalEff=-1,halfWidth=-1,nDataObs=-1,bg=-1,bgErr=-1):
    self.coupling = coupling
    self.mass = mass
    self.totalXSec = totalXSec
    self.totalEff = totalEff
    self.halfWidth = halfWidth
    self.nDataObs = nDataObs
    self.nBackground = bg
    self.nBackgroundErr = bgErr
    self.expLimit = -1
    self.expLimitOneSigmaHigh = -1
    self.expLimitOneSigmaLow = -1
    self.expLimitTwoSigmaHigh = -1
    self.expLimitTwoSigmaLow = -1
    self.obsLimit = -1

  def AddLimitResult(self,lr):
    self.expLimit = lr.GetExpectedLimit()
    self.expLimitOneSigmaHigh = lr.GetOneSigmaHighRange()
    self.expLimitOneSigmaLow = lr.GetOneSigmaLowRange()
    self.expLimitTwoSigmaHigh = lr.GetTwoSigmaHighRange()
    self.expLimitTwoSigmaLow = lr.GetTwoSigmaLowRange()
    self.obsLimit = lr.GetObservedLimit()

  def Print(self):
    print "Coupling : " , self.coupling
    print "Mass: ",self.mass
    print "TotalXSection: " , self.totalXSec
    print "TotalEff: " , self.totalEff
    print "HalfWidth: " , self.halfWidth
    print "NDataObs: " , self.nDataObs
    print "NBackground: " , self.nBackground , " +/- " , self.nBackgroundErr
    print "Expected Limit:" , self.expLimit , " + " , self.expLimitOneSigmaHigh , " - " , self.expLimitOneSigmaLow
    print "Expected Limit2SigmaBounds: + " , self.expLimitTwoSigmaHigh , " - " , self.expLimitTwoSigmaLow
    print "Observed Limit: " , self.obsLimit
    print
    print


  def Write(self,file):
    file.write("Coupling: " + str(self.coupling) + "\n")
    file.write("Mass: " + str(self.mass) + "\n")
    file.write("TotalXSection: " + str(self.totalXSec) + "\n")
    file.write("TotalEff: " + str(self.totalEff) + "\n")
    file.write("NDataObs: " + str(self.nDataObs) + "\n")
    file.write("NBackground: " + str(self.nBackground) + "\n")
    file.write("NBackgroundErr: " + str(self.nBackgroundErr) + "\n")
    file.write("ExpectedLimit: " + str(self.expLimit) + "\n")
    file.write("ExpectedLimitOneSigmaHigh: " + str(self.expLimitOneSigmaHigh) + "\n")
    file.write("ExpectedLimitOneSigmaLow: " + str(self.expLimitOneSigmaLow) + "\n")
    file.write("ExpectedLimitTwoSigmaHigh: " + str(self.expLimitTwoSigmaHigh) + "\n")
    file.write("ExpectedLimitTwoSigmaLow: " + str(self.expLimitTwoSigmaLow) + "\n")
    file.write("ObservedLimit: " + str(self.obsLimit) + "\n\n")

  def LatexTableLine(self,lumi):
    latexLine=self.coupling+'&'+str(self.mass)+str(round(self.totalEff,2))+'&'
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    latexLine+=str(int(expectedSignalEvents))
    latexLine+=str(round(self.nBackground,2))+'$\pm$'+str(round(self.nBackground,2))+'&'
    latexLine+=str(int(self.nDataObs))+'&'
    latexLine+='%.1E'%float(self.totalXSec)+'&'
    latexLine+='%.1E'%self.expLimit+'&'
    latexLine+='%.1E'%self.obsLimit+'&'
    latexLine+='\\\\'
    return latexLine


  def StringTableLine(self,lumi):
    tableString=string.ljust(str(self.coupling),9)+string.ljust(str(int(self.mass)),6)+string.ljust('%0.3f'%self.totalEff,8)
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    tableString+=string.center(str(round(expectedSignalEvents,2)),9)
    backExpString = '%0.2f'%self.nBackground+'+/-'+'%0.2f'%self.nBackgroundErr
    tableString+=string.center(backExpString,18)+string.center(str(self.nDataObs),10)
    tableString+=string.center('%.1E'%float(self.totalXSec),10)+string.center(str(round(self.expLimit,6)),10)
    tableString+=string.center(str(round(self.obsLimit,5)),12)
    return tableString
