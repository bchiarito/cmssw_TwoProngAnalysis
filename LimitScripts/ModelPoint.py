#!/usr/bin/env python

#
# Define ModelPoint class
#
# Seth I. Cooper, U. Alabama
# November 21 2012


class ModelPoint:
  def __init__(self,coupling=-1,mass=-1,totalXSec=-1,totalEff=-1,nDataObs=-1,bg=-1,bgErr=-1):
    self.coupling = coupling
    self.mass = mass
    self.totalXSec = totalXSec
    self.totalEff = totalEff
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
