#!/usr/bin/env python

#
# Define ModelPoint class
#
# Seth I. Cooper, U. Alabama
# November 21 2012

import string

class ModelPoint:
  def __init__(self, *args, **kwargs):
    self.coupling             = kwargs.get('coupling',None)
    self.mass                 = kwargs.get('mass',None)
    self.totalXSec            = kwargs.get('totalXSec',None)
    self.totalEff             = kwargs.get('totalEff',None)
    self.halfWidth            = kwargs.get('halfWidth',None)
    self.optMassWindowLow     = kwargs.get('optMassWindowLow',None)
    self.optMassWindowHigh    = kwargs.get('optMassWindowHigh',None)
    self.nDataObs             = kwargs.get('nDataObs',None)
    self.nBackground          = kwargs.get('nBg',None)
    self.nBackgroundErr       = kwargs.get('nBgErr',None)
    self.expLimit             = kwargs.get('expLimit',None)
    self.expLimitOneSigmaHigh = kwargs.get('expLimitOneSigmaHigh',None)
    self.expLimitOneSigmaLow  = kwargs.get('expLimitOneSigmaLow',None)
    self.expLimitTwoSigmaHigh = kwargs.get('expLimitTwoSigmaHigh',None)
    self.expLimitTwoSigmaLow  = kwargs.get('expLimitTwoSigmaLow',None)
    self.obsLimit             = kwargs.get('obsLimit',None)
    self.fileName             = kwargs.get('fileName',None)

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
    print "TotalEff: %.5f"%self.totalEff
    print "HalfWidth: " , self.halfWidth
    print "OptMassWindow: ", self.optMassWindowLow,"-",self.optMassWindowHigh
    print "NDataObs: " , self.nDataObs
    print "NBackground: %.5f"%self.nBackground , " +/- %.5f"%self.nBackgroundErr
    print "Filename:",self.fileName
    print "Expected Limit:" , self.expLimit , " + " , self.expLimitOneSigmaHigh , " - " , self.expLimitOneSigmaLow
    print "Expected Limit2SigmaBounds: + " , self.expLimitTwoSigmaHigh , " - " , self.expLimitTwoSigmaLow
    print "Observed Limit: " , self.obsLimit


  def Write(self,file):
    file.write("Coupling: " + str(self.coupling) + "\n")
    file.write("Mass: " + str(self.mass) + "\n")
    file.write("TotalXSection: " + str(self.totalXSec) + "\n")
    file.write("TotalEff: " + str(self.totalEff) + "\n")
    file.write("HalfWidth: " + str(self.halfWidth) + "\n")
    file.write("OptMassWindowLow: " + str(self.optMassWindowLow) + "\n")
    file.write("OptMassWindowHigh: " + str(self.optMassWindowHigh) + "\n")
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
    latexLine='\t\t'+str(self.coupling)+' & '+str(int(self.mass))+' & '
    #print 'numsigmas = ',numsigmas
    #print 'halfWidth = ',self.halfWidth
    #print 'mass = ',self.mass
    #print 'self.mass - numsigmas * self.halfWidth',(self.mass - numsigmas * self.halfWidth)
    #print 'minMass = ',minMass
    #print 'maxMass = ',maxMass
    if not self.optMassWindowLow is None:
      latexLine+=str(int(self.optMassWindowLow))+' to '+str(int(self.optMassWindowHigh))+' & '
    else:
      latexLine+=str(int(self.mass-3*self.halfWidth))+' to '+str(int(self.mass+3*self.halfWidth))+' & '
    latexLine+='%.2f'%self.totalEff+' & '
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    latexLine+='%.2f'%float(expectedSignalEvents)+' & '
    latexLine+='%.4f'%self.nBackground+' $\\pm$ '+'%.4f'%self.nBackgroundErr+' & '
    latexLine+=str(int(self.nDataObs))#+'&'
    #latexLine+='%.1E'%float(self.totalXSec)+'&'
    #latexLine+='%.1E'%self.expLimit+'&'
    #latexLine+='%.1E'%self.obsLimit+'&'
    latexLine+=' \\\\'
    return latexLine


  def StringTableLine(self,lumi):
    tableString=string.ljust(str(self.coupling),9)+string.ljust(str(int(self.mass)),6)+string.ljust('%0.3f'%self.totalEff,8)
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    tableString+=string.center(str(round(expectedSignalEvents,2)),9)
    #backExpString = '%0.4f'%self.nBackground+'+/-'+'%0.4f'%self.nBackgroundErr
    #backExpString = str(self.nBackground)+'+/-'+str(self.nBackgroundErr)
    if(self.nBackground < 1):
      backExpString = '%s' % float('%.2g' % self.nBackground)+'+/-'+'%s'%float('%.2g' % self.nBackgroundErr)
    else:
      backExpString = '%0.2f'%self.nBackground+'+/-'+'%0.2f'%self.nBackgroundErr
    tableString+=string.center(backExpString,18)+string.center(str(self.nDataObs),10)
    tableString+=string.center('%.1E'%float(self.totalXSec),10)+string.center(str(round(self.expLimit,6)),10)
    tableString+=string.center(str(round(self.obsLimit,5)),12)
    return tableString


  def TwikiTableLine(self,lumi):
    tableString='|'+string.ljust(str(self.coupling),9)+'|'+string.ljust(str(int(self.mass)),6)+'|'
    if not self.optMassWindowLow is None:
      tableString+=str(int(self.optMassWindowLow))+' to '+str(int(self.optMassWindowHigh))
    else:
      tableString+=str(int(self.mass-3*self.halfWidth))+' to '+str(int(self.mass+3*self.halfWidth))
    tableString+='|'+string.ljust('%0.3f'%self.totalEff,8)+'|'
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    tableString+=string.center(str(round(expectedSignalEvents,2)),9)+'|'
    if(self.nBackground < 1):
      backExpString = '%s' % float('%.2g' % self.nBackground)+'+/-'+'%s'%float('%.2g' % self.nBackgroundErr)
    else:
      backExpString = '%0.2f'%self.nBackground+'+/-'+'%0.2f'%self.nBackgroundErr
    tableString+=string.center(backExpString,18)+'|'+string.center(str(self.nDataObs),10)+'|'
    tableString+=string.center(str(round(self.expLimit,6)),10)+'|'
    tableString+=string.center(str(round(self.obsLimit,5)),12)+'|'
    return tableString


def ReadFromFile(file):
  outputModelPoints = []
  for line in file:
    if len(line) < 2:
      continue
    try:
      value = float(line.split(': ')[1])
    except ValueError:
      value = None
    if "Coupling:" in line:
      mp = ModelPoint()
      mp.coupling = value
    elif "Mass:" in line:
      mp.mass = value
    elif "TotalXSection:" in line:
      mp.totalXSec = value
    elif "TotalEff:" in line:
      mp.totalEff = value
    elif "HalfWidth:" in line:
      mp.halfWidth = value
    elif "OptMassWindowLow:" in line:
      mp.optMassWindowLow = value
    elif "OptMassWindowHigh:" in line:
      mp.optMassWindowHigh = value
    elif "NDataObs:" in line:
      mp.nDataObs = value
    elif "NBackground:" in line and not "Err" in line:
      mp.nBackground = value
    elif "NBackgroundErr:" in line:
      mp.nBackgroundErr = value
    elif "ExpectedLimit:" in line and not "Sigma" in line:
      mp.expLimit = value
    elif "ExpectedLimitOneSigmaHigh:" in line:
      mp.expLimitOneSigmaHigh = value
    elif "ExpectedLimitOneSigmaLow:" in line:
      mp.expLimitOneSigmaLow = value
    elif "ExpectedLimitTwoSigmaHigh:" in line:
      mp.expLimitTwoSigmaHigh = value
    elif "ExpectedLimitTwoSigmaLow:" in line:
      mp.expLimitTwoSigmaLow = value
    elif "ObservedLimit:" in line:
      mp.obsLimit = value
      outputModelPoints.append(mp)
  return outputModelPoints





