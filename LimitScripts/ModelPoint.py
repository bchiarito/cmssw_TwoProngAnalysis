#!/usr/bin/env python

#
# Define ModelPoint class
#
# Seth I. Cooper, U. Alabama
# November 21 2012
#  update for individual channels
#    July-August 2015

# important note: when modifying this class, ComputeLimit.C must also be modified
# to pass in and write out any new data members!
# NB: No longer needed with move to combine for limits

import string

class ModelPoint:
  def __init__(self, *args, **kwargs):
    self.coupling               = kwargs.get('coupling',-1)
    self.mass                   = kwargs.get('mass',-1)
    self.channel                = kwargs.get('channel','None')
    self.totalXSec              = kwargs.get('totalXSec',-1)
    self.kFactor                = kwargs.get('kFactor',1.0)
    self.acceptance             = kwargs.get('acceptance',0.0)
    self.acceptanceMassWindow   = kwargs.get('acceptanceMassWindow',0.0)
    self.totalSignalEvents      = kwargs.get('totalSignalEvents',0.0)
    self.preMWEff               = kwargs.get('preMWEff',0.0)
    self.massPeak               = kwargs.get('massPeak',-1)
    self.optSSBValue            = kwargs.get('optSSBValue',-1)
    self.totalEff               = kwargs.get('totalEff',-1)
    self.totalEffErrStat        = kwargs.get('totalEffErrStat',0.0)
    self.totalEffErrSyst        = kwargs.get('totalEffErrSyst',0.0)
    self.totalEffMScaleSystUp   = kwargs.get('totalEffMScaleSystUp',0.0)
    self.totalEffMScaleSystDown = kwargs.get('totalEffMScaleSystDown',0.0)
    self.totalEffMResSystUp     = kwargs.get('totalEffMResSystUp',0.0)
    self.totalEffMResSystDown   = kwargs.get('totalEffMResSystDown',0.0)
    self.totalEffPileupSystUp   = kwargs.get('totalEffPileupSystUp',0.0)
    self.totalEffPileupSystDown = kwargs.get('totalEffPileupSystDown',0.0)
    self.halfWidth              = kwargs.get('halfWidth',-1)
    self.optMassWindowLow       = kwargs.get('optMassWindowLow',-1)
    self.optMassWindowHigh      = kwargs.get('optMassWindowHigh',-1)
    self.nDataObs               = kwargs.get('nDataObs',-1)
    self.nBackground            = kwargs.get('nBg',-1)
    self.nBackgroundErrStat     = kwargs.get('nBgErrStat',0.0)
    self.nBackgroundErrSyst     = kwargs.get('nBgErrSyst',0.0)
    self.expLimit               = kwargs.get('expLimit',-1)
    self.expLimitOneSigmaHigh   = kwargs.get('expLimitOneSigmaHigh',-1)
    self.expLimitOneSigmaLow    = kwargs.get('expLimitOneSigmaLow',-1)
    self.expLimitTwoSigmaHigh   = kwargs.get('expLimitTwoSigmaHigh',-1)
    self.expLimitTwoSigmaLow    = kwargs.get('expLimitTwoSigmaLow',-1)
    self.obsLimit               = kwargs.get('obsLimit',-1)
    self.obsLimitErr            = kwargs.get('obsLimitErr',0.0)
    self.fileName               = kwargs.get('fileName','None')
    self.bgFileName             = kwargs.get('bgFileName','None')
    self.dFileName             = kwargs.get('dFileName','None')

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
    print "Channel: ",self.channel
    print "TotalXSection: " , self.totalXSec
    print "KFactor: " , self.kFactor
    print "Acceptance: " , self.acceptance
    print "AcceptanceMassWindow: " , self.acceptanceMassWindow
    print "Total Signal Events: " , self.totalSignalEvents
    print "Efficiency before mass window:",self.preMWEff
    print "Mass peak:",self.massPeak
    print "Opt. s/sqrt{s+b} value:",self.optSSBValue
    print "TotalEff: %.5f"%self.totalEff , " +/- %.5f (stat)"%self.totalEffErrStat, " +/- %.5f (syst)"%self.totalEffErrSyst
    print "TotalEffMScaleSystUp: %.5f"%self.totalEffMScaleSystUp
    print "TotalEffMScaleSystDown: %.5f"%self.totalEffMScaleSystDown
    print "TotalEffMResSystUp: %.5f"%self.totalEffMResSystUp
    print "TotalEffMResSystDown: %.5f"%self.totalEffMResSystDown
    print "TotalEffPileupSystUp: %.5f"%self.totalEffPileupSystUp
    print "TotalEffPileupSystDown: %.5f"%self.totalEffPileupSystDown
    print "HalfWidth: " , self.halfWidth
    print "OptMassWindow: ", self.optMassWindowLow,"-",self.optMassWindowHigh
    print "NDataObs: " , self.nDataObs
    print "NBackground: %.5f"%self.nBackground , " +/- %.5f (stat)"%self.nBackgroundErrStat, " +/- %.5f (syst)"%self.nBackgroundErrSyst
    print "Filename:",self.fileName
    print "BG MC Filename:",self.bgFileName
    print "Data Filename:",self.dFileName
    print "Expected Limit:" , self.expLimit , " + " , self.expLimitOneSigmaHigh , " - " , self.expLimitOneSigmaLow
    print "Expected Limit2SigmaBounds: + " , self.expLimitTwoSigmaHigh , " - " , self.expLimitTwoSigmaLow
    print "Observed Limit: %.5f"%self.obsLimit, " +/- %.5f "%self.obsLimitErr


  def Write(self,file):
    file.write("Coupling: " + str(self.coupling) + "\n")
    file.write("Mass: " + str(self.mass) + "\n")
    file.write("Channel: " + self.channel + "\n")
    file.write("FileName: " + self.fileName + "\n")
    file.write("BG FileName: " + self.bgFileName + "\n")
    file.write("Data FileName: " + self.dFileName + "\n")
    file.write("TotalXSection: " + str(self.totalXSec) + "\n")
    file.write("KFactor: " + str(self.kFactor) + "\n")
    file.write("Acceptance: " + str(self.acceptance) + "\n")
    file.write("AcceptanceMassWindow: " + str(self.acceptanceMassWindow) + "\n")
    file.write("TotalSignalEvents: " + str(self.totalSignalEvents) + "\n")
    file.write("PreMWEff: " + str(self.preMWEff) + "\n")
    file.write("MassPeak: " + str(self.massPeak) + "\n")
    file.write("OptSSBValue: " + str(self.optSSBValue) + "\n")
    file.write("TotalEff: " + str(self.totalEff) + "\n")
    file.write("TotalEffErrStat: " + str(self.totalEffErrStat) + "\n")
    file.write("TotalEffErrSyst: " + str(self.totalEffErrSyst) + "\n")
    file.write("TotalEffMScaleSystUp: " + str(self.totalEffMScaleSystUp) + "\n")
    file.write("TotalEffMScaleSystDown: " + str(self.totalEffMScaleSystDown) + "\n")
    file.write("TotalEffMResSystUp: " + str(self.totalEffMResSystUp) + "\n")
    file.write("TotalEffMResSystDown: " + str(self.totalEffMResSystDown) + "\n")
    file.write("TotalEffPileupSystUp: " + str(self.totalEffPileupSystUp) + "\n")
    file.write("TotalEffPileupSystDown: " + str(self.totalEffPileupSystDown) + "\n")
    file.write("HalfWidth: " + str(self.halfWidth) + "\n")
    file.write("OptMassWindowLow: " + str(self.optMassWindowLow) + "\n")
    file.write("OptMassWindowHigh: " + str(self.optMassWindowHigh) + "\n")
    file.write("NDataObs: " + str(self.nDataObs) + "\n")
    file.write("NBackground: " + str(self.nBackground) + "\n")
    file.write("NBackgroundErrStat: " + str(self.nBackgroundErrStat) + "\n")
    file.write("NBackgroundErrSyst: " + str(self.nBackgroundErrSyst) + "\n")
    file.write("ExpectedLimit: " + str(self.expLimit) + "\n")
    file.write("ExpectedLimitOneSigmaHigh: " + str(self.expLimitOneSigmaHigh) + "\n")
    file.write("ExpectedLimitOneSigmaLow: " + str(self.expLimitOneSigmaLow) + "\n")
    file.write("ExpectedLimitTwoSigmaHigh: " + str(self.expLimitTwoSigmaHigh) + "\n")
    file.write("ExpectedLimitTwoSigmaLow: " + str(self.expLimitTwoSigmaLow) + "\n")
    file.write("ObservedLimit: " + str(self.obsLimit) + "\n")
    file.write("ObservedLimitErr: " + str(self.obsLimitErr) + "\n\n")


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
    latexLine+='%.4f'%self.nBackground+' $\\pm$ '+'%.4f'%self.nBackgroundErrStat+' (stat) $\\pm$ '+'%.4f'%self.nBackgroundErrSyst+' (syst) & '
    latexLine+=str(int(self.nDataObs))#+'&'
    #latexLine+='%.1E'%float(self.totalXSec)+'&'
    #latexLine+='%.1E'%self.expLimit+'&'
    #latexLine+='%.1E'%self.obsLimit+'&'
    latexLine+=' \\\\'
    return latexLine


  def StringTableLine(self,lumi):
    tableString=string.ljust(str(self.coupling),6)+string.ljust(str(int(self.mass)),6)
    tableString+=string.ljust(self.channel,10)
    massWindowString=str(int(self.optMassWindowLow))+'-'+str(int(self.optMassWindowHigh))
    tableString+=string.ljust(massWindowString,11)
    tableString+=string.ljust('%0.3f'%self.totalEff,8)
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    tableString+=string.ljust(str(round(expectedSignalEvents,3)),10)
    #backExpString = '%0.4f'%self.nBackground+'+/-'+'%0.4f'%self.nBackgroundErr
    #backExpString = str(self.nBackground)+'+/-'+str(self.nBackgroundErr)
    if(self.nBackground < 1):
      backExpString = '%s' % float('%.2g' % self.nBackground)+'+/-'+'%s'%float('%.2g' % self.nBackgroundErrStat)+' (stat) +/-'+'%s'%float('%.2g' % self.nBackgroundErrSyst)+' (syst)'
    else:
      backExpString = '%0.2f'%self.nBackground+'+/-'+'%0.2f'%self.nBackgroundErrStat+' (stat) +/-'+'%0.2f'%self.nBackgroundErrSyst+' (syst)'
    tableString+=string.ljust(backExpString,41)+string.center(str(self.nDataObs),10)
    tableString+=string.center('%.1E'%float(self.totalXSec),10)
    try:
      tableString+=string.center(str(round(self.expLimit,6)),10)
      tableString+=string.center(str(round(self.obsLimit,6)),12)
    except TypeError:
      pass
    return tableString


  def TwikiTableLine(self,lumi):
    tableString='|'+string.ljust(str(self.coupling),9)+'|'+string.ljust(str(int(self.mass)),6)+'|'
    tableString+=string.ljust(self.channel,6)+'|'
    if not self.optMassWindowLow is None:
      tableString+=str(int(self.optMassWindowLow))+' to '+str(int(self.optMassWindowHigh))
    else:
      tableString+=str(int(self.mass-3*self.halfWidth))+' to '+str(int(self.mass+3*self.halfWidth))
    tableString+='|'+string.ljust('%0.3f'%self.totalEff,8)+'|'
    expectedSignalEvents = lumi * self.totalXSec * self.totalEff
    tableString+=string.center(str(round(expectedSignalEvents,2)),9)+'|'
    if(self.nBackground < 1):
      backExpString = '%s' % float('%.2g' % self.nBackground)+'+/-'+'%s'%float('%.2g' % self.nBackgroundErrStat)+' (stat) +/-'+'%s'%float('%.2g' % self.nBackgroundErrSyst)+' (syst)'
    else:
      backExpString = '%0.2f'%self.nBackground+'+/-'+'%0.2f'%self.nBackgroundErrStat+' (stat) +/-'+'%0.2f'%self.nBackgroundErrSyst+' (syst)'
    tableString+=string.center(backExpString,18)+'|'+string.center(str(self.nDataObs),10)+'|'
    try:
      tableString+=string.center(str(round(self.expLimit,6)),10)+'|'
      tableString+=string.center(str(round(self.obsLimit,6)),12)+'|'
    except TypeError:
      pass
    return tableString


def ReadFromLines(lines):
  outputModelPoints = []
  for line in lines:
    line.strip('\n')
    if len(line) < 2:
      continue
    try:
      value = float(line.split(': ')[1].rstrip())
    except ValueError:
      try:
        value = str(line.split(': ')[1]).rstrip()
      except ValueError:
        value = None
    if "Coupling:" in line:
      mp = ModelPoint()
      mp.coupling = value
    elif "Mass:" in line:
      mp.mass = value
    elif "Channel" in line:
      mp.channel = value.strip()
    elif "FileName" in line:
      if "BG" in line:
        mp.bgFileName = value
      elif "Data" in line:
        mp.dFileName = value
      else:
        mp.fileName = value
    elif "TotalXSection:" in line:
      mp.totalXSec = value
    elif "KFactor:" in line:
      mp.kFactor = value
    elif "Acceptance:" in line:
      mp.acceptance = value
    elif "AcceptanceMassWindow:" in line:
      mp.acceptanceMassWindow = value
    elif "TotalSignalEvents:" in line:
      mp.totalSignalEvents = value
    elif "PreMWEff:" in line:
      mp.preMWEff = value
    elif "MassPeak:" in line:
      mp.massPeak = value
    elif "OptSSBValue:" in line:
      mp.optSSBValue = value
    elif "TotalEff:" in line:
      mp.totalEff = value
    elif "TotalEffErrStat:" in line:
      mp.totalEffErrStat = value
    elif "TotalEffErrSyst:" in line:
      mp.totalEffErrSyst = value
    elif "TotalEffMScaleSystUp:" in line:
      mp.totalEffMScaleSystUp = value
    elif "TotalEffMScaleSystDown:" in line:
      mp.totalEffMScaleSystDown = value
    elif "TotalEffMResSystUp:" in line:
      mp.totalEffMResSystUp = value
    elif "TotalEffMResSystDown:" in line:
      mp.totalEffMResSystDown = value
    elif "TotalEffPileupSystUp:" in line:
      mp.totalEffPileupSystUp = value
    elif "TotalEffPileupSystDown:" in line:
      mp.totalEffPileupSystDown = value
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
    elif "NBackgroundErrStat:" in line:
      mp.nBackgroundErrStat = value
    elif "NBackgroundErrSyst:" in line:
      mp.nBackgroundErrSyst = value
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
    elif "ObservedLimitErr" in line:
      mp.obsLimitErr = value
      outputModelPoints.append(mp) # since this is the last var written for this MP
    elif "ObservedLimit:" in line:
      mp.obsLimit = value
  return outputModelPoints


def ReadFromFile(file):
  lines = file.readlines()
  return ReadFromLines(lines)



