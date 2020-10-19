# 
#
# Read from xml and insert into database using
# PopCon 
#
# This is a template, generate real test using
#
# sed 's/EcalGainRatios/your-record/g' testTemplate.py > testyourrecord.py
#
# Stefano Argiro', $Id: testEcalGainRatios.py,v 1.1 2008/11/14 15:46:03 argiro Exp $
#
#

import FWCore.ParameterSet.Config as cms
import os, sys, imp, re
import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("TEST")

options = VarParsing.VarParsing('standard')
options.register('xmlFile',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "xmlFile")
options.register('sqliteFile',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "sqliteFile")
options.register('firstRun',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 "firstRun")
options.register('tag',
                 "",
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "tag")

options.parseArguments()
#print options
print "---> Input xmlFile    : ",options.xmlFile 
print "---> Output sqliteFile: ",options.sqliteFile 
print "---> firstRun         : ",options.firstRun
print "---> tag              : ",options.tag

process.MessageLogger=cms.Service("MessageLogger",
                              destinations=cms.untracked.vstring("cout"),
                              cout=cms.untracked.PSet(
                              )
)
process.load("CondCore.CondDB.CondDB_cfi")

process.CondDB.connect = cms.string('sqlite_file:'+str(options.sqliteFile))


process.source = cms.Source("EmptyIOVSource",
    timetype = cms.string('runnumber'),
    firstValue = cms.uint64(1),
    lastValue  = cms.uint64(1),
    interval = cms.uint64(1)
)

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    process.CondDB,
    timetype = cms.untracked.string('runnumber'),
    toPut = cms.VPSet(cms.PSet(
        record = cms.string('EcalTimeCalibConstantsRcd'),
        tag = cms.string(str(options.tag))
         )),
    logconnect= cms.untracked.string('sqlite_file:logtestEcalTimeCalib.db')                                     
)

process.mytest = cms.EDAnalyzer("EcalTimeCalibConstantsAnalyzer",
    record = cms.string('EcalTimeCalibConstantsRcd'),
    loggingOn= cms.untracked.bool(True),
    SinceAppendMode=cms.bool(True),
    Source=cms.PSet(
    xmlFile = cms.untracked.string(str(options.xmlFile)),
    since = cms.untracked.int64(int(options.firstRun))
    )                            
)

process.p = cms.Path(process.mytest)




