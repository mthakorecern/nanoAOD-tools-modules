#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class Hbbtagger(Module):
    def __init__(self, filename, year="2024", isData=False):
        self.year = str(year)
        self.isData = isData
        self.isMC = not isData
        self.filename = filename

    
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("FatJet_globpart_bbvsqcd", "F")

    def analyze(self, event):
        fatjets = Collection(event, "FatJet")
        if (event.ngood_FatJets>0):
            Fatjet_globpart_bbvsqcd = (fatjets[event.index_gFatJets[0]].globalParT3_Xbb)/(fatjets[event.index_gFatJets[0]].globalParT3_Xbb + fatjets[event.index_gFatJets[0]].globalParT3_QCD)
            self.out.fillBranch("FatJet_globpart_bbvsqcd",Fatjet_globpart_bbvsqcd)
        else:
            # print("Number of Fatjets in..Nominal...is ZERO")
            self.out.fillBranch("FatJet_globpart_bbvsqcd",-99.99)
        
        return True
            


         




