#!/usr/bin/env python3
import os
import json
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class XSWeightOnly(Module):
    def __init__(self, filename, year=2016, isData=False):
        self.year = year
        self.isData = isData
        self.isMC = not isData
        self.filename = filename
        print("ACTIVATE: XSWeightOnly for file:", self.filename, " year:", self.year)

        # Define integrated luminosity in pb^-1 internally
        if self.year == "2024":
            self.LHCLumi = 109.08e3
        else:
            self.LHCLumi = 1.0
        print("Set Lumi =", self.LHCLumi, " (pb^-1 units stored internally)")

    def beginJob(self):
        if self.isMC:
            self.xsJSON = "2024_Samples.json"
            with open(self.xsJSON, "r") as xsjson:
                self.xsInfo = json.load(xsjson)
            print("Loaded XS JSON from:", self.xsJSON)

            self.countsJSON = "2024_weight.json"
            with open(self.countsJSON, "r") as countsjson:
                self.countsInfo = json.load(countsjson)
            print("Loaded gen weight JSON from:", self.countsJSON)

    def endJob(self):
        print("Finished job for year", self.year)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.isMC:
            print("BeginFile called. Processing input file:", inputFile.GetName())

            filenameForJson = (inputFile.GetName().strip().split("/")[-1]).split(".")[-2]
            self.filename = filenameForJson

            self.XS = self.xsInfo[self.filename]["XS"]  # convert pb to fb
            self.totalNumberOfEvents = self.countsInfo[self.filename]

            print("Dataset:", self.filename)
            print("XS =", self.XS, " fb")
            print("Sum of gen weights =", self.totalNumberOfEvents)

            # Branches for XS-only normalization
            self.out.branch("xsWeight", "F")
            self.out.branch("xs", "F")
            self.out.branch("lhclumi", "F")
            self.out.branch("lhclumi_pb", "F")
            self.out.branch("sumofgenwts", "F")

    def analyze(self, event):
        if self.isMC:
            xsWeight = (self.XS * self.LHCLumi * event.genWeight) / self.totalNumberOfEvents

            # Debugging per event (optional, can be commented out for large jobs)
            # print("Event genWeight:", event.genWeight, " xsWeight:", xsWeight)

            self.out.fillBranch("xsWeight", xsWeight)
            self.out.fillBranch("xs", self.XS)
            self.out.fillBranch("lhclumi", self.LHCLumi * 1e-3)  # fb^-1
            self.out.fillBranch("lhclumi_pb", self.LHCLumi)       # pb^-1
            self.out.fillBranch("sumofgenwts", self.totalNumberOfEvents)
        return True
