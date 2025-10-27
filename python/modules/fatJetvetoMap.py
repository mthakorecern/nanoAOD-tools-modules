from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import numpy as np
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class fatJetVMAP(Module):
    def __init__(self, json_JVMAP, corrName=None, veto_map_name="jetvetomap"):
        self.veto_map_name = veto_map_name
        self.evaluator_JVMAP = correctionlib.CorrectionSet.from_file(json_JVMAP)
        self.evaluator_VETO = self.evaluator_JVMAP[corrName]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Flag_FatJetVetoed", "O", title="Event veto flag from FatJet Veto Map")

    def fixPhi(self, phi):
        return np.arctan2(np.sin(phi), np.cos(phi))

    def analyze(self, event):
        fatjets = Collection(event, "FatJet")
        veto_flag = False

        for fatjet in fatjets:
            if fatjet.pt > 170 and (fatjet.jetId & 2) and (fatjet.chEmEF + fatjet.neEmEF) < 0.9:
                phi = self.fixPhi(fatjet.phi)
                if self.evaluator_VETO.evaluate(self.veto_map_name, fatjet.eta, phi) > 0:
                    veto_flag = True
                    break

        self.out.fillBranch("Flag_FatJetVetoed", veto_flag)
        return True
