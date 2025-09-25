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
        epsilon = 1e-6
        if phi > np.pi:
            phi = np.pi - epsilon
        elif phi < -np.pi:
            phi = -np.pi + epsilon
        return phi

    def analyze(self, event):
        fatjets = Collection(event, "FatJet")
        muons   = Collection(event, "Muon")
        veto_flag = False

        for fatjet in fatjets:
            if (fatjet.pt > 170
                and (fatjet.jetId == 2 or fatjet.jetId == 6)  # requires FatJet_jetId branch
                and (fatjet.chEmEF + fatjet.neEmEF) < 0.9):

                fatjet_p4 = ROOT.TLorentzVector()
                fatjet_p4.SetPtEtaPhiM(fatjet.pt, fatjet.eta, fatjet.phi, fatjet.mass)

                overlap = False
                for mu in muons:
                    mu_p4 = ROOT.TLorentzVector()
                    mu_p4.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.mass)
                    if fatjet_p4.DeltaR(mu_p4) < 0.2:
                        overlap = True
                        break
                if overlap:
                    continue

                phi = self.fixPhi(fatjet.phi)
                veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, fatjet.eta, phi)
                if veto_map_value > 0:
                    veto_flag = True
                    break

        self.out.fillBranch("Flag_FatJetVetoed", veto_flag)
        return True
