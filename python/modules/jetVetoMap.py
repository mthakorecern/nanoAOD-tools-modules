"""Add jet veto map flag to jets.
See example in test/example_jetVeto.py for usage.
"""
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

class jetVMAP(Module):
    def __init__(self, json_JVMAP, corrName=None, veto_map_name = "jetvetomap"):
        """Module to apply veto maps that veto out events with important jets in the "hot" or 
        "cold" zones. These maps should be applied similarly both on Data and MC, to keep the
        phase-spaces equal. 
        Parameters:
        """
        self.veto_map_name = veto_map_name
        self.evaluator_JVMAP= correctionlib.CorrectionSet.from_file(json_JVMAP)
        self.evaluator_VETO = self.evaluator_JVMAP[corrName]
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Flag_JetVetoed", "O", title="Event veto flag from Jet Veto Map")
    
    def fixPhi(self, phi):
        epsilon = 1e-6  # Small offset to avoid boundary issues
        if phi > np.pi:
            #print(f"phi {phi} is greater than pi. Setting phi to pi - epsilon.")
            phi = np.pi - epsilon
        elif phi < -np.pi:
            #print(f"phi {phi} is less than -pi. Setting phi to -pi + epsilon.")
            phi = -np.pi + epsilon
        return phi

    def analyze(self, event):
        
        '''nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        - jets that don't overlap with PF muon (dR < 0.2)
        '''     
        jets  = Collection(event, "Jet")
        muons = Collection(event, "Muon")
        veto_flag = False

        for jet in jets:
            if (jet.pt > 15 
                and (jet.jetId == 2 or jet.jetId == 6) 
                and (jet.chEmEF + jet.neEmEF) < 0.9):

                # Build jet TLorentzVector
                jet_p4 = ROOT.TLorentzVector()
                jet_p4.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.mass)

                # Check overlap with all muons
                overlap = False
                for mu in muons:
                    mu_p4 = ROOT.TLorentzVector()
                    mu_p4.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.mass)
                    if jet_p4.DeltaR(mu_p4) < 0.2:
                        overlap = True
                        break

                if overlap:
                    continue  # skip this jet, since it overlaps with a muon

                # If jet passes selection and does not overlap → apply veto map
                phi = self.fixPhi(jet.phi)
                veto_map_value = self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi)
                if veto_map_value > 0:
                    veto_flag = True
                    break

        self.out.fillBranch("Flag_JetVetoed", veto_flag)
        return True
