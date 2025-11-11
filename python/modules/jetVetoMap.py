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
        return np.arctan2(np.sin(phi), np.cos(phi))

    def analyze(self, event):
        
        '''nominal “loose selection”
        - jet pT > 15 GeV
        - tight jet ID
        - jet EM fraction (charged + neutral) < 0.9
        '''     
        jets  = Collection(event, "Jet")
        veto_flag = False

        for jet in jets:
            if (jet.pt > 15 and (jet.jetId & 4) and (jet.chEmEF + jet.neEmEF) < 0.9):
                phi = self.fixPhi(jet.phi)
                if self.evaluator_VETO.evaluate(self.veto_map_name, jet.eta, phi) > 0:
                    veto_flag = True
                    break

        self.out.fillBranch("Flag_JetVetoed", veto_flag)
        return True
