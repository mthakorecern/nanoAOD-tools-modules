"""
Add FatJet_jetId variable, following the style of jetId module for AK4 jets.
See jsonpog-integration jetID JSON files for supported keys.
"""

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import correctionlib
from array import array

class fatJetId(Module):
    def __init__(self, json, jetType="AK8PUPPI"):
        """Module to determine FatJetID variables (passTight, passTightLepVeto),
        packed in FatJet_jetId branch like in NanoAODv12.
        
        Parameters:
        - json: jetID json file (from jsonpog-integration)
        - jetType: default "AK8PUPPI" (alternatives may exist if JSON provides)
        """
        self.evaluator = correctionlib.CorrectionSet.from_file(json)
        self.key_tight = f"{jetType}_Tight"
        self.key_tightLeptonVeto = f"{jetType}_TightLeptonVeto"

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        # Add branch for FatJet_jetId, same encoding: bit2 = tight, bit3 = tightLepVeto
        self.out.branch("FatJet_jetId", "b", lenVar="nFatJet",
                        title="FatJet ID flag: bit2 is tight, bit3 is tightLepVeto (recomputed using JSON)")

    def analyze(self, event):
        fatjets = Collection(event, "FatJet")

        fatjet_Ids = array('B', event.nFatJet * [0])  # UChar_t is 'B'
        for ijet, jet in enumerate(fatjets):
            multiplicity = jet.chMultiplicity + jet.neMultiplicity

            passTight = self.evaluator[self.key_tight].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )

            passTightLepVeto = self.evaluator[self.key_tightLeptonVeto].evaluate(
                jet.eta,
                jet.chHEF,
                jet.neHEF,
                jet.chEmEF,
                jet.neEmEF,
                jet.muEF,
                jet.chMultiplicity,
                jet.neMultiplicity,
                multiplicity
            )

            fatjet_Ids[ijet] = int(passTight) * 2 + int(passTightLepVeto) * 4

        self.out.fillBranch("FatJet_jetId", fatjet_Ids)
        return True
