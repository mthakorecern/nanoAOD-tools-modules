#!/usr/bin/env python3
import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from correctionlib import _core
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

class BTagWeightTool:
    """
    BTagWeightTool for AK4 jets using the UParT tagger and the `UParTAK4_kinfit` correction.
    -------------------------------------------------------------------------------
    Correction: UParTAK4_kinfit
    -------------------------------------------------------------------------------
    - Purpose: Per-jet scale factors (SFs) for b/c jets in fixed working point b-tagging for the 2024_Summer24 campaign using the UParT algorithm.
    - Inputs:
        1. systematic     (str): 'central', 'fsrdef', 'isrdef', 'hdamp', 'jer', 'jes', 'mass', 'tune', 'statistic'
        2. working_point  (str): 'L', 'M', 'T', 'XT', 'XXT'
        3. flavor         (int): 0 (light), 4 (charm), 5 (bottom)
        4. abseta         (float): absolute pseudorapidity of the jet
        5. pt            (float): transverse momentum of the jet
    - Output: scale factor (float)
    -------------------------------------------------------------------------------
    Correction: UParTAK4_wp_values
    -------------------------------------------------------------------------------
    - Purpose: Provides discriminator thresholds for UParT working points.
    - Inputs:
        1. working_point (str): 'L', 'M', 'T', 'XT', 'XXT'
    - Output: discriminator threshold (float)
    -------------------------------------------------------------------------------
    """

    def __init__(self, wp="medium", sigmabc="central", sigmalight="central", year="2024"):
        """
        Initialize the UParT BTagWeightTool.

        Parameters
        ----------
        wp : str
            Working point. Accepted: 'loose', 'medium', 'tight', 'XT', 'XXT'
        sigmabc : str
            Systematic variation for b/c jets. Options: 'central', 'fsrdef', 'isrdef', 'hdamp', 'jer', 'jes', 'mass', 'tune', 'statistic'
        sigmalight : str
            Systematic variation for light jets. Same options as above.
        year : str
            Year label for JSON path resolution (default: '2024').
        """
        self.tagger = "UParT"
        self.wp = wp[0].upper()  
        self.sigmabc = sigmabc
        self.sigmalight = sigmalight
        self.year = year

        jsonPath = os.path.join(
            os.environ["CMSSW_BASE"],
            "src/jsonpog-integration/POG/BTV/btagging_preliminary.json.gz"
        )
        if not os.path.exists(jsonPath):
            raise RuntimeError(f"[BTagWeightTool] JSON not found: {jsonPath}")

        # Load correctionlib objects
        self.cset = _core.CorrectionSet.from_file(jsonPath)
        self.corr = self.cset["UParTAK4_kinfit"]
        self.wp_values = self.cset["UParTAK4_wp_values"]
        self.wp_threshold = self.wp_values.evaluate(self.wp)

        print(
            f"[BTagWeightTool] Loaded UParTAK4_kinfit for {year} | WP={self.wp}, "
            f"sigmabc={self.sigmabc}, sigmalight={self.sigmalight}"
        )

    # --------------------------------------------------------------------------
    # Per-jet SF evaluation
    # --------------------------------------------------------------------------
    def getSF(self, jet, sys=None):
        """
        Compute the per-jet scale factor for a single AK4 jet.

        Parameters
        ----------
        jet : NanoAODTools Jet
            Jet object with attributes: pt, eta, hadronFlavour.
        sys : str, optional
            Override systematic label for this evaluation.

        Returns
        -------
        float
            Scale factor for the jet.
        """
        systematic = self.sigmabc if abs(jet.hadronFlavour) in [4, 5] else self.sigmalight
        if sys is not None:
            systematic = sys

        return self.corr.evaluate(
            systematic,
            self.wp,
            jet.hadronFlavour,
            abs(jet.eta),
            jet.pt
        )

    # --------------------------------------------------------------------------
    # Per-event SF weight (product of per-jet SFs)
    # --------------------------------------------------------------------------
    def getWeight(self, jets, sys=None):
        """
        Compute the per-event b-tagging scale factor (product over all jets).

        Parameters
        ----------
        jets : list of Jet
            List of NanoAODTools Jet objects.
        sys : str, optional
            Override the systematic globally for all jets.

        Returns
        -------
        float
            Event-level b-tagging scale factor.
        """
        if not jets:
            return 1.0

        weight = 1.0
        for jet in jets:
            if jet.pt < 20 or abs(jet.eta) > 2.5:
                continue
            if abs(jet.hadronFlavour) not in [0, 4, 5]:
                continue
            weight *= self.getSF(jet, sys=sys)
        return weight

    # --------------------------------------------------------------------------
    # Convenience: directly compute from event
    # --------------------------------------------------------------------------
    def getEventWeight(self, event, jetCollectionName="Jet", sys=None):
        """
        Compute the event-level scale factor directly from a NanoAODTools event.

        Parameters
        ----------
        event : Event 
        jetCollectionName : str
            Jet collection name (default: "Jet").
        sys : str, optional
            Override systematic label.

        Returns
        -------
        float
            Event-level SF weight.
        """
        jets = Collection(event, jetCollectionName)
        return self.getWeight(jets, sys=sys)

    # --------------------------------------------------------------------------
    # Discriminator threshold check
    # --------------------------------------------------------------------------
    def isTagged(self, jet):
        """
        Check whether a jet passes the UParT discriminator threshold for this WP.

        Parameters
        ----------
        jet : Jet
            NanoAODTools Jet object.

        Returns
        -------
        bool
            True if the jet is tagged.
        """
        return jet.btagUParTAK4B > self.wp_threshold
