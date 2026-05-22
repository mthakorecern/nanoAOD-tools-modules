#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection


class TruthMatchFlags(Module):
    """
    Standalone NanoAODTools module to produce truth-origin flags for:
      - Electron
      - Muon
      - Tau
      - boostedTau

    Produced branches
    -----------------
    Electron:
      Electron_truthMatched
      Electron_truthFromTau
      Electron_truthFromHiggsTau
      Electron_truthTauAncestorIdx
      Electron_truthHiggsAncestorIdx

    Muon:
      Muon_truthMatched
      Muon_truthFromTau
      Muon_truthFromHiggsTau
      Muon_truthTauAncestorIdx
      Muon_truthHiggsAncestorIdx

    Tau:
      Tau_truthMatched
      Tau_truthIsHadronic
      Tau_truthFromHiggs
      Tau_truthHiggsAncestorIdx

    boostedTau:
      boostedTau_truthMatched
      boostedTau_truthIsHadronic
      boostedTau_truthFromHiggs
      boostedTau_truthHiggsAncestorIdx
    """

    def __init__(self, isData=False, debug=False, higgs_pdgid=25):
        self.isData = isData
        self.isMC = not isData
        self.debug = debug
        self.higgs_pdgid = abs(higgs_pdgid)

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Electron branches
        self.out.branch("Electron_truthMatched", "I", lenVar="nElectron")
        self.out.branch("Electron_truthFromTau", "I", lenVar="nElectron")
        self.out.branch("Electron_truthFromHiggsTau", "I", lenVar="nElectron")
        self.out.branch("Electron_truthTauAncestorIdx", "I", lenVar="nElectron")
        self.out.branch("Electron_truthHiggsAncestorIdx", "I", lenVar="nElectron")

        # Muon branches
        self.out.branch("Muon_truthMatched", "I", lenVar="nMuon")
        self.out.branch("Muon_truthFromTau", "I", lenVar="nMuon")
        self.out.branch("Muon_truthFromHiggsTau", "I", lenVar="nMuon")
        self.out.branch("Muon_truthTauAncestorIdx", "I", lenVar="nMuon")
        self.out.branch("Muon_truthHiggsAncestorIdx", "I", lenVar="nMuon")

        # HPS Tau branches
        self.out.branch("Tau_truthMatched", "I", lenVar="nTau")
        self.out.branch("Tau_truthIsHadronic", "I", lenVar="nTau")
        self.out.branch("Tau_truthFromHiggs", "I", lenVar="nTau")
        self.out.branch("Tau_truthHiggsAncestorIdx", "I", lenVar="nTau")

        # boostedTau branches
        self.out.branch("boostedTau_truthMatched", "I", lenVar="nboostedTau")
        self.out.branch("boostedTau_truthIsHadronic", "I", lenVar="nboostedTau")
        self.out.branch("boostedTau_truthFromHiggs", "I", lenVar="nboostedTau")
        self.out.branch("boostedTau_truthHiggsAncestorIdx", "I", lenVar="nboostedTau")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    # -----------------------------
    # Helper methods
    # -----------------------------
    def _valid_index(self, idx, coll):
        return 0 <= idx < len(coll)

    def _mother_idx(self, genparts, idx):
        if not self._valid_index(idx, genparts):
            return -1
        mom = genparts[idx].genPartIdxMother
        return mom if self._valid_index(mom, genparts) else -1

    def _walk_to_first_ancestor(self, genparts, start_idx, target_pdgid):
        """
        Walk upward from start_idx until the first ancestor with abs(pdgId)==target_pdgid is found.
        Returns ancestor index or -1.
        """
        target_pdgid = abs(target_pdgid)
        visited = set()
        idx = start_idx

        while self._valid_index(idx, genparts) and idx not in visited:
            visited.add(idx)

            if abs(genparts[idx].pdgId) == target_pdgid:
                return idx

            idx = self._mother_idx(genparts, idx)

        return -1

    def _walk_to_higgs_ancestor(self, genparts, start_idx):
        """
        Walk upward from start_idx until a Higgs ancestor is found.
        Returns Higgs ancestor index or -1.
        """
        visited = set()
        idx = start_idx

        while self._valid_index(idx, genparts) and idx not in visited:
            visited.add(idx)

            if abs(genparts[idx].pdgId) == self.higgs_pdgid:
                return idx

            idx = self._mother_idx(genparts, idx)

        return -1

    # -----------------------------
    # Electron helpers
    # -----------------------------
    def _electron_tau_ancestor_idx(self, ele, genparts):
        if not self.isMC:
            return -1
        if getattr(ele, "genPartIdx", -1) < 0:
            return -1
        if getattr(ele, "genPartFlav", 0) != 15:
            return -1

        # Start from matched gen electron and walk upward to first tau ancestor
        return self._walk_to_first_ancestor(genparts, ele.genPartIdx, 15)

    def _electron_higgs_ancestor_idx(self, ele, genparts):
        tau_idx = self._electron_tau_ancestor_idx(ele, genparts)
        if tau_idx < 0:
            return -1
        return self._walk_to_higgs_ancestor(genparts, tau_idx)

    # -----------------------------
    # Muon helpers
    # -----------------------------
    def _muon_tau_ancestor_idx(self, mu, genparts):
        if not self.isMC:
            return -1
        if getattr(mu, "genPartIdx", -1) < 0:
            return -1
        if getattr(mu, "genPartFlav", 0) != 15:
            return -1

        # Start from matched gen muon and walk upward to first tau ancestor
        return self._walk_to_first_ancestor(genparts, mu.genPartIdx, 15)

    def _muon_higgs_ancestor_idx(self, mu, genparts):
        tau_idx = self._muon_tau_ancestor_idx(mu, genparts)
        if tau_idx < 0:
            return -1
        return self._walk_to_higgs_ancestor(genparts, tau_idx)

    # -----------------------------
    # Tau / boostedTau helpers
    # -----------------------------
    def _tau_higgs_ancestor_idx(self, tau, genparts):
        """
        For reco taus matched to status==2 taus.
        Requires hadronic tau flavour (genPartFlav == 5).
        """
        if not self.isMC:
            return -1
        if getattr(tau, "genPartIdx", -1) < 0:
            return -1
        if getattr(tau, "genPartFlav", 0) != 5:
            return -1

        return self._walk_to_higgs_ancestor(genparts, tau.genPartIdx)

    # -----------------------------
    # Main analyze
    # -----------------------------
    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        boosted_taus = Collection(event, "boostedTau")
        genparts = Collection(event, "GenPart") if self.isMC else []

        # -----------------------------
        # Prepare output containers
        # -----------------------------
        electron_truthMatched = []
        electron_truthFromTau = []
        electron_truthFromHiggsTau = []
        electron_truthTauAncestorIdx = []
        electron_truthHiggsAncestorIdx = []

        muon_truthMatched = []
        muon_truthFromTau = []
        muon_truthFromHiggsTau = []
        muon_truthTauAncestorIdx = []
        muon_truthHiggsAncestorIdx = []

        tau_truthMatched = []
        tau_truthIsHadronic = []
        tau_truthFromHiggs = []
        tau_truthHiggsAncestorIdx = []

        boostedTau_truthMatched = []
        boostedTau_truthIsHadronic = []
        boostedTau_truthFromHiggs = []
        boostedTau_truthHiggsAncestorIdx = []

        # -----------------------------
        # Electrons
        # -----------------------------
        for i, ele in enumerate(electrons):
            matched = int(self.isMC and getattr(ele, "genPartIdx", -1) >= 0)
            tau_anc = self._electron_tau_ancestor_idx(ele, genparts) if self.isMC else -1
            higgs_anc = self._electron_higgs_ancestor_idx(ele, genparts) if self.isMC else -1

            from_tau = int(tau_anc >= 0)
            from_higgs_tau = int(higgs_anc >= 0)

            electron_truthMatched.append(matched)
            electron_truthFromTau.append(from_tau)
            electron_truthFromHiggsTau.append(from_higgs_tau)
            electron_truthTauAncestorIdx.append(int(tau_anc))
            electron_truthHiggsAncestorIdx.append(int(higgs_anc))

            if self.debug and from_higgs_tau:
                print(
                    "[TruthMatchFlags] Electron[{0}] genPartIdx={1} genPartFlav={2} "
                    "tauAnc={3} higgsAnc={4}".format(
                        i,
                        getattr(ele, "genPartIdx", -1),
                        getattr(ele, "genPartFlav", -1),
                        tau_anc,
                        higgs_anc,
                    )
                )

        # -----------------------------
        # Muons
        # -----------------------------
        for i, mu in enumerate(muons):
            matched = int(self.isMC and getattr(mu, "genPartIdx", -1) >= 0)
            tau_anc = self._muon_tau_ancestor_idx(mu, genparts) if self.isMC else -1
            higgs_anc = self._muon_higgs_ancestor_idx(mu, genparts) if self.isMC else -1

            from_tau = int(tau_anc >= 0)
            from_higgs_tau = int(higgs_anc >= 0)

            muon_truthMatched.append(matched)
            muon_truthFromTau.append(from_tau)
            muon_truthFromHiggsTau.append(from_higgs_tau)
            muon_truthTauAncestorIdx.append(int(tau_anc))
            muon_truthHiggsAncestorIdx.append(int(higgs_anc))

            if self.debug and from_higgs_tau:
                print(
                    "[TruthMatchFlags] Muon[{0}] genPartIdx={1} genPartFlav={2} "
                    "tauAnc={3} higgsAnc={4}".format(
                        i,
                        getattr(mu, "genPartIdx", -1),
                        getattr(mu, "genPartFlav", -1),
                        tau_anc,
                        higgs_anc,
                    )
                )

        # -----------------------------
        # HPS Taus
        # -----------------------------
        for i, tau in enumerate(taus):
            matched = int(self.isMC and getattr(tau, "genPartIdx", -1) >= 0)
            is_hadronic = int(self.isMC and getattr(tau, "genPartFlav", 0) == 5)
            higgs_anc = self._tau_higgs_ancestor_idx(tau, genparts) if self.isMC else -1
            from_higgs = int(higgs_anc >= 0)

            tau_truthMatched.append(matched)
            tau_truthIsHadronic.append(is_hadronic)
            tau_truthFromHiggs.append(from_higgs)
            tau_truthHiggsAncestorIdx.append(int(higgs_anc))

            if self.debug and from_higgs:
                print(
                    "Tau[{0}] genPartIdx={1} genPartFlav={2} "
                    "higgsAnc={3}".format(
                        i,
                        getattr(tau, "genPartIdx", -1),
                        getattr(tau, "genPartFlav", -1),
                        higgs_anc,
                    )
                )

        # -----------------------------
        # boostedTaus
        # -----------------------------
        for i, btau in enumerate(boosted_taus):
            matched = int(self.isMC and getattr(btau, "genPartIdx", -1) >= 0)
            is_hadronic = int(self.isMC and getattr(btau, "genPartFlav", 0) == 5)
            higgs_anc = self._tau_higgs_ancestor_idx(btau, genparts) if self.isMC else -1
            from_higgs = int(higgs_anc >= 0)

            boostedTau_truthMatched.append(matched)
            boostedTau_truthIsHadronic.append(is_hadronic)
            boostedTau_truthFromHiggs.append(from_higgs)
            boostedTau_truthHiggsAncestorIdx.append(int(higgs_anc))

            if self.debug and from_higgs:
                print(
                    "boostedTau[{0}] genPartIdx={1} genPartFlav={2} "
                    "higgsAnc={3}".format(
                        i,
                        getattr(btau, "genPartIdx", -1),
                        getattr(btau, "genPartFlav", -1),
                        higgs_anc,
                    )
                )

        # -----------------------------
        # Fill branches
        # -----------------------------
        self.out.fillBranch("Electron_truthMatched", electron_truthMatched)
        self.out.fillBranch("Electron_truthFromTau", electron_truthFromTau)
        self.out.fillBranch("Electron_truthFromHiggsTau", electron_truthFromHiggsTau)
        self.out.fillBranch("Electron_truthTauAncestorIdx", electron_truthTauAncestorIdx)
        self.out.fillBranch("Electron_truthHiggsAncestorIdx", electron_truthHiggsAncestorIdx)

        self.out.fillBranch("Muon_truthMatched", muon_truthMatched)
        self.out.fillBranch("Muon_truthFromTau", muon_truthFromTau)
        self.out.fillBranch("Muon_truthFromHiggsTau", muon_truthFromHiggsTau)
        self.out.fillBranch("Muon_truthTauAncestorIdx", muon_truthTauAncestorIdx)
        self.out.fillBranch("Muon_truthHiggsAncestorIdx", muon_truthHiggsAncestorIdx)

        self.out.fillBranch("Tau_truthMatched", tau_truthMatched)
        self.out.fillBranch("Tau_truthIsHadronic", tau_truthIsHadronic)
        self.out.fillBranch("Tau_truthFromHiggs", tau_truthFromHiggs)
        self.out.fillBranch("Tau_truthHiggsAncestorIdx", tau_truthHiggsAncestorIdx)

        self.out.fillBranch("boostedTau_truthMatched", boostedTau_truthMatched)
        self.out.fillBranch("boostedTau_truthIsHadronic", boostedTau_truthIsHadronic)
        self.out.fillBranch("boostedTau_truthFromHiggs", boostedTau_truthFromHiggs)
        self.out.fillBranch("boostedTau_truthHiggsAncestorIdx", boostedTau_truthHiggsAncestorIdx)

        return True



#############################################

'''
Leading selected electron from Higgs tau: ngood_Electrons > 0 && Electron_truthFromHiggsTau[index_gElectrons[0]] == 1

Leading selected muon from Higgs tau: ngood_Muons > 0 && Muon_truthFromHiggsTau[index_gMuons[0]] == 1

Leading selected HPS tau from Higgs: ngood_Taus > 0 && Tau_truthFromHiggs[index_gTaus[0]] == 1

Leading selected boosted tau from Higgs: ngood_boostedTaus > 0 && boostedTau_truthFromHiggs[index_gboostedTaus[0]] == 1

Electron and tau from the same Higgs: 
ngood_Electrons > 0 &&
ngood_Taus > 0 &&
Electron_truthFromHiggsTau[index_gElectrons[0]] == 1 &&
Tau_truthFromHiggs[index_gTaus[0]] == 1 &&
Electron_truthHiggsAncestorIdx[index_gElectrons[0]] == Tau_truthHiggsAncestorIdx[index_gTaus[0]]


Electron and boosted tau from the same Higgs:

ngood_Electrons > 0 &&
ngood_boostedTaus > 0 &&
Electron_truthFromHiggsTau[index_gElectrons[0]] == 1 &&
boostedTau_truthFromHiggs[index_gboostedTaus[0]] == 1 &&
Electron_truthHiggsAncestorIdx[index_gElectrons[0]] == boostedTau_truthHiggsAncestorIdx[index_gboostedTaus[0]]

Muon and tau from the same Higgs:
ngood_Muons > 0 &&
ngood_Taus > 0 &&
Muon_truthFromHiggsTau[index_gMuons[0]] == 1 &&
Tau_truthFromHiggs[index_gTaus[0]] == 1 &&
Muon_truthHiggsAncestorIdx[index_gMuons[0]] == Tau_truthHiggsAncestorIdx[index_gTaus[0]]

Muon and boosted tau from the same Higgs:
ngood_Muons > 0 &&
ngood_boostedTaus > 0 &&
Muon_truthFromHiggsTau[index_gMuons[0]] == 1 &&
boostedTau_truthFromHiggs[index_gboostedTaus[0]] == 1 &&
Muon_truthHiggsAncestorIdx[index_gMuons[0]] == boostedTau_truthHiggsAncestorIdx[index_gboostedTaus[0]]

'''