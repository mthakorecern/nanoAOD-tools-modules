#!/usr/bin/env python3
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NATModules.modules.tools import deltaR


class LeptonTauBranches(Module):
    """
    Extend NanoAOD with global lepton branches and per-tau/per-boostedTau matched
    electron and muon branches. Provides Leading, SubLeading, and SubSubLeading.
    """

    def __init__(self):
        self.drEle = 0.3
        self.drMu = 0.4
        self.m_ele = 0.000511
        self.m_mu = 0.10566

    def sort_by_pt(self, leptons):
        return sorted(leptons, key=lambda l: l.pt, reverse=True)

    def match_leptons_to_tau(self, tau, leptons, max_dr):
        """Return up to 3 leptons within deltaR < max_dr sorted by pT."""
        matched = []
        for lep in leptons:
            dr = deltaR(tau.eta, tau.phi, lep.eta, lep.phi)
            if dr < max_dr:
                matched.append((lep, dr))
        matched.sort(key=lambda x: x[0].pt, reverse=True)
        return matched[:3]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # Global event-level branches
        for prefix in ["Electron", "Muon"]:
            for rank in ["Leading", "SubLeading", "SubSubLeading"]:
                self.out.branch(f"{rank}{prefix}Pt", "F", title=f"pT of {rank} {prefix} in event")
                self.out.branch(f"{rank}{prefix}Eta", "F", title=f"eta of {rank} {prefix} in event")
                self.out.branch(f"{rank}{prefix}Phi", "F", title=f"phi of {rank} {prefix} in event")
                self.out.branch(f"{rank}{prefix}M", "F", title=f"mass of {rank} {prefix}")
                self.out.branch(f"{rank}{prefix}CorrIso", "F", title=f"corrected isolation of {rank} {prefix}")

        # Tau and boostedTau matched branches
        for tauType, lenVar in [("Tau", "nTau"), ("boostedTau", "nboostedTau")]:
            for prefix in ["Electron", "Muon"]:
                for rank in ["Leading", "SubLeading", "SubSubLeading"]:
                    self.out.branch(f"{tauType}_{rank}{prefix}Pt", "F", lenVar=lenVar,
                                    title=f"pT of {rank} {prefix} matched to {tauType}")
                    self.out.branch(f"{tauType}_{rank}{prefix}Eta", "F", lenVar=lenVar,
                                    title=f"eta of {rank} {prefix} matched to {tauType}")
                    self.out.branch(f"{tauType}_{rank}{prefix}Phi", "F", lenVar=lenVar,
                                    title=f"phi of {rank} {prefix} matched to {tauType}")
                    self.out.branch(f"{tauType}_{rank}{prefix}M", "F", lenVar=lenVar,
                                    title=f"mass of {rank} {prefix} matched to {tauType}")
                    self.out.branch(f"{tauType}_{rank}{prefix}CorrIso", "F", lenVar=lenVar,
                                    title=f"corrected isolation of {rank} {prefix} matched to {tauType}")
                    self.out.branch(f"{tauType}_{rank}{prefix}delR", "F", lenVar=lenVar,
                                    title=f"deltaR between {tauType} and {rank} {prefix}")

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        btaus = Collection(event, "boostedTau")

        # Global event-level electrons
        sorted_ele = self.sort_by_pt(electrons)
        for i, rank in enumerate(["Leading", "SubLeading", "SubSubLeading"]):
            if i < len(sorted_ele):
                ele = sorted_ele[i]
                self.out.fillBranch(f"{rank}ElectronPt", ele.pt)
                self.out.fillBranch(f"{rank}ElectronEta", ele.eta)
                self.out.fillBranch(f"{rank}ElectronPhi", ele.phi)
                self.out.fillBranch(f"{rank}ElectronM", self.m_ele)
                self.out.fillBranch(f"{rank}ElectronCorrIso", ele.pfRelIso03_all)
            else:
                for v in ["Pt", "Eta", "Phi", "M", "CorrIso"]:
                    self.out.fillBranch(f"{rank}Electron{v}", -99.0)

        # Global event-level muons
        sorted_mu = self.sort_by_pt(muons)
        for i, rank in enumerate(["Leading", "SubLeading", "SubSubLeading"]):
            if i < len(sorted_mu):
                mu = sorted_mu[i]
                self.out.fillBranch(f"{rank}MuonPt", mu.pt)
                self.out.fillBranch(f"{rank}MuonEta", mu.eta)
                self.out.fillBranch(f"{rank}MuonPhi", mu.phi)
                self.out.fillBranch(f"{rank}MuonM", self.m_mu)
                self.out.fillBranch(f"{rank}MuonCorrIso", mu.pfRelIso04_all)
            else:
                for v in ["Pt", "Eta", "Phi", "M", "CorrIso"]:
                    self.out.fillBranch(f"{rank}Muon{v}", -99.0)

        # Tau-electron matches
        tau_e_branches = {f"{rank}{v}": [] for rank in ["Leading","SubLeading","SubSubLeading"]
                          for v in ["ElectronPt","ElectronEta","ElectronPhi","ElectronM",
                                    "ElectronCorrIso","ElectrondelR"]}
        for tau in taus:
            matches = self.match_leptons_to_tau(tau, electrons, self.drEle)
            for j, rank in enumerate(["Leading","SubLeading","SubSubLeading"]):
                if j < len(matches):
                    ele, dr = matches[j]
                    tau_e_branches[f"{rank}ElectronPt"].append(ele.pt)
                    tau_e_branches[f"{rank}ElectronEta"].append(ele.eta)
                    tau_e_branches[f"{rank}ElectronPhi"].append(ele.phi)
                    tau_e_branches[f"{rank}ElectronM"].append(self.m_ele)
                    tau_e_branches[f"{rank}ElectronCorrIso"].append(ele.pfRelIso03_all)
                    tau_e_branches[f"{rank}ElectrondelR"].append(dr)
                else:
                    for v in ["Pt","Eta","Phi","M","CorrIso","delR"]:
                        tau_e_branches[f"{rank}Electron{v}"].append(-99.0)
        for name, arr in tau_e_branches.items():
            self.out.fillBranch(f"Tau_{name}", arr)

        # Tau-muon matches
        tau_mu_branches = {f"{rank}{v}": [] for rank in ["Leading","SubLeading","SubSubLeading"]
                           for v in ["MuonPt","MuonEta","MuonPhi","MuonM",
                                     "MuonCorrIso","MuondelR"]}
        for tau in taus:
            matches = self.match_leptons_to_tau(tau, muons, self.drMu)
            for j, rank in enumerate(["Leading","SubLeading","SubSubLeading"]):
                if j < len(matches):
                    mu, dr = matches[j]
                    tau_mu_branches[f"{rank}MuonPt"].append(mu.pt)
                    tau_mu_branches[f"{rank}MuonEta"].append(mu.eta)
                    tau_mu_branches[f"{rank}MuonPhi"].append(mu.phi)
                    tau_mu_branches[f"{rank}MuonM"].append(self.m_mu)
                    tau_mu_branches[f"{rank}MuonCorrIso"].append(mu.pfRelIso04_all)
                    tau_mu_branches[f"{rank}MuondelR"].append(dr)
                else:
                    for v in ["Pt","Eta","Phi","M","CorrIso","delR"]:
                        tau_mu_branches[f"{rank}Muon{v}"].append(-99.0)
        for name, arr in tau_mu_branches.items():
            self.out.fillBranch(f"Tau_{name}", arr)

        # boostedTau-electron matches
        btau_e_branches = {f"{rank}{v}": [] for rank in ["Leading","SubLeading","SubSubLeading"]
                           for v in ["ElectronPt","ElectronEta","ElectronPhi","ElectronM",
                                     "ElectronCorrIso","ElectrondelR"]}
        for btau in btaus:
            matches = self.match_leptons_to_tau(btau, electrons, self.drEle)
            for j, rank in enumerate(["Leading","SubLeading","SubSubLeading"]):
                if j < len(matches):
                    ele, dr = matches[j]
                    btau_e_branches[f"{rank}ElectronPt"].append(ele.pt)
                    btau_e_branches[f"{rank}ElectronEta"].append(ele.eta)
                    btau_e_branches[f"{rank}ElectronPhi"].append(ele.phi)
                    btau_e_branches[f"{rank}ElectronM"].append(self.m_ele)
                    btau_e_branches[f"{rank}ElectronCorrIso"].append(ele.pfRelIso03_all)
                    btau_e_branches[f"{rank}ElectrondelR"].append(dr)
                else:
                    for v in ["Pt","Eta","Phi","M","CorrIso","delR"]:
                        btau_e_branches[f"{rank}Electron{v}"].append(-99.0)
        for name, arr in btau_e_branches.items():
            self.out.fillBranch(f"boostedTau_{name}", arr)

        # boostedTau-muon matches
        btau_mu_branches = {f"{rank}{v}": [] for rank in ["Leading","SubLeading","SubSubLeading"]
                            for v in ["MuonPt","MuonEta","MuonPhi","MuonM",
                                      "MuonCorrIso","MuondelR"]}
        for btau in btaus:
            matches = self.match_leptons_to_tau(btau, muons, self.drMu)
            for j, rank in enumerate(["Leading","SubLeading","SubSubLeading"]):
                if j < len(matches):
                    mu, dr = matches[j]
                    btau_mu_branches[f"{rank}MuonPt"].append(mu.pt)
                    btau_mu_branches[f"{rank}MuonEta"].append(mu.eta)
                    btau_mu_branches[f"{rank}MuonPhi"].append(mu.phi)
                    btau_mu_branches[f"{rank}MuonM"].append(self.m_mu)
                    btau_mu_branches[f"{rank}MuonCorrIso"].append(mu.pfRelIso04_all)
                    btau_mu_branches[f"{rank}MuondelR"].append(dr)
                else:
                    for v in ["Pt","Eta","Phi","M","CorrIso","delR"]:
                        btau_mu_branches[f"{rank}Muon{v}"].append(-99.0)
        for name, arr in btau_mu_branches.items():
            self.out.fillBranch(f"boostedTau_{name}", arr)

        return True
