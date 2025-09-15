"""Add branches for fatjet JES and JER corrections with regrouped uncertainties.
See example in test/example_jetCorr.py for usage.
"""

###
# Useful CMStalk thread on implementing jet corrections: 
# https://cms-talk.web.cern.ch/t/jes-for-2022-re-reco-cde-and-prompt-fg/32873/3
# Minimal demo from JME POG:
# https://github.com/cms-jet/JECDatabase/blob/master/scripts/JERC2JSON/minimalDemo.py
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import numpy as np
import awkward as ak
import correctionlib


# -----------------------------
# Coffea-style rand_gauss
# -----------------------------
def rand_gauss(item):
    seeds = (
        ak.typetracer.length_one_if_typetracer(item).to_numpy()[[0, -1]].view("i4")
    )
    randomstate = np.random.Generator(np.random.PCG64(seeds))

    def getfunction(layout, depth, **kwargs):
        if isinstance(layout, ak.contents.NumpyArray) or not isinstance(
            layout, (ak.contents.Content,)
        ):
            return ak.contents.NumpyArray(
                randomstate.normal(size=len(layout)).astype(np.float32)
            )
        return None

    out = ak.transform(
        getfunction,
        ak.typetracer.length_zero_if_typetracer(item),
        behavior=item.behavior,
    )
    if ak.backend(item) == "typetracer":
        out = ak.Array(
            out.layout.to_typetracer(forget_length=True), behavior=out.behavior
        )
    return out


# -----------------------------
# Coffea-style jer_smear
# -----------------------------
def jer_smear(pt_gen, jetPt, etaJet, jer, jer_rand_gauss, jersf, variation, forceStochastic):
    pt_gen = pt_gen if not forceStochastic else None
    if forceStochastic:
        pt_gen = ak.without_parameters(ak.zeros_like(jetPt))

    jersmear = jer * jer_rand_gauss
    jersf = jersf[:, variation]

    deltaPtRel = (jetPt - pt_gen) / jetPt
    doHybrid = (pt_gen > 0) & (np.abs(deltaPtRel) < 3 * jer)

    detSmear = 1 + (jersf - 1) * deltaPtRel
    stochSmear = 1 + np.sqrt(np.maximum(jersf**2 - 1, 0)) * jersmear

    min_jet_pt = (1e-2) / np.cosh(etaJet)
    min_jet_pt_corr = min_jet_pt / jetPt

    smearfact = ak.where(doHybrid, detSmear, stochSmear)
    smearfact = ak.where((smearfact * jetPt) < min_jet_pt, min_jet_pt_corr, smearfact)
    return smearfact


# -----------------------------
# FatJet JERC Module
# -----------------------------
class fatJetJERC(Module):
    def __init__(self,
                 json_JERC,
                 json_JERsmear, 
                 L1Key=None,
                 L2Key=None,
                 L3Key=None,
                 L2L3Key=None,
                 scaleTotalKey=None,
                 smearKey=None,
                 JERKey=None,
                 JERsfKey=None,
                 overwritePt=False, 
                 usePhiDependentJEC=False, 
                 useRunDependentJEC=False,
                 forceStochastic=False):
        """Correct fatjets following JME POG recommendations."""

        self.overwritePt = overwritePt
        self.usePhiDependentJEC = usePhiDependentJEC
        self.useRunDependentJEC = useRunDependentJEC
        self.forceStochastic = forceStochastic

        # Load correction JSONs
        self.evaluator_JERC = correctionlib.CorrectionSet.from_file(json_JERC)
        self.evaluator_jer = correctionlib.CorrectionSet.from_file(json_JERsmear)

        self.evaluator_L1 = self.evaluator_JERC[L1Key]
        self.evaluator_L2 = self.evaluator_JERC[L2Key]
        self.evaluator_L3 = self.evaluator_JERC[L3Key]
        self.evaluator_L2L3 = self.evaluator_JERC[L2L3Key]

        self.is_mc = False
        self.evaluator_JERsmear = None
        self.evaluator_JER = None
        self.evaluator_JERsf = None
        self.evaluator_JES = None
        self.jes_sources = []
        self.use_json_smear = False

        if smearKey is not None:  # JER is MC-only
            self.evaluator_JERsmear = self.evaluator_jer[smearKey] if smearKey in self.evaluator_jer else None
            self.use_json_smear = self.evaluator_JERsmear is not None
            self.evaluator_JER = self.evaluator_JERC[JERKey]
            self.evaluator_JERsf = self.evaluator_JERC[JERsfKey]
            self.evaluator_JES = self.evaluator_JERC[scaleTotalKey]
            self.is_mc = True

            # JES regrouped uncertainty sources
            self.jes_sources = [
                k for k in self.evaluator_JERC.keys()
                if "MC_Regrouped_" in k and k.endswith("_AK4PFPuppi")
            ]

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.overwritePt:
            self.out.branch("FatJet_uncorrected_pt", "F", lenVar="nFatJet")
            self.out.branch("FatJet_uncorrected_mass", "F", lenVar="nFatJet", limitedPrecision=12)

        self.out.branch("FatJet_pt_nom", "F", lenVar="nFatJet")
        self.out.branch("FatJet_mass_nom", "F", lenVar="nFatJet", limitedPrecision=12)

        if self.is_mc:
            self.out.branch("FatJet_pt_jesTotalUp", "F", lenVar="nFatJet")
            self.out.branch("FatJet_pt_jesTotalDown", "F", lenVar="nFatJet")
            self.out.branch("FatJet_mass_jesTotalUp", "F", lenVar="nFatJet", limitedPrecision=12)
            self.out.branch("FatJet_mass_jesTotalDown", "F", lenVar="nFatJet", limitedPrecision=12)
            self.out.branch("FatJet_pt_jerUp", "F", lenVar="nFatJet")
            self.out.branch("FatJet_pt_jerDown", "F", lenVar="nFatJet")
            self.out.branch("FatJet_mass_jerUp", "F", lenVar="nFatJet", limitedPrecision=12)
            self.out.branch("FatJet_mass_jerDown", "F", lenVar="nFatJet", limitedPrecision=12)

            for src in self.jes_sources:
                name = src.split("MC_")[1].replace("_AK8PFPuppi", "")
                self.out.branch(f"FatJet_pt_jes{name}Up", "F", lenVar="nFatJet")
                self.out.branch(f"FatJet_pt_jes{name}Down", "F", lenVar="nFatJet")
                self.out.branch(f"FatJet_mass_jes{name}Up", "F", lenVar="nFatJet")
                self.out.branch(f"FatJet_mass_jes{name}Down", "F", lenVar="nFatJet")

    def analyze(self, event):
        jets = Collection(event, "FatJet")
        if self.is_mc:
            gen_jets = Collection(event, "GenJet")
            gen_jets_pt = np.array([g.pt for g in gen_jets])
            gen_jets_eta = np.array([g.eta for g in gen_jets])
            gen_jets_phi = np.array([g.phi for g in gen_jets])

        pt_corr, mass_corr = [], []
        pt_uncorr, mass_uncorr = [], []
        pt_smear_up, pt_smear_dn = [], []
        mass_smear_up, mass_smear_dn = [], []
        pt_scale_up, pt_scale_dn = [], []
        mass_scale_up, mass_scale_dn = [], []
        pt_sources_up, pt_sources_dn = {}, {}
        mass_sources_up, mass_sources_dn = {}, {}

        for jet in jets:
            pt_raw = jet.pt * (1 - jet.rawFactor)
            mass_raw = jet.mass * (1 - jet.rawFactor)

            pt_L1 = pt_raw * self.evaluator_L1.evaluate(jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll)
            pt_L2 = pt_L1 * self.evaluator_L2.evaluate(jet.eta, jet.phi, pt_L1) if self.usePhiDependentJEC else self.evaluator_L2.evaluate(jet.eta, pt_L1) * pt_L1
            pt_L3 = pt_L2 * self.evaluator_L3.evaluate(jet.eta, pt_L2)
            pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(float(event.run), jet.eta, pt_L3) if self.useRunDependentJEC else self.evaluator_L2L3.evaluate(jet.eta, pt_L3) * pt_L3

            JEC = pt_JEC / pt_raw
            mass_JEC = mass_raw * JEC

            if self.is_mc:
                JER = self.evaluator_JER.evaluate(jet.eta, pt_JEC, event.Rho_fixedGridRhoFastjetAll)
                JERsf = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "nom")
                JERsf_up = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "up")
                JERsf_dn = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "down")

                # gen match
                pt_gen = -1
                for gpt, geta, gphi in zip(gen_jets_pt, gen_jets_eta, gen_jets_phi):
                    if np.sqrt((jet.eta - geta)**2 + (jet.phi - gphi)**2) < 0.2:
                        pt_gen = gpt
                        break

                if self.use_json_smear:
                    JERsmear = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf)
                    JERsmear_up = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_up)
                    JERsmear_dn = self.evaluator_JERsmear.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_dn)
                else:
                    rand_val = float(rand_gauss(np.array([jet.pt]))[0])
                    jersf_array = np.array([[JERsf, JERsf_up, JERsf_dn]], dtype=np.float32)
                    JERsmear = float(jer_smear(pt_gen, pt_JEC, jet.eta, JER, rand_val, jersf_array, 0, self.forceStochastic))
                    JERsmear_up = float(jer_smear(pt_gen, pt_JEC, jet.eta, JER, rand_val, jersf_array, 1, self.forceStochastic))
                    JERsmear_dn = float(jer_smear(pt_gen, pt_JEC, jet.eta, JER, rand_val, jersf_array, 2, self.forceStochastic))

                # horns issue: force scaling method if no gen match and 2.5<|eta|<3.0
                if pt_gen < 0 and (2.5 < abs(jet.eta) < 3.0):
                    JERsmear_nominal = 1.0
                    JERsmear_up = JERsmear_up / JERsmear
                    JERsmear_dn = JERsmear_dn / JERsmear
                else:
                    JERsmear_nominal = JERsmear

                pt_corr.append(pt_JEC * JERsmear_nominal)
                mass_corr.append(mass_JEC * JERsmear_nominal)
                pt_smear_up.append(pt_JEC * JERsmear_up)
                pt_smear_dn.append(pt_JEC * JERsmear_dn)
                mass_smear_up.append(mass_JEC * JERsmear_up)
                mass_smear_dn.append(mass_JEC * JERsmear_dn)

                # JES total
                JESunc = self.evaluator_JES.evaluate(jet.eta, pt_JEC)
                pt_scale_up.append(pt_JEC * (1 + JESunc))
                pt_scale_dn.append(pt_JEC * (1 - JESunc))
                mass_scale_up.append(mass_JEC * (1 + JESunc))
                mass_scale_dn.append(mass_JEC * (1 - JESunc))

                # regrouped sources
                for src in self.jes_sources:
                    unc = self.evaluator_JERC[src].evaluate(jet.eta, pt_JEC)
                    name = src.split("MC_")[1].replace("_AK4PFPuppi", "")
                    pt_sources_up.setdefault(name, []).append(pt_JEC * (1 + unc))
                    pt_sources_dn.setdefault(name, []).append(pt_JEC * (1 - unc))
                    mass_sources_up.setdefault(name, []).append(mass_JEC * (1 + unc))
                    mass_sources_dn.setdefault(name, []).append(mass_JEC * (1 - unc))
            else:
                pt_corr.append(pt_JEC)
                mass_corr.append(mass_JEC)

            pt_uncorr.append(pt_raw)
            mass_uncorr.append(mass_raw)

        # fill branches
        if self.overwritePt:
            self.out.fillBranch("FatJet_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("FatJet_uncorrected_mass", mass_uncorr)

        self.out.fillBranch("FatJet_pt_nom", pt_corr)
        self.out.fillBranch("FatJet_mass_nom", mass_corr)

        if self.is_mc:
            self.out.fillBranch("FatJet_pt_jerUp", pt_smear_up)
            self.out.fillBranch("FatJet_pt_jerDown", pt_smear_dn)
            self.out.fillBranch("FatJet_mass_jerUp", mass_smear_up)
            self.out.fillBranch("FatJet_mass_jerDown", mass_smear_dn)
            self.out.fillBranch("FatJet_pt_jesTotalUp", pt_scale_up)
            self.out.fillBranch("FatJet_pt_jesTotalDown", pt_scale_dn)
            self.out.fillBranch("FatJet_mass_jesTotalUp", mass_scale_up)
            self.out.fillBranch("FatJet_mass_jesTotalDown", mass_scale_dn)

            for name in pt_sources_up:
                self.out.fillBranch(f"FatJet_pt_jes{name}Up", pt_sources_up[name])
                self.out.fillBranch(f"FatJet_pt_jes{name}Down", pt_sources_dn[name])
                self.out.fillBranch(f"FatJet_mass_jes{name}Up", mass_sources_up[name])
                self.out.fillBranch(f"FatJet_mass_jes{name}Down", mass_sources_dn[name])

        return True
