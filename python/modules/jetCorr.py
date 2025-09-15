"""Add branches for jet energy scale (JES) and resolution (JER) corrections."""

###
# Useful CMStalk thread on how implementing jet corrections: https://cms-talk.web.cern.ch/t/jes-for-2022-re-reco-cde-and-prompt-fg/32873/3
# Minimal demo provided from JME POG: https://github.com/cms-jet/JECDatabase/blob/master/scripts/JERC2JSON/minimalDemo.py
# TODO:
# - Implementing JES uncertainties. TBD which scheme we want to implement
###

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import numpy as np
import correctionlib

class jetJERC(Module):
    def __init__(self,
                json_JERC,             # Path to JERC JSON
                json_JERsmear=None,    # Path to optional JER smearing JSON
                L1Key=None,            # Key for L1 correction
                L2Key=None,            # Key for L2 correction
                L3Key=None,            # Key for L3 correction
                L2L3Key=None,          # Key for residual correction
                scaleTotalKey=None,    # JES total uncertainty key
                smearKey=None,         # Smear formula key (optional, for future)
                JERKey=None,           # JER resolution key
                JERsfKey=None,         # JER scale factor key
                overwritePt=False,     # Whether to overwrite Jet_pt
                usePhiDependentJEC=False,  # Switch for phi-dependent JEC
                useRunDependentJEC=False,   # Switch for run-dependent JEC
                isMC=False):
        """Correct jets following recommendations of JME POG.
        Parameters:
            json_JERC: full path of json file with JERC corrections
            json_JERsmear: full path of json file with smearing terms for JER
            L1Key: key for L1 corrections
            L2Key: key for L2 corrections
            L3Key: key for L3 corrections
            L2L3Key: key for residual corrections
            scaleTotalKey : key for JES uncertainties
            smearKey: key for smearing formula (None for Data)
            JERKey: key for JER (None for Data)
            JERsfKey: key for JER scale factor (None for Data)
            overwritePt: replace value in the pt branch, and store the old one as "uncorrected_pt"
        """
        self.overwritePt = overwritePt  # Save flag for overwriting pt
        self.usePhiDependentJEC = usePhiDependentJEC  # Save flag for phi-dependent correction
        self.useRunDependentJEC = useRunDependentJEC  # Save flag for run-dependent correction

        # Load JERC JSON
        self.evaluator_JERC = correctionlib.CorrectionSet.from_file(json_JERC)

        # JES evaluator (total uncertainty only)
        self.evaluator_JES = self.evaluator_JERC[scaleTotalKey] if scaleTotalKey else None
        self.jes_total = scaleTotalKey  # Save JES total key

        # Collect all JES per-source keys (filtering out total)
        self.jes_sources = [
            k for k in self.evaluator_JERC.keys()
            if k.startswith("Summer24Prompt24_V1_MC_")
            and k.endswith("_AK4PFPuppi")
            and not k.endswith("_Total_AK4PFPuppi")
            and "L1FastJet" not in k
            and "L2Relative" not in k
            and "L2L3Residual" not in k
        ]
        
        # JER evaluators (resolution + scale factors)
        self.evaluator_JER   = self.evaluator_JERC[JERKey]   if JERKey   else None
        self.evaluator_JERsf = self.evaluator_JERC[JERsfKey] if JERsfKey else None
        
        # Smear evaluator (optional JSON, usually not provided yet)
        self.evaluator_JERsmear = None
        self.evaluator_SmearFormula = None

        if json_JERsmear:
            self.evaluator_JERsmear = correctionlib.CorrectionSet.from_file(json_JERsmear)
            if smearKey and smearKey in self.evaluator_JERsmear.keys():
                self.evaluator_SmearFormula = self.evaluator_JERsmear[smearKey]

        # Flags to track availability
        self.has_jes = bool(self.evaluator_JES)
        self.has_jer = bool(self.evaluator_JER and self.evaluator_JERsf)
        self.has_smear = bool(self.evaluator_SmearFormula)
        self.is_mc = isMC  
        
        ## Starting from Run 3, the use of PUPPI jets eliminates the need for L1 corrections. To maintain compatibility with Run 2 scripts, a dummy file is provided.
        self.evaluator_L1 = self.evaluator_JERC[L1Key]

        ## For the next step, there is a mismatch between the terminology in the twiki and the tags in correctionlib
        ## MC-truth = L2 + L3
        self.evaluator_L2 = self.evaluator_JERC[L2Key]
        self.evaluator_L3 = self.evaluator_JERC[L3Key]
        
        ## L2L3residuals should be applied only to data
        ## A dummy file with all entries equal to 1 is provided for MC to have a common script for both data and MC
        self.evaluator_L2L3 = self.evaluator_JERC[L2L3Key]

 
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        # Nominal corrected branches (always present)
        self.out.branch("Jet_pt_nom", "F", lenVar="nJet", title="Nominal corrected pt")
        self.out.branch("Jet_mass_nom", "F", lenVar="nJet", title="Nominal corrected mass", limitedPrecision=12)

        self.out.branch("Jet_JECfactor", "F", lenVar="nJet", title="JEC correction factor (pt_corr/pt_raw)")

        if self.overwritePt:
        # If overwriting pt/mass, keep raw values too
            self.out.branch("Jet_uncorrected_pt", "F", lenVar="nJet", title="original (uncorrected) pT")
            self.out.branch("Jet_uncorrected_mass", "F", lenVar="nJet", title="original (uncorrected) mass", limitedPrecision=12)
       


        
        # If MC: add JEC/JER/JES branches
        if self.is_mc:
            
            # JER scale factor branches
            if self.has_jer:
                self.out.branch("Jet_JERsf_nom", "F", lenVar="nJet", title="Nominal JER scale factor")
                self.out.branch("Jet_JERsf_up", "F", lenVar="nJet", title="JER scale factor up")
                self.out.branch("Jet_JERsf_down", "F", lenVar="nJet", title="JER scale factor down")
                self.out.branch("Jet_pt_jerUp", "F", lenVar="nJet", title="JER up variation")
                self.out.branch("Jet_pt_jerDown", "F", lenVar="nJet", title="JER down variation")
                self.out.branch("Jet_mass_jerUp", "F", lenVar="nJet", title="JER up variation")
                self.out.branch("Jet_mass_jerDown", "F", lenVar="nJet", title="JER down variation")

            # JER smear branches (only if smearing JSON exists)
            if self.has_smear:
                self.out.branch("Jet_JERsmear_nom", "F", lenVar="nJet", title="Nominal JER smear factor")
                self.out.branch("Jet_JERsmear_up", "F", lenVar="nJet", title="JER smear factor up")
                self.out.branch("Jet_JERsmear_down", "F", lenVar="nJet", title="JER smear factor down")


            # JES uncertainties (total and per source)
            if self.has_jes:
                # Global JES
                self.out.branch("Jet_pt_jesTotalUp", "F", lenVar="nJet", title="JES total up variation")
                self.out.branch("Jet_pt_jesTotalDown", "F", lenVar="nJet", title="JES total down variation")
                self.out.branch("Jet_mass_jesTotalUp", "F", lenVar="nJet", title="JES total up variation")
                self.out.branch("Jet_mass_jesTotalDown", "F", lenVar="nJet", title="JES total down variation")
                # Per-source JES
                for src in self.jes_sources:
                    name = src.split("MC_")[1].replace("_AK4PFPuppi","")  # e.g. "AbsoluteScale"
                    self.out.branch(f"Jet_pt_jes{name}Up", "F", lenVar="nJet", title=f"JES {name} up variation")
                    self.out.branch(f"Jet_pt_jes{name}Down", "F", lenVar="nJet", title=f"JES {name} down variation")
                    self.out.branch(f"Jet_mass_jes{name}Up", "F", lenVar="nJet", title=f"JES {name} up variation")
                    self.out.branch(f"Jet_mass_jes{name}Down", "F", lenVar="nJet", title=f"JES {name} down variation")

         
    def fixPhi(self, phi)   :# Utility to wrap phi into [-pi, pi]
        if phi > np.pi:
            phi -= 2*np.pi
        elif phi < -np.pi:
            phi += 2*np.pi
        return phi

    def analyze(self, event):
        jets = Collection(event, "Jet")
        if self.is_mc: ## genJet info is necessary for JER
            gen_jets = Collection(event, "GenJet")
            gen_jets_pt = np.array([gen_jet.pt for gen_jet in gen_jets])
            gen_jets_eta = np.array([gen_jet.eta for gen_jet in gen_jets])
            gen_jets_phi = np.array([gen_jet.phi for gen_jet in gen_jets])

        pt_corr = []
        pt_uncorr = []
        mass_corr = []
        mass_uncorr = []
        pt_smear_up = []
        pt_smear_dn = []
        mass_smear_up = []
        mass_smear_dn = []
        pt_scale_up = []
        pt_scale_dn = []
        mass_scale_up = []
        mass_scale_dn = []
        jec_factors = []
        jer_sf_nom_vals, jer_sf_up_vals, jer_sf_dn_vals = [], [], []
        jer_smear_nom_vals, jer_smear_up_vals, jer_smear_dn_vals = [], [], []
        
        pt_sources_up, pt_sources_dn = {}, {}
        mass_sources_up, mass_sources_dn = {}, {}

        for jet in jets:
            #### JEC ####
            ## To be applied to both data and MC
            ## Jet in NanoAOD are already corrected
            ## The correction should be removed and the latest one available should be applied
            pt_raw = jet.pt * (1 - jet.rawFactor)
            mass_raw = jet.mass * (1 - jet.rawFactor)

            # Apply step-wise JEC
            ## The three steps of JEC corrections are provided separately
            # --- L1FastJet ---
            pt_L1 = pt_raw * self.evaluator_L1.evaluate(
                jet.area, jet.eta, pt_raw, event.Rho_fixedGridRhoFastjetAll
            )

            # --- L2Relative ---
            pt_L2 = pt_L1 * self.evaluator_L2.evaluate(
                jet.eta, jet.phi, pt_L1
            )

            # --- L3Absolute ---
            pt_L3 = pt_L2 * self.evaluator_L3.evaluate(
                jet.eta, pt_L2
            )

            # --- Residual (Data only, or dummy=1.0 for MC) ---
            if self.useRunDependentJEC:
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(
                    float(event.run), jet.eta, pt_L3
                )
            else:
                pt_JEC = pt_L3 * self.evaluator_L2L3.evaluate(
                    jet.eta, pt_L3
                )

            # Compute JEC factor and apply to mass
            JEC = pt_JEC / pt_raw
            mass_JEC = mass_raw * JEC
            jec_factors.append(JEC)


            if self.is_mc:
                #### JER ####
                ## Hybrid method is implemented [https://cms-jerc.web.cern.ch/JER/#smearing-procedures]
                ## except in the eta region 2.5 < |eta| < 3.0 where the Scaling method is implemented
                ## Scaling method was introduced following JME recommendations to handle the 'Horns issue'
                ## Ref: https://gitlab.cern.ch/cms-jetmet/coordination/coordination/-/issues/113
                ## TODO: This should be removed once a JSON-level fix is implemented.
                
               
                # --- JES uncertainties ---
                # --- Global JES ---
                JESuncert_total = self.evaluator_JERC[self.jes_total].evaluate(jet.eta, pt_JEC)
                pt_JES_total_up = pt_JEC * (1 + JESuncert_total)
                pt_JES_total_dn = pt_JEC * (1 - JESuncert_total)
                pt_scale_up.append(pt_JES_total_up)
                pt_scale_dn.append(pt_JES_total_dn)
                
                mass_JES_total_up = mass_JEC * (1 + JESuncert_total)
                mass_JES_total_dn = mass_JEC * (1 - JESuncert_total)
                mass_scale_up.append(mass_JES_total_up)
                mass_scale_dn.append(mass_JES_total_dn)

                # --- Per-source JES ---
                for src in self.jes_sources:
                    unc = self.evaluator_JERC[src].evaluate(jet.eta, pt_JEC)
                    name = src.split("MC_")[1].replace("_AK4PFPuppi","")
                    # Append to dictionaries so you can fill later
                    if name not in pt_sources_up:
                        pt_sources_up[name], pt_sources_dn[name] = [], []
                        mass_sources_up[name], mass_sources_dn[name] = [], []
                    pt_sources_up[name].append(pt_JEC * (1 + unc))
                    pt_sources_dn[name].append(pt_JEC * (1 - unc))
                    mass_sources_up[name].append(mass_JEC * (1 + unc))
                    mass_sources_dn[name].append(mass_JEC * (1 - unc))

                # --- JER corrections ---
                if self.has_jer:
                    JER = self.evaluator_JER.evaluate(jet.eta, pt_JEC, event.Rho_fixedGridRhoFastjetAll)
                    ## GenMatching with genJet
                    delta_eta = jet.eta - gen_jets_eta
                    fixPhi = np.vectorize(self.fixPhi, otypes=[float])
                    delta_phi = fixPhi(jet.phi - gen_jets_phi)
                    pt_gen = np.where((np.abs(pt_JEC - gen_jets_pt) < 3 * pt_JEC * JER) & (np.sqrt(delta_eta**2 + delta_phi**2)<0.2), gen_jets_pt, -1.0)
                    pt_gen = pt_gen[pt_gen > 0][0] if np.any(pt_gen > 0) else -1. ## If no gen-matching, simply -1
                    
                    # Get JER scale factors
                    JERsf = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "nom")
                    JERsf_up = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "up")
                    JERsf_dn = self.evaluator_JERsf.evaluate(jet.eta, jet.pt, "down")
                    jer_sf_nom_vals.append(JERsf)
                    jer_sf_up_vals.append(JERsf_up)
                    jer_sf_dn_vals.append(JERsf_dn)

                    # Apply smearing if JSON exists
                    if self.has_smear:
                        JERsmear = self.evaluator_SmearFormula.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf)
                        
                        JERsmear_up = self.evaluator_SmearFormula.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_up)
                        
                        JERsmear_dn = self.evaluator_SmearFormula.evaluate(pt_JEC, jet.eta, pt_gen, event.Rho_fixedGridRhoFastjetAll, event.event, JER, JERsf_dn)
                        
                        jer_smear_nom_vals.append(JERsmear)
                        jer_smear_up_vals.append(JERsmear_up)
                        jer_smear_dn_vals.append(JERsmear_dn)
                        

                        if pt_gen < 0 and (2.5 < abs(jet.eta) < 3):
                            JERsmear_nominal = 1.0
                            JERsmear_up = JERsmear_up / JERsmear
                            JERsmear_dn = JERsmear_dn / JERsmear
                        else:
                            JERsmear_nominal = JERsmear

                    else:
                        # No smearing JSON: fallback
                        JERsmear_nominal = 1.0
                        JERsmear_up = 1.0
                        JERsmear_dn = 1.0
                        

                    pt_corr.append(pt_JEC * JERsmear_nominal)
                    mass_corr.append(mass_JEC * JERsmear_nominal)
                    pt_smear_up.append(pt_JEC * JERsmear_up)
                    pt_smear_dn.append(pt_JEC * JERsmear_dn)
                    mass_smear_up.append(mass_JEC * JERsmear_up)
                    mass_smear_dn.append(mass_JEC * JERsmear_dn)
                else:
                    # No JER available → just JEC
                    pt_corr.append(pt_JEC)
                    mass_corr.append(mass_JEC)
                
                # Always keep uncorrected values
                pt_uncorr.append(pt_raw)
                mass_uncorr.append(mass_raw)
 
            else:
                # Data: just apply JEC
                pt_corr.append(pt_JEC)
                pt_uncorr.append(pt_raw)
                mass_corr.append(mass_JEC)
                mass_uncorr.append(mass_raw)

        if self.overwritePt:
            self.out.fillBranch("Jet_uncorrected_pt", pt_uncorr)
            self.out.fillBranch("Jet_uncorrected_mass", mass_uncorr)
        
        self.out.fillBranch("Jet_pt_nom", pt_corr)
        self.out.fillBranch("Jet_mass_nom", mass_corr)
        self.out.fillBranch("Jet_JECfactor", jec_factors)
    
        if self.is_mc and self.has_jes:
            # Fill JES total variations

            # Global JES
            self.out.fillBranch("Jet_pt_jesTotalUp", pt_scale_up)
            self.out.fillBranch("Jet_pt_jesTotalDown", pt_scale_dn)
            self.out.fillBranch("Jet_mass_jesTotalUp", mass_scale_up)
            self.out.fillBranch("Jet_mass_jesTotalDown", mass_scale_dn)


            # Fill JER scale factors
            if self.has_jer:
                self.out.fillBranch("Jet_JERsf_nom", jer_sf_nom_vals)
                self.out.fillBranch("Jet_JERsf_up", jer_sf_up_vals)
                self.out.fillBranch("Jet_JERsf_down", jer_sf_dn_vals)
            
            # Fill JER smearing if available
            if self.has_smear:
                self.out.fillBranch("Jet_JERsmear_nom", jer_smear_nom_vals)
                self.out.fillBranch("Jet_JERsmear_up", jer_smear_up_vals)
                self.out.fillBranch("Jet_JERsmear_down", jer_smear_dn_vals)

            # Fill JES per-source uncertainties
            for name in pt_sources_up:
                self.out.fillBranch(f"Jet_pt_jes{name}Up", pt_sources_up[name])
                self.out.fillBranch(f"Jet_pt_jes{name}Down", pt_sources_dn[name])
                self.out.fillBranch(f"Jet_mass_jes{name}Up", mass_sources_up[name])
                self.out.fillBranch(f"Jet_mass_jes{name}Down", mass_sources_dn[name])

            # Fill JER pt/mass variations
            if self.has_jer:    
                # JER variations
                self.out.fillBranch("Jet_pt_jerUp", pt_smear_up)
                self.out.fillBranch("Jet_pt_jerDown", pt_smear_dn)
                self.out.fillBranch("Jet_mass_jerUp", mass_smear_up)
                self.out.fillBranch("Jet_mass_jerDown", mass_smear_dn)

        return True
