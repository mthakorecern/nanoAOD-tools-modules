#!/usr/bin/env python3

import os, math, json, gzip
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

import numpy as np
import correctionlib._core as correction
from typing import Tuple

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection

applyOnMET = True

def phi_mpi_pi(dphi: float) -> float:
    """Wrap an angle difference into the range [-pi, pi].
    Used to compute delta phi consistently."""
    return math.atan2(math.sin(dphi), math.cos(dphi))

def deltaR(eta1: float, phi1: float, eta2: float, phi2: float) -> float:
    """Compute deltaR separation between two objects.
    Uses eta difference and phi difference wrapped into [-pi, pi]."""
    dEta = eta1 - eta2
    dPhi = phi_mpi_pi(phi1 - phi2)
    return math.sqrt(dEta*dEta + dPhi*dPhi)

def representative_run_number(year: str, run: int) -> float:
    """Return a representative run number for given year.
    Used for evaluating residual corrections on Data."""
    if year == "2023Pre": return 367080.0
    if year == "2023Post": return 369803.0
    if year == "2024": return 379412.0
    return float(run)

def normal_from_seeds(seedA: int, seedB: int) -> float:
    """Generate a reproducible Gaussian random number.
    Uses seeds from run, lumi, event, and jet index."""
    UINT64_MAX = np.uint64(0xFFFFFFFFFFFFFFFF)
    a = (np.uint64(seedA) ^ np.uint64(0x9E3779B97F4A7C15)) * np.uint64(0xBF58476D1CE4E5B9) & UINT64_MAX
    b = (np.uint64(seedB) ^ np.uint64(0x94D049BB133111EB)) * np.uint64(0xD2B74407B1CE6E93) & UINT64_MAX
    u1 = (a.astype(np.float64) + 1.0) / (UINT64_MAX.astype(np.float64) + 2.0)
    u2 = (b.astype(np.float64) + 1.0) / (UINT64_MAX.astype(np.float64) + 2.0)
    return math.sqrt(-2.0*math.log(u1)) * math.cos(2.0*math.pi*u2)

def build_jet_seeds(run: int, lumi: int, event: int, jet_index: int) -> Tuple[int,int]:
    """Build integer seeds for random smearing of jets.
    Seeds are derived from run, lumi, event, and jet index."""
    seedA = (np.uint64(run) << np.uint64(1)) ^ np.uint64(lumi) ^ np.uint64(jet_index)
    seedB = (np.uint64(event) << np.uint64(1)) ^ np.uint64(jet_index * 1315423911)
    return int(seedA), int(seedB)

class CorrectionRefs:
    """Helper class to load and store correction references.
    Provides access to L1, L2, L3, residuals, JER resolution and scale factors."""
    def __init__(self,
                tagNameL1FastJet,
                tagNameL2Relative,
                tagNameL3Absolute,
                tagNameL2L3Residual,
                tagNamePtResolution,
                tagNameJerScaleFactor,
                jercjson: str,
                isData: bool):
        
        self.cs = correction.CorrectionSet.from_file(jercjson)
        
        #JES 
        self.cL1 = self.cs[tagNameL1FastJet]
        self.cL2 = self.cs[tagNameL2Relative]
        self.cL3 = self.cs[tagNameL3Absolute]
        
        # Since for Monte Carlo: L1 + MC-truth &&  Data: L1 + MC-truth + L2L3Residuals
        if isData and tagNameL2L3Residual:
            self.cL2L3Res = self.cs[tagNameL2L3Residual]
        else:
            self.cL2L3Res = None        
        
        #JER based
        self.cReso = self.cs[tagNamePtResolution] if (not isData and tagNamePtResolution) else None
        self.cJerSF = self.cs[tagNameJerScaleFactor] if (not isData and tagNameJerScaleFactor) else None


## Mapping: branch -> JSON key
BRANCH_TO_JSON = {
    "jer": {
        "PtResolution": "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
        "ScaleFactor": "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
    },
    ## Single Total - 1
    "jesTotal": "Summer24Prompt24_V1_MC_Total_AK4PFPuppi",
    
    ## Regrouped - 11 Sources
    "jesAbsolute": "Summer24Prompt24_V1_MC_Regrouped_Absolute_AK4PFPuppi",
    "jesAbsolute_2024": "Summer24Prompt24_V1_MC_Regrouped_Absolute_2024_AK4PFPuppi",
    "jesBBEC1": "Summer24Prompt24_V1_MC_Regrouped_BBEC1_AK4PFPuppi",
    "jesBBEC1_2024": "Summer24Prompt24_V1_MC_Regrouped_BBEC1_2024_AK4PFPuppi",
    "jesEC2": "Summer24Prompt24_V1_MC_Regrouped_EC2_AK4PFPuppi",
    "jesEC2_2024": "Summer24Prompt24_V1_MC_Regrouped_EC2_2024_AK4PFPuppi",
    "jesFlavorQCD": "Summer24Prompt24_V1_MC_Regrouped_FlavorQCD_AK4PFPuppi",
    "jesHF": "Summer24Prompt24_V1_MC_Regrouped_HF_AK4PFPuppi",
    "jesHF_2024": "Summer24Prompt24_V1_MC_Regrouped_HF_2024_AK4PFPuppi",
    "jesRelativeBal": "Summer24Prompt24_V1_MC_Regrouped_RelativeBal_AK4PFPuppi",
    "jesRelativeSample_2024": "Summer24Prompt24_V1_MC_Regrouped_RelativeSample_2024_AK4PFPuppi",
    
    ##Scanning tags for Data and Monte Carlo
    "L1FastJet": {
        "MC": "Summer24Prompt24_V1_MC_L1FastJet_AK4PFPuppi",
        "DATA": "Summer24Prompt24_V1_DATA_L1FastJet_AK4PFPuppi",
    },
    "L2Relative": {
        "MC": "Summer24Prompt24_V1_MC_L2Relative_AK4PFPuppi",
        "DATA": "Summer24Prompt24_V1_DATA_L2Relative_AK4PFPuppi",
    },
    "L3Absolute": {
        "MC": "Summer24Prompt24_V1_MC_L3Absolute_AK4PFPuppi",
        "DATA": "Summer24Prompt24_V1_DATA_L3Absolute_AK4PFPuppi",
    },
    "L2L3Residual": {
        "MC": "Summer24Prompt24_V1_MC_L2L3Residual_AK4PFPuppi",
        "DATA": "Summer24Prompt24_V1_DATA_L2L3Residual_AK4PFPuppi",
    },
}

def get_event_collections(event, isData: bool):
    return {
        "run": event.run,
        "lumi": event.luminosityBlock,
        "evt": event.event,
        "rho": event.Rho_fixedGridRhoFastjetAll,
        "met": event.PuppiMET_pt,
        "metphi": event.PuppiMET_phi,
        "rawPuppiMET": getattr(event, "RawPuppiMET_pt", -1),
        "rawPuppiMETphi": getattr(event, "RawPuppiMET_phi", -99),
        "jets": Collection(event, "Jet"),
        "fatjets": Collection(event, "FatJet"),
        "genjets": Collection(event, "GenJet") if not isData else [],
        "genjetsAK8": Collection(event, "GenJetAK8") if not isData else [],
    }

def branch_name_from_sys(sys: str, year_unc: str) -> str:
    """Return a branch suffix string for systematics.
    Uses nominal if sys is empty, otherwise returns sys."""
    if not sys: return "nom"
    return sys

class ApplyJercAll(Module):
    def __init__(self, year, isData, jercjson, era=None, year_unc="2024"):
        self.year = year
        self.isData = bool(isData)
        self.isMC   = not self.isData
        self.era = era
        self.year_unc = year_unc
        self.jercjson = jercjson

        print("ApplyJercAll initialized with:")
        print("  year =", self.year)
        print("  isData =", self.isData)
        print("  era =", self.era)
        print("  jercjson =", self.jercjson)

        # Correction refs (pass explicit keys)
        self.refsAK4 = CorrectionRefs(
            tagNameL1FastJet=BRANCH_TO_JSON["L1FastJet"]["DATA" if self.isData else "MC"],
            tagNameL2Relative=BRANCH_TO_JSON["L2Relative"]["DATA" if self.isData else "MC"],
            tagNameL3Absolute=BRANCH_TO_JSON["L3Absolute"]["DATA" if self.isData else "MC"],
            tagNameL2L3Residual=BRANCH_TO_JSON["L2L3Residual"]["DATA"] if self.isData else None,
            tagNamePtResolution=BRANCH_TO_JSON["jer"]["PtResolution"] if self.isMC else None,
            tagNameJerScaleFactor=BRANCH_TO_JSON["jer"]["ScaleFactor"] if self.isMC else None,
            jercjson=jercjson,
            isData=self.isData,
        )
        self.refsAK8 = self.refsAK4  # reusing AK4 keys for AK8
        self.csAK4 = self.refsAK4.cs

        
        """Build the list of systematic variations for jet corrections.
        Includes nominal, all JES (Up/Down) and JER (Up/Down) variations."""
        self.jet_pt_systematics = [""]

        if self.isMC:
            for base, full in BRANCH_TO_JSON.items():
                if base.startswith("jes") and isinstance(full, str):
                    for var in ("Up", "Down"):
                        self.jet_pt_systematics.append(f"{base}{var}")
            for var in ("Up", "Down"):
                self.jet_pt_systematics.append(f"jer{var}")

        print("Systematic variations to be produced:")
        for sys in self.jet_pt_systematics:
            print("  ", branch_name_from_sys(sys, self.year_unc))


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        print("Starting new file:", inputFile.GetName())
        
        for sys in self.jet_pt_systematics:
            suffix = branch_name_from_sys(sys, self.year_unc)
            #print("Creating branches with suffix:", suffix)
            self.out.branch(f"Jet_pt_{suffix}", "F", lenVar="nJet")
            self.out.branch(f"Jet_mass_{suffix}", "F", lenVar="nJet")

            self.out.branch(f"FatJet_pt_{suffix}", "F", lenVar="nFatJet")
            self.out.branch(f"FatJet_mass_{suffix}", "F", lenVar="nFatJet")
            self.out.branch(f"FatJet_msoftdrop_{suffix}", "F", lenVar="nFatJet")
            if applyOnMET:
                self.out.branch(f"PuppiMET_pt_{suffix}", "F")
                self.out.branch(f"PuppiMET_phi_{suffix}", "F")
        
        self.out.branch("FatJet_globalParT3Xbb_mass", "F", lenVar="nFatJet")
        #print("Created branch FatJet_globalParT3Xbb_mass")


    def analyze(self, event):
        ev = get_event_collections(event, self.isData)
        run,lumi,evt = ev["run"], ev["lumi"], ev["evt"]
        rho = ev["rho"]
        jets,fatjets = ev["jets"], ev["fatjets"]
        genjets, genjetsAK8 = ev["genjets"], ev["genjetsAK8"]

        #print("Processing event:", evt, "run:", run, "lumi:", lumi)


        # GlobalParT3 masses
        """Compute GlobalParT3 corrected fatjet masses.
        Use globalParT3_massCorrX2p if available, otherwise fill with -999."""
        gp3_masses = []
        
        for fj in fatjets:
            if hasattr(fj,"globalParT3_massCorrX2p"):
                raw_mass = fj.mass * (1.0 - fj.rawFactor)
                gp3_masses.append(fj.globalParT3_massCorrX2p * raw_mass)
            else: gp3_masses.append(-999.0)
        
        self.out.fillBranch("FatJet_globalParT3Xbb_mass", gp3_masses)
        #print("Filled FatJet_globalParT3Xbb_mass with", len(gp3_masses), "entries")

        
        
        # Variations loop
        for sys in self.jet_pt_systematics:
            suffix = branch_name_from_sys(sys, self.year_unc)
            #print("Applying corrections for systematics:", suffix)
            jet_corr_pts, jet_corr_mass = [], []
            fat_corr_pts, fat_corr_mass, fat_corr_msoft = [], [], []

            # --- AK4 jets ---
            for j in jets:
                """
                AK4 jet correction chain:
                  - Start from raw pt (pt_raw = pt * (1 - rawFactor))
                  - Apply L1FastJet
                  - Apply L2Relative
                  - Apply:
                      * Data -> L2L3Residual
                      * MC   -> L3Absolute
                  - Then apply JES systematics (Up/Down variations)
                  - Then apply JER systematics (only in MC, using genJet match or smearing)
                """
                pt_raw = j.pt * (1.0 - j.rawFactor)
                mass_raw = j.mass * (1.0 - j.rawFactor)

                
                c1 = self.refsAK4.cL1.evaluate(j.area,j.eta,pt_raw,rho)
                pt = pt_raw * c1
                
                c2 = self.refsAK4.cL2.evaluate(j.eta,j.phi,pt)
                pt *= c2

                if self.isData and self.refsAK4.cL2L3Res:
                    run_res = representative_run_number(self.year,run)
                    cR = self.refsAK4.cL2L3Res.evaluate(run_res,j.eta,pt)
                    pt *= cR
                else:
                    # Apply L3Absolute for MC
                    pt *= self.refsAK4.cL3.evaluate(j.eta, pt)

                pt_corr = pt

                # JES syst
                if sys.startswith("jes") and sys.replace("Up","").replace("Down","") in self.csAK4:
                    var = "Up" if sys.endswith("Up") else "Down"
                    key = sys.replace("Up","").replace("Down","")
                    #print("    Applying JES key =", key, "variation =", var)
                    
                    scale = self.csAK4[key].evaluate(j.eta,pt_corr)
                    pt_corr *= (1 + scale) if var=="Up" else (1 - scale)

                # JER syst
                if sys.startswith("jer") and self.isMC:
                    jer_var = "up" if sys.endswith("Up") else "down"
                    reso = self.refsAK4.cReso.evaluate(j.eta,pt_corr,rho)
                    sf = self.refsAK4.cJerSF.evaluate(j.eta,pt_corr,jer_var)
                
                    gpt,is_match=0.0,False
                
                    if hasattr(j,"genJetIdx") and j.genJetIdx>=0 and j.genJetIdx<len(genjets):
                        g = genjets[j.genJetIdx]
                        gpt = g.pt
                        if deltaR(j.eta, j.phi, g.eta, g.phi) < 0.2: is_match=True
                
                    if is_match:
                        pt_corr = max(0.0, 1.0 + (sf - 1.0) * (pt_corr - gpt)/pt_corr) *pt_corr
                
                    else:
                        seedA,seedB = build_jet_seeds(run,lumi,evt,j._index)
                        z = normal_from_seeds(seedA,seedB)
                        sigma = math.sqrt(max(sf*sf-1.0,0.0)) * reso
                        pt_corr *= max(0.0, 1.0 + z*sigma)

                mass_corr = mass_raw * (pt_corr / pt_raw) if pt_raw > 0 else mass_raw
                jet_corr_pts.append(pt_corr)
                jet_corr_mass.append(mass_corr)


            # --- AK8 jets ---
            for fj in fatjets:
                """
                AK8 jet correction chain:
                  - Start from raw pt (pt_raw = pt * (1 - rawFactor))
                  - Apply L1FastJet
                  - Apply L2Relative
                  - Apply:
                      * Data -> L2L3Residual
                      * MC   -> L3Absolute
                  - Then apply JES systematics (Up/Down variations)
                  - Then apply JER systematics (only in MC, using genJetAK8 match or smearing)
                  - Also propagate correction to softdrop mass
                """
                pt_raw=fj.pt * (1.0-fj.rawFactor)
                mass_raw = fj.mass * (1.0 - fj.rawFactor)
                msoft_raw = fj.msoftdrop * (1.0 - fj.rawFactor)


                
                c1 = self.refsAK4.cL1.evaluate(fj.area,fj.eta,pt_raw,rho)
                pt = pt_raw * c1
                
                c2 = self.refsAK4.cL2.evaluate(fj.eta,fj.phi,pt)
                pt *= c2

                if self.isData and self.refsAK4.cL2L3Res:
                    run_res = representative_run_number(self.year,run)
                    cR = self.refsAK4.cL2L3Res.evaluate(run_res,fj.eta,pt)
                    pt *= cR
                
                else:
                    # Apply L3Absolute for MC
                    pt *= self.refsAK4.cL3.evaluate(fj.eta, pt)

                pt_corr = pt

                if sys.startswith("jes") and sys.replace("Up","").replace("Down","") in self.csAK4:
                    var = "Up" if sys.endswith("Up") else "Down"
                    key = sys.replace("Up","").replace("Down","")
                    #print("    Applying JES key =", key, "variation =", var)
                    scale = self.csAK4[key].evaluate(fj.eta,pt_corr)
                    pt_corr *= (1 + scale) if var == "Up" else ( 1 - scale)

                if sys.startswith("jer") and self.isMC:
                    jer_var = "up" if sys.endswith("Up") else "down"
                    reso= self.refsAK4.cReso.evaluate(fj.eta,pt_corr,rho)
                    sf= self.refsAK4.cJerSF.evaluate(fj.eta,pt_corr,jer_var)
                
                    gpt, is_match = 0.0, False
                
                    if hasattr(fj,"genJetAK8Idx") and fj.genJetAK8Idx>=0 and fj.genJetAK8Idx<len(genjetsAK8):
                        g = genjetsAK8[fj.genJetAK8Idx]
                        gpt = g.pt
                        if deltaR(fj.eta,fj.phi,g.eta,g.phi)<0.4: is_match=True
                
                    if is_match:
                        pt_corr = max(0.0, 1.0 + (sf - 1.0) * (pt_corr-gpt) / pt_corr) * pt_corr
                
                    else:
                        seedA,seedB = build_jet_seeds(run,lumi,evt,fj._index)
                        z = normal_from_seeds(seedA,seedB)
                        sigma = math.sqrt(max(sf * sf - 1.0, 0.0)) * reso
                        pt_corr *= max(0.0, 1.0 + z*sigma)
                
                mass_corr  = mass_raw  * (pt_corr / pt_raw) if pt_raw > 0 else mass_raw
                msoft_corr = msoft_raw * (pt_corr / pt_raw) if pt_raw > 0 else msoft_raw

                fat_corr_pts.append(pt_corr)
                fat_corr_mass.append(mass_corr)
                fat_corr_msoft.append(msoft_corr)

            suffix = branch_name_from_sys(sys,self.year_unc)
            self.out.fillBranch(f"Jet_pt_{suffix}",jet_corr_pts)
            self.out.fillBranch(f"Jet_mass_{suffix}", jet_corr_mass)
            self.out.fillBranch(f"FatJet_pt_{suffix}",fat_corr_pts)
            self.out.fillBranch(f"FatJet_mass_{suffix}", fat_corr_mass)
            self.out.fillBranch(f"FatJet_msoftdrop_{suffix}",fat_corr_msoft)
            #print("  Filled Jet_pt_", suffix, "with", len(jet_corr_pts), "jets")
            #print("  Filled FatJet_pt_", suffix, "with", len(fat_corr_pts), "fatjets")

            if applyOnMET:

                """
                PuppiMET Type-1 correction chain:
                  - Start from raw PuppiMET (preferred in 2022+; fallback to nominal PuppiMET otherwise).
                  - Recompute MET_x, MET_y vectors.

                  For each AK4 jet:
                    * Undo rawFactor -> pt_raw.
                    * Apply L1FastJet correction.
                    * Apply muon subtraction (Jet_muonSubtrFactor).
                    * Apply L2Relative correction (eta/phi-dependent for 2023/2024).
                    * Apply:
                        - Data -> L2L3Residual (run-dependent).
                        - MC   -> L3Absolute (done earlier).
                    * Define this as nominal corrected pt (pt_corr_nom).
                    * Compare to systematic-shifted pt (pt_corr from JES/JER loop).
                    * Compute Deltapt = (pt_corr - pt_corr_nom).
                    * Update MET_x, MET_y: subtract Deltapt * cos(phi), Deltapt * sin(phi).

                  - After looping over jets, recompute:
                        PuppiMET_pt  = sqrt(met_x^2 + met_y^2)
                        PuppiMET_phi = atan2(met_y, met_x)

                  - These corrected values are stored in branches:
                        PuppiMET_pt_<sys>, PuppiMET_phi_<sys>
                  where <sys> can be 'nom' or JES/JER variation tags.
                """
                # --- Start from raw PuppiMET (preferred in 2022+, fallback otherwise)
                base_met_pt  = ev["rawPuppiMET"] if ev["rawPuppiMET"] > 0 else ev["met"]
                base_met_phi = ev["rawPuppiMETphi"] if ev["rawPuppiMETphi"] > -90 else ev["metphi"]
                met_x = base_met_pt * math.cos(base_met_phi)
                met_y = base_met_pt * math.sin(base_met_phi)

                # Loop over jets and recompute Type-1 corrections
                for j, pt_corr in zip(jets, jet_corr_pts):
                    # Raw jet pt after undoing JEC
                    pt_raw = j.pt * (1.0 - j.rawFactor)

                    # L1 correction
                    c1 = self.refsAK4.cL1.evaluate(j.area, j.eta, pt_raw, rho)
                    pt_l1 = pt_raw * c1

                    # Subtract muons (NanoAOD: Jet_muonSubtrFactor)
                    if hasattr(j, "muonSubtrFactor"):
                        pt_l1 *= (1.0 - j.muonSubtrFactor)

                    # L2Relative (sometimes phi-dependent in 2023/2024)
                    c2 = self.refsAK4.cL2.evaluate(j.eta, j.phi, pt_l1)
                    pt_corr_nom = pt_l1 * c2

                    # Residuals (only Data, year-dependent)
                    if self.isData and self.refsAK4.cL2L3Res:
                        run_res = representative_run_number(self.year, run)
                        cR = self.refsAK4.cL2L3Res.evaluate(run_res, j.eta, pt_corr_nom)
                        pt_corr_nom *= cR

                    # Now apply JES/JER systematic shift (already in pt_corr)
                    # Deltapt = new (systematic) - nominal corrected pt
                    delta_pt = pt_corr - pt_corr_nom

                    # Type-1 MET propagation
                    met_x -= delta_pt * math.cos(j.phi)
                    met_y -= delta_pt * math.sin(j.phi)

                met_pt_corr  = math.sqrt(met_x**2 + met_y**2)
                met_phi_corr = math.atan2(met_y, met_x)
                self.out.fillBranch(f"PuppiMET_pt_{suffix}", float(met_pt_corr))
                self.out.fillBranch(f"PuppiMET_phi_{suffix}", float(met_phi_corr))
                #print("  Filled PuppiMET branches with pt =", met_pt_corr, "phi =", met_phi_corr)


        return True

