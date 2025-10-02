#!/usr/bin/env python3

# NanoAODTools-based module version of JERC 

import os, math, json
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from dataclasses import dataclass
from functools import lru_cache
import gzip
from typing import Dict, List, Optional, Tuple

import numpy as np
import correctionlib._core as correction

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor



# ---------------------------
# Global toggles (same as before)
# ---------------------------
applyOnlyOnAK4 = False
applyOnlyOnAK8 = False
applyOnAK4AndAK8 = True
applyOnMET = True
isDebug = True


def printDebug(*args, indent: int = 0, enable: bool = None) -> None:
    """Conditionally print debug information with indentation."""
    if enable is None:
        enable = isDebug
    if not enable:
        return

    prefix = " " * max(indent, 0)
    message = " ".join(str(arg) for arg in args)
    print(prefix + message if prefix else message)

def debug_object(obj_name: str, obj, indent: int = 2, enable: bool = None):
    """
    Print basic kinematics of an object collection (NanoAODTools Collection).
    Example: jets, fatjets, etc.
    """
    if enable is None:
        enable = isDebug
    if not enable:
        return
    
    printDebug(f"{obj_name}: n = {len(obj)}", indent=indent)
    for i, o in enumerate(obj):
        if hasattr(o, "pt") and hasattr(o, "eta") and hasattr(o, "phi"):
            printDebug(f"  idx={i} pt={o.pt:.2f} eta={o.eta:.2f} phi={o.phi:.2f}", indent=indent+2)
        else:
            printDebug(f"  idx={i} {o}", indent=indent+2)

def debug_met(met_obj, indent: int = 2, enable: bool = None):
    """Print MET object values."""
    if enable is None:
        enable = isDebug
    if not enable:
        return

    printDebug(f"MET_pt={met_obj.pt:.2f}, phi={met_obj.phi:.2f}", indent=indent)


@lru_cache(maxsize=None)
def load_json_cached(path: str) -> dict:
    """Load and cache a JSON or JSON.GZ config file."""
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:   # "rt" = text mode
            return json.load(f)
    else:
        with open(path, "r") as f:
            return json.load(f)


class CorrCache:
    """
    Cache correctionlib correction sets (JEC/JER/JVM).
    Reuses them across all events once per year/isData/era.
    """
    _cfg = {}   # kind -> parsed tag JSON (AK4/AK8)
    _refs = {}  # (kind, year, isData, era) -> CorrectionRefs

    @classmethod
    def get_cfg(cls, kind: str, jercjson: str) -> dict:
        if kind not in cls._cfg:
            path = jercjson if os.path.exists(jercjson) else jercjson + ".gz"
            cls._cfg[kind] = load_json_cached(path)
        return cls._cfg[kind]

    @classmethod
    def get_refs(cls, kind: str, year: str, isData: bool, era: Optional[str], jercjson: str):
        key = (kind, year, bool(isData), era or "")
        if key not in cls._refs:
            cfg = cls.get_cfg(kind, jercjson)
            if kind == "AK8" and year not in cfg:
                cfg = cls.get_cfg("AK4", jercjson)  # fallback
            tags = get_tag_names(cfg, year, isData, era)
            cls._refs[key] = CorrectionRefs(tags, jercjson)
        return cls._refs[key]



def phi_mpi_pi(dphi: float) -> float:
    """
    Wrap delta-phi into [-pi, pi].
    Works event-by-event in NanoAODTools.
    """
    return math.atan2(math.sin(dphi), math.cos(dphi))


def deltaR(eta1: float, phi1: float, eta2: float, phi2: float) -> float:
    """
    Compute DeltaR between two objects.
    """
    dEta = eta1 - eta2
    dPhi = phi_mpi_pi(phi1 - phi2)
    return math.sqrt(dEta * dEta + dPhi * dPhi)


def deltaR_obj(obj1, obj2) -> float:
    """
    Convenience wrapper: directly compute ΔR between NanoAODTools objects.
    Example: jet, fatjet, muon, etc.
    """
    return deltaR(obj1.eta, obj1.phi, obj2.eta, obj2.phi)

def eval_numeric(corr, *values) -> float:
    """
    Evaluate a correctionlib node with numeric arguments.
    Works on single jets/objects (event-by-event).
    """
    try:
        return float(corr.evaluate(*[float(v) for v in values]))
    except Exception as e:
        printDebug(f"[CorrectionLib] Numeric evaluation failed: {e}", indent=2)
        return 1.0


def eval_with_string_last(corr, *num_vals, cat: str = "nom") -> float:
    """
    Evaluate correctionlib node where categorical argument is LAST.
    Example: (eta, pt, 'Up'/'Down').
    """
    try:
        return float(corr.evaluate(*[float(v) for v in num_vals], cat))
    except Exception as e:
        printDebug(f"[CorrectionLib] Eval with string last failed: {e}", indent=2)
        return 1.0


def eval_with_string_first(corr, cat: str, *num_vals) -> float:
    """
    Evaluate correctionlib node where categorical argument is FIRST.
    Example: ('map_key', eta, phi).
    """
    try:
        return float(corr.evaluate(cat, *[float(v) for v in num_vals]))
    except Exception as e:
        printDebug(f"[CorrectionLib] Eval with string first failed: {e}", indent=2)
        return 1.0

def representative_run_number(year: str, run: int) -> float:
    """Representative run number for residual corrections."""
    if year == "2023Pre":
        return 367080.0
    if year == "2023Post":
        return 369803.0
    if year == "2024":
        return 379412.0
    return float(run)


def has_phi_dependent_L2(year: str) -> bool:
    """Some years require phi-dependent L2 corrections."""
    return year in ("2023Post", "2024")


def requires_run_based_residual(year: str) -> bool:
    """Years where residuals depend on run number."""
    return year in ("2023Pre", "2023Post", "2024")


def uses_puppi_met(year: str) -> bool:
    """Check if PuppiMET is the default MET type for this year."""
    return year in {"2022Pre","2022Post","2023Pre","2023Post","2024"}


def is_mc(isData: bool) -> bool:
    """Inverse of isData flag."""
    return not isData


def normal_from_seeds(seedA: int, seedB: int) -> float:
    """
    Deterministic per-jet Gaussian using Box–Muller transform.
    Based on hashed seeds from run/lumi/event/jet index.
    """
    UINT64_MAX = np.uint64(0xFFFFFFFFFFFFFFFF)

    a = (np.uint64(seedA) ^ np.uint64(0x9E3779B97F4A7C15)) * np.uint64(0xBF58476D1CE4E5B9) & UINT64_MAX
    b = (np.uint64(seedB) ^ np.uint64(0x94D049BB133111EB)) * np.uint64(0xD2B74407B1CE6E93) & UINT64_MAX

    # map to (0,1) avoiding 0
    u1 = (a.astype(np.float64) + 1.0) / (UINT64_MAX.astype(np.float64) + 2.0)
    u2 = (b.astype(np.float64) + 1.0) / (UINT64_MAX.astype(np.float64) + 2.0)

    # Box–Muller
    return math.sqrt(-2.0 * math.log(u1)) * math.cos(2.0 * math.pi * u2)


def build_jet_seeds(run: int, lumi: int, event: int, jet_index: int) -> Tuple[int, int]:
    """
    Build deterministic seeds for a jet in this event.
    Used for stochastic JER smearing.
    """
    seedA = (np.uint64(run) << np.uint64(1)) ^ np.uint64(lumi) ^ np.uint64(jet_index)
    seedB = (np.uint64(event) << np.uint64(1)) ^ np.uint64(jet_index * 1315423911)
    return int(seedA), int(seedB)


@dataclass
class Tags:
    """
    Container for JERC JSON tag names.
    These come from the year/era config JSON.
    """
    jercJsonPath: str
    tagNameL1FastJet: str
    tagNameL2Relative: str
    tagNameL3Absolute: str
    # tagNameL2Residual: str
    tagNameL2L3Residual: str
    tagNamePtResolution: str
    tagNameJerScaleFactor: str


def get_tag_names(base: dict, year: str, isData: bool, era: Optional[str]) -> Tags:
    """
    Build correct tag names from the Summer24Prompt24 convention.
    Uses flat keys like Summer24Prompt24_V1_MC_* or Summer24Prompt24_V1_DATA_*.
    """
    prefix = f"Summer{year}Prompt{year}_V1_{'DATA' if isData else 'MC'}"

    return Tags(
        jercJsonPath=None,  # already opened by CorrCache
        tagNameL1FastJet    = f"{prefix}_L1FastJet_AK4PFPuppi",
        tagNameL2Relative   = f"{prefix}_L2Relative_AK4PFPuppi",
        tagNameL3Absolute   = f"{prefix}_L3Absolute_AK4PFPuppi",
        # tagNameL2Residual   = f"{prefix}_L2Residual_AK4PFPuppi",
        tagNameL2L3Residual = f"{prefix}_L2L3Residual_AK4PFPuppi",
        # JER tags are still Summer23… until Summer24 JRV JSON is released
        tagNamePtResolution = "Summer23BPixPrompt23_RunD_JRV1_MC_PtResolution_AK4PFPuppi",
        tagNameJerScaleFactor = "Summer23BPixPrompt23_RunD_JRV1_MC_ScaleFactor_AK4PFPuppi",
    )



class CorrectionRefs:
    """
    Wrap correctionlib CorrectionSet and expose key nodes
    (L1, L2, Residuals, JER resolution and scale factors).
    """
    def __init__(self, tags: Tags, jercjson: str):
        # load from CorrCache’s parsed JSON instead of tags.jercJsonPath
        self.cs = correction.CorrectionSet.from_file(jercjson)
        self.cL1    = self.cs[tags.tagNameL1FastJet]
        self.cL2    = self.cs[tags.tagNameL2Relative]
        self.cL2L3Res = self.cs[tags.tagNameL2L3Residual]
        self.cReso  = self.cs[tags.tagNamePtResolution]
        self.cJerSF = self.cs[tags.tagNameJerScaleFactor]

@dataclass
class JerBin:
    """Region definition for JER smearing systematics."""
    label: str
    etaMin: float
    etaMax: float
    ptMin: float
    ptMax: float


def jes_nominal(jet, rho: float, year: str, isData: bool, run: int, refs: CorrectionRefs):
    """
    Compute the **nominal JES-corrected pT** for a given jet.

    Steps applied in order:
      1. Undo raw factor → recover uncorrected jet pT.
      2. Apply **L1FastJet correction**: area-based pileup subtraction.
      3. Apply **L2Relative correction**: relative η-dependent calibration.
         - For 2023Post/2024, includes phi-dependent corrections.
      4. Apply **L2Residual/L2L3Residual**:
         - For **Data only**: run-dependent residual JES correction.
         - For MC: skipped.

    Returns
    -------
    pt_corr : float
        Fully corrected jet pT after JES.
    pt_after_l1 : float
        Jet pT after only L1FastJet correction (before higher levels).
    """
    pt_raw = jet.pt * (1.0 - jet.rawFactor)

    # L1FastJet(area, eta, pt, rho)
    c1 = eval_numeric(refs.cL1, jet.area, jet.eta, pt_raw, rho)
    pt = pt_raw * c1
    pt_after_l1 = pt

    # L2Relative
    c2 = eval_numeric(refs.cL2, jet.eta, jet.phi, pt)
    pt *= c2

    # Residuals for Data
    if isData:
        if requires_run_based_residual(year):
            run_for_res = representative_run_number(year, run)
            cR = eval_numeric(refs.cL2L3Res, run_for_res, jet.eta, pt)
        else:
            cR = eval_numeric(refs.cL2L3Res, jet.eta, pt)
        pt *= cR

    return pt, pt_after_l1


def jes_syst_shift(pt_nom: float, eta: float, refs_cs: correction.CorrectionSet, syst_name: str, var: str) -> float:
    """
    Apply a **JES systematic variation** (Up/Down) to a nominal corrected pT.

    Parameters
    ----------
    pt_nom : float
        JES-corrected nominal jet pT.
    eta : float
        Jet pseudorapidity.
    refs_cs : CorrectionSet
        CorrectionLib set containing JES uncertainty sources.
    syst_name : str
        Name of JES source (e.g. 'jesTotal', 'jesAbsolute_2024').
    var : str
        Either 'Up' or 'Down'.

    Returns
    -------
    float
        Shifted pT value according to the JES uncertainty variation.
    """
    syst = refs_cs[syst_name]
    scale = eval_numeric(syst, eta, pt_nom)
    return pt_nom * (1.0 + scale) if var == "Up" else pt_nom * (1.0 - scale)

def apply_jes_variation(pt_nom: float, eta: float, refs_cs: correction.CorrectionSet,
                        sys: str, year_unc: str) -> float:
    """
    Apply JES systematic variation if it exists in the CorrectionSet.
    sys is like 'jesTotalUp', 'jesAbsoluteDown', 'jesHF_2024Up', etc.
    """
    if not sys or not refs_cs:
        return pt_nom

    # Extract base name and var
    var = "Up" if sys.endswith("Up") else "Down"
    syst_name = sys.replace("Up", "").replace("Down", "")

    # Handle year-specific JES (jesAbsolute_2024, jesHF_2024, etc.)
    if "_"+year_unc in syst_name:
        base = syst_name.split(f"_{year_unc}")[0]
        key = f"{base}_{year_unc}"
    else:
        key = syst_name

    if key not in refs_cs:
        printDebug(f"[apply_jes_variation] JES {key} not found in CorrectionSet")
        return pt_nom

    scale = eval_numeric(refs_cs[key], eta, pt_nom)
    return pt_nom * (1.0 + scale) if var == "Up" else pt_nom * (1.0 - scale)



def jer_smear(jet, rho: float, year: str, refs: CorrectionRefs,
              var: str, run: int, lumi: int, event: int,
              genjets, mindr: float, region: Optional[JerBin]):
    """
    Apply **Jet Energy Resolution (JER) smearing** to a jet.

    Logic:
      - First, obtain the nominal jet pT and JER resolution.
      - Compute JER scale factor (SF) for given η, pT, rho.
      - If jet matches a GenJet within ΔR < mindr and momentum balance criteria:
          → Apply "matched" correction: scale by (1 + (SF-1)*(pT-reco-pT-gen)/pT).
      - Otherwise:
          → Apply "stochastic smearing" using deterministic Gaussian seeds
            derived from run/lumi/event/jet index.

    Parameters
    ----------
    var : str
        'nom' for nominal, 'up' or 'down' for uncertainty variations.
    region : JerBin or None
        Defines η/pT region for region-dependent JER systematics.

    Returns
    -------
    float
        Smeared jet pT.
    """
    pt_in = jet.pt
    eta = jet.eta
    phi = jet.phi

    # Region gating
    if region:
        in_reg = (abs(eta) >= region.etaMin and abs(eta) < region.etaMax and
                  pt_in >= region.ptMin and pt_in < region.ptMax)
        use_var = var if in_reg else "nom"
    else:
        use_var = var

    # Resolution
    reso = eval_numeric(refs.cReso, eta, pt_in, rho)

    # Scale factor
    if uses_puppi_met(year):
        sf = eval_with_string_last(refs.cJerSF, eta, pt_in, cat=use_var)
    else:
        sf = eval_with_string_last(refs.cJerSF, eta, cat=use_var)

    # GenJet matching
    gpt, is_match = 0.0, False
    if hasattr(jet, "genJetIdx") and jet.genJetIdx >= 0 and jet.genJetIdx < len(genjets):
        gen = genjets[jet.genJetIdx]
        gpt = gen.pt
        dRmatch = deltaR(eta, phi, gen.eta, gen.phi)
        if dRmatch < mindr and abs(pt_in - gpt) < 3.0 * reso * pt_in:
            is_match = True

    # Matched correction
    if is_match:
        corr = max(0.0, 1.0 + (sf - 1.0) * (pt_in - gpt) / pt_in)
    else:
        # Stochastic correction (use seeds)
        seedA, seedB = build_jet_seeds(run, lumi, event, jet._index)  # _index = per-event jet index
        z = normal_from_seeds(seedA, seedB)
        sigma = math.sqrt(max(sf * sf - 1.0, 0.0)) * reso
        corr = max(0.0, 1.0 + z * sigma)

    return pt_in * corr

def apply_jer_variation(jet, rho: float, year: str, refs: CorrectionRefs,
                        sys: str, run: int, lumi: int, evt: int, genjets,
                        mindr: float, region: Optional[JerBin] = None) -> float:
    """
    Wrapper for applying JER systematics.

    Parameters
    ----------
    sys : str
        '', 'jerUp', or 'jerDown'.
    region : JerBin, optional
        Region bin for η/pT dependent systematics.

    Returns
    -------
    float
        Jet pT after JER nominal or systematic variation.
    """
    if not sys or sys == "":  # nominal
        return jet.pt

    jer_var = "up" if sys.endswith("Up") else "down"
    return jer_smear(jet, rho, year, refs, jer_var,
                     run, lumi, evt, genjets,
                     mindr=mindr, region=region)

def recompute_met(event_met_pt, event_met_phi, jets, jet_corr_pts):
    """
    Recompute MET_x and MET_y after applying jet corrections.
    Uses standard Type-1 propagation: subtract difference between
    corrected and raw jet pT from MET.
    """
    met_x = event_met_pt * math.cos(event_met_phi)
    met_y = event_met_pt * math.sin(event_met_phi)

    for j, pt_corr in zip(jets, jet_corr_pts):
        # ΔpT = (corrected - raw)
        delta_pt = pt_corr - j.pt
        met_x -= delta_pt * math.cos(j.phi)
        met_y -= delta_pt * math.sin(j.phi)

    met_pt_corr = math.sqrt(met_x**2 + met_y**2)
    met_phi_corr = math.atan2(met_y, met_x)
    return met_pt_corr, met_phi_corr



def get_syst_tag_names(year_json: dict, year: str) -> Dict[str, List[Tuple[str, str]]]:
    out = {}
    for key in year_json.keys():
        if key.startswith(f"Summer{year}Prompt{year}_V1_MC_"):
            base = key.replace(f"Summer{year}Prompt{year}_V1_MC_", "")
            out.setdefault("default", []).append((key, base))
    return out


def get_jer_sets(year_json: dict, year: str) -> Dict[str, List[JerBin]]:
    """
    Extract JER systematic bins from the year JSON.
    Returns dict: setName -> [JerBin(...), ...]
    """
    out = {}
    if year not in year_json:
        return out

    j = year_json[year].get("ForUncertaintyJER", {})
    for setName, obj in j.items():
        bins = []
        for label, arr in obj.items():
            if isinstance(arr, list) and len(arr) == 4:
                bins.append(JerBin(label, arr[0], arr[1], arr[2], arr[3]))
        out[setName] = bins
    return out


def correct_fatjet_softdropmass(fj, pt_corr_nom) -> float:
    """
    Recompute the corrected softdrop mass consistently with JES/JER corrections.

    Logic:
      - Undo rawFactor on pt and msoftdrop.
      - Apply the same scaling factor used for pt correction:
        scale = pt_corr_nom / pt_raw.
      - Ensures msoftdrop follows the same JES/JER corrections
        as the jet pt.

    Returns
    -------
    float
        Corrected softdrop mass for this fatjet.
    """
    msoftdrop_raw = fj.msoftdrop * (1.0 - fj.rawFactor)
    pt_raw = fj.pt * (1.0 - fj.rawFactor)

    # Use pt_corr_nom/pt_raw as base JES*JER*unc scaling
    scale = pt_corr_nom / pt_raw if pt_raw > 0 else 1.0
    return msoftdrop_raw * scale

def compute_globalParT3_mass(fatjets) -> list:
    """
    Compute globalParT3Xbb-corrected fatjet mass.

    Formula:
        gp3_mass = jet.globalParT3_massCorrX2p * (jet.mass * (1 - jet.rawFactor))

    - Uses nominal fatjet mass with rawFactor undone.
    - Appends only if the jet has attribute `globalParT3_massCorrX2p`.

    Parameters
    ----------
    fatjets : Collection
        NanoAODTools fatjet collection for the event.

    Returns
    -------
    list of float
        Corrected globalParT3 masses for each fatjet in the event.
    """
    gp3_masses = []
    for fj in fatjets:
        if hasattr(fj, "globalParT3_massCorrX2p"):
            raw_mass = fj.mass * (1.0 - fj.rawFactor)
            gp3_mass = fj.globalParT3_massCorrX2p * raw_mass
            gp3_masses.append(gp3_mass)
        else:
            gp3_masses.append(-999.0)  # fallback if missing
    return gp3_masses



def get_event_collections(event, isData: bool):
    """
    Return NanoAODTools Collections and scalars for an event.
    This replaces uproot-based array loading.
    """
    # Event-level info
    run   = event.run
    lumi  = event.luminosityBlock
    evt   = event.event
    rho   = event.Rho_fixedGridRhoFastjetAll

    # MET
    met   = event.MET_pt
    metphi= event.MET_phi
    rawMET     = getattr(event, "RawMET_pt", -1)
    rawMETphi  = getattr(event, "RawMET_phi", -99)
    rawPuppiMET= getattr(event, "RawPuppiMET_pt", -1)
    rawPuppiMETphi = getattr(event, "RawPuppiMET_phi", -99)

    # Jets
    jets    = Collection(event, "Jet")
    fatjets = Collection(event, "FatJet")

    # GenJets (only for MC)
    genjets   = Collection(event, "GenJet") if not isData else []
    genjetsAK8= Collection(event, "GenJetAK8") if not isData else []

    return {
        "run": run,
        "lumi": lumi,
        "evt": evt,
        "rho": rho,
        "met": met,
        "metphi": metphi,
        "rawMET": rawMET,
        "rawMETphi": rawMETphi,
        "rawPuppiMET": rawPuppiMET,
        "rawPuppiMETphi": rawPuppiMETphi,
        "jets": jets,
        "fatjets": fatjets,
        "genjets": genjets,
        "genjetsAK8": genjetsAK8,
    }

def branch_name_from_sys(sys: str, year_unc: str) -> str:
    """
    Map sys string into exact NanoAOD-style branch suffix.
    """
    if sys == "" or sys is None:
        return "nom"

    # TES passthrough (if added later)
    if sys in ("tesUp", "tesDown"):
        return sys

    # JER global
    if sys.startswith("jer_"):
        return "jerUp" if sys.endswith("up") else "jerDown"

    # JES variations already come as "jesSomethingUp/Down"
    return sys


class ApplyJerc(Module):
    def __init__(self, year, isData, era=None,
                 syst_kind="Nominal",
                 jes_base_tag=None,      # short, for branch naming
                 jes_full_tag_ak4=None,  # full, for correctionlib
                 jes_full_tag_ak8=None,
                 jes_var=None,
                 jer_bin=None,
                 jer_var=None,
                 year_unc="2024",
                 jercjson=None):
        self.year = year
        self.isData = isData
        self.era = era
        self.syst_kind = syst_kind
        self.jes_base_tag = jes_base_tag
        self.jes_full_tag_ak4 = jes_full_tag_ak4
        self.jes_full_tag_ak8 = jes_full_tag_ak8
        self.jes_var = jes_var
        self.jer_bin = jer_bin
        self.jer_var = jer_var
        self.year_unc = year_unc
        self.jercjson = jercjson


        self.refsAK4 = CorrCache.get_refs("AK4", year, isData, era, jercjson)
        self.refsAK8 = CorrCache.get_refs("AK8", year, isData, era, jercjson)

        # Cache correction sets too (avoid repeated .cs lookups)
        self.csAK4 = self.refsAK4.cs
        self.csAK8 = self.refsAK8.cs

        if syst_kind == "Nominal":
            self.jet_pt_systematics = [""]
        elif syst_kind == "JES" and jes_base_tag and jes_var:
            self.jet_pt_systematics = [f"{jes_base_tag}{jes_var}"]  # branch suffix
        elif syst_kind == "JER" and jer_bin and jer_var:
            self.jet_pt_systematics = [f"jer_{jer_bin.label}_{jer_var}"]
        else:
            self.jet_pt_systematics = [""]
            printDebug(f"[WARNING] Unknown syst_kind={syst_kind}, falling back to nominal")


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        # Define branches: one for each systematic
        for sys in self.jet_pt_systematics:
            suffix = branch_name_from_sys(sys, self.year_unc)
            self.out.branch(f"Jet_pt_{suffix}", "F", lenVar="nJet")
            self.out.branch(f"FatJet_pt_{suffix}", "F", lenVar="nFatJet")
            self.out.branch(f"FatJet_msoftdrop_{suffix}", "F", lenVar="nFatJet")


            if applyOnMET:
                self.out.branch(f"MET_pt_{suffix}", "F")
                self.out.branch(f"MET_phi_{suffix}", "F")

        self.out.branch("FatJet_globalParT3Xbb_mass", "F", lenVar="nFatJet")


    def analyze(self, event):
        ev = get_event_collections(event, self.isData)
        run, lumi, evt = ev["run"], ev["lumi"], ev["evt"]
        rho = ev["rho"]
        jets, fatjets = ev["jets"], ev["fatjets"]
        genjets, genjetsAK8 = ev["genjets"], ev["genjetsAK8"]

        gp3_masses = compute_globalParT3_mass(fatjets)
        self.out.fillBranch("FatJet_globalParT3Xbb_mass", gp3_masses)



        for sys in self.jet_pt_systematics:
            jet_corr_pts, fat_corr_pts = [], []

            # --- AK4 jets ---
            for j in jets:
                pt_corr, _ = jes_nominal(j, rho, self.year, self.isData, run, self.refsAK4)

                if self.syst_kind == "JES" and self.jes_full_tag_ak4:
                    pt_corr = jes_syst_shift(pt_corr, j.eta, self.csAK4,
                                            self.jes_full_tag_ak4, self.jes_var)

                if self.syst_kind == "JER" and is_mc(self.isData):
                    pt_corr = jer_smear(j, rho, self.year, self.refsAK4,
                                        self.jer_var, run, lumi, evt, genjets,
                                        mindr=0.2, region=self.jer_bin)

                jet_corr_pts.append(pt_corr)

            # --- AK8 jets ---
            fat_corr_msoft = []
            for fj in fatjets:
                pt_corr, _ = jes_nominal(fj, rho, self.year, self.isData, run, self.refsAK8)

                if self.syst_kind == "JES" and self.jes_full_tag_ak8:
                    pt_corr = jes_syst_shift(pt_corr, fj.eta,self.csAK8,
                                            self.jes_full_tag_ak8, self.jes_var)

                if self.syst_kind == "JER" and is_mc(self.isData):
                    pt_corr = jer_smear(fj, rho, self.year, self.refsAK8,
                                        self.jer_var, run, lumi, evt, genjetsAK8,
                                        mindr=0.4, region=self.jer_bin)

                msoft_corr = correct_fatjet_softdropmass(fj, pt_corr)

                fat_corr_msoft.append(msoft_corr)
                fat_corr_pts.append(pt_corr)



            suffix = branch_name_from_sys(sys, self.year_unc)
            self.out.fillBranch(f"Jet_pt_{suffix}", jet_corr_pts)
            self.out.fillBranch(f"FatJet_pt_{suffix}", fat_corr_pts)
            self.out.fillBranch(f"FatJet_msoftdrop_{suffix}", fat_corr_msoft)


            if applyOnMET:
                base_met = ev["rawPuppiMET"] if uses_puppi_met(self.year) else ev["rawMET"]
                base_met_phi = ev["rawPuppiMETphi"] if uses_puppi_met(self.year) else ev["rawMETphi"]
                met_pt_corr, met_phi_corr = recompute_met(base_met, base_met_phi, jets, jet_corr_pts)
                self.out.fillBranch(f"MET_pt_{suffix}", float(met_pt_corr))
                self.out.fillBranch(f"MET_phi_{suffix}", float(met_phi_corr))

        return True




def build_jerc_modules(year: str, isData: bool, era: Optional[str] = None, jercjson: str = None):
    """
    Build a list of ApplyJerc modules for Nominal, JES, JER systematics.
    Each module carries syst_kind + JES/JER tags explicitly.
    """
    modules = []

    # --- Nominal ---
    modules.append(
        ApplyJerc(year, isData, syst_kind="Nominal", year_unc=year, jes_full_tag_ak4=None,
                  jes_full_tag_ak8=None, jercjson=jercjson)
    )

    # --- JES systematics (only for MC) ---
    if is_mc(isData):
        cfgAK4 = CorrCache.get_cfg("AK4", jercjson)
        cfgAK8 = CorrCache.get_cfg("AK8", jercjson)

        sAK4 = get_syst_tag_names(cfgAK4, year)
        sAK8 = get_syst_tag_names(cfgAK8, year) if year in cfgAK8 else sAK4

        base_to_ak4 = {base: full for setName, pairs in sAK4.items() for (full, base) in pairs}
        base_to_ak8 = {base: full for setName, pairs in sAK8.items() for (full, base) in pairs}
        common_bases = sorted(set(base_to_ak4) & set(base_to_ak8))

        for base in common_bases:  # base = short name, full = correctionlib key
            for var in ("Up", "Down"):
                modules.append(
                    ApplyJerc(year, isData,
                                    syst_kind="JES",
                                    jes_base_tag=base,
                                    jes_full_tag_ak4=base_to_ak4[base],
                                    jes_full_tag_ak8=base_to_ak8[base],
                                    jes_var=var,
                                    jercjson=jercjson)
                )
    # --- JER systematics (only for MC) ---
    if is_mc(isData):
        for var in ("up", "down"):
            modules.append(
                ApplyJerc(year, isData, 
                                syst_kind="JER",
                                jer_var=var,
                                jercjson=jercjson)
            )
    return modules
