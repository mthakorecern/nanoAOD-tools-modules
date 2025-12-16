#!/usr/bin/env python3
import os
import json
import re
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module


class XSWeightOnly(Module):
    def __init__(self, filename, year="2024", isData=False):
        self.year = str(year)
        self.isData = isData
        self.isMC = not isData
        self.filename = filename
        print(f"[INIT] XSWeightOnly for file: {self.filename} | Year: {self.year}")

        # Integrated luminosity (pb^-1)
        self.LHCLumi = 109.08e3 if self.year == "2024" else 1.0
        print(f"[INFO] Set Lumi = {self.LHCLumi} pb^-1")

        # List of JSON files
        self.xs_json_list = [
            f"{self.year}_DIB_samples.json",
            f"{self.year}_DY_PT_Binned.json",
            # f"{self.year}_DY_samples.json",
            f"{self.year}_QCD_samples.json",
            f"{self.year}_Radion_samples.json",
            f"{self.year}_STop_samples.json",
            f"{self.year}_TTbar_samples.json",
            f"{self.year}_WJets_PT_Binned.json",
            # f"{self.year}_WJets_samples.json",
            


        ]
        self.counts_json_list = [
            # f"{self.year}_new_processed_samples.json",
            f"{self.year}_CRAB_processed.json",
  

        ]

        # Containers for data
        self.xsInfo = {}
        self.countsInfo = {}

    def beginJob(self):
        if self.isMC:
            # Load all XS JSONs
            for js in self.xs_json_list:
                if os.path.exists(js):
                    with open(js, "r") as f:
                        data = json.load(f)
                        self.xsInfo.update(data)
                        print(f"[LOAD] XS JSON loaded: {js} ({len(data)} entries)")
                else:
                    print(f"[WARN] Missing XS JSON: {js}")

            # Load all genweight JSONs
            for js in self.counts_json_list:
                if os.path.exists(js):
                    with open(js, "r") as f:
                        data = json.load(f)
                        self.countsInfo.update(data)
                        print(f"[LOAD] GenWeight JSON loaded: {js} ({len(data)} entries)")
                else:
                    print(f"[WARN] Missing GenWeight JSON: {js}")

    def endJob(self):
        print(f"[DONE] Finished job for year {self.year}")

    def clean_basename(self, filename):
        """
        Get dataset base name, both with and without index.
        Example:
            QCD-HT1000to1200_12.root -> 
                base_no_index = QCD-HT1000to1200
                base_with_index = QCD-HT1000to1200_12
        """
        base = os.path.basename(filename).replace(".root", "")
        match = re.match(r"(.+)_([0-9]+)$", base)
        if match:
            base_no_index = match.group(1)
            base_with_index = base
        else:
            base_no_index = base
            base_with_index = base
        return base_no_index, base_with_index

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if not self.isMC:
            return

        full_input_name = inputFile.GetName()
        print(f"[FILE] Processing input file: {full_input_name}")

        base_no_idx, base_with_idx = self.clean_basename(full_input_name)
        print(f"[INFO] Cleaned names -> base_no_idx: {base_no_idx}, base_with_idx: {base_with_idx}")

        # --- Cross-section lookup ---
        if base_no_idx in self.xsInfo:
            self.XS = self.xsInfo[base_no_idx]["XS"]
        else:
            raise KeyError(f"[ERROR] Cross-section not found for dataset: {base_no_idx}")

        # --- GenWeight lookup ---
        if base_with_idx in self.countsInfo:
            self.totalNumberOfEvents = self.countsInfo[base_with_idx]
            found_key = base_with_idx
        elif base_no_idx in self.countsInfo:
            self.totalNumberOfEvents = self.countsInfo[base_no_idx]
            found_key = base_no_idx
        else:
            raise KeyError(f"[ERROR] GenWeight sum not found for dataset: {base_with_idx} or {base_no_idx}")

        print(f"[INFO] Dataset: {base_no_idx}")
        print(f"       XS = {self.XS} pb")
        print(f"       Sum of gen weights = {self.totalNumberOfEvents} (from key: {found_key})")

        # Output branches
        self.out.branch("xsWeight", "F")
        self.out.branch("xs", "F")
        self.out.branch("lhclumi", "F")
        self.out.branch("lhclumi_pb", "F")
        self.out.branch("sumofgenwts", "F")

    def analyze(self, event):
        if not self.isMC:
            return True

        xsWeight = (self.XS * self.LHCLumi * event.genWeight) / self.totalNumberOfEvents

        self.out.fillBranch("xsWeight", xsWeight)
        self.out.fillBranch("xs", self.XS)
        self.out.fillBranch("lhclumi", self.LHCLumi * 1e-3)  # fb^-1
        self.out.fillBranch("lhclumi_pb", self.LHCLumi)
        self.out.fillBranch("sumofgenwts", self.totalNumberOfEvents)
        return True
