import os
import ROOT
from ROOT import TH1, TFile

ROOT.PyConfig.IgnoreCommandLineOptions = True
import argparse
import glob
import multiprocessing as mp
from re import search, IGNORECASE

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

def get2016EWKW(pt):
    weight = 1
    if pt < 170:
        weight = 0.94269
    elif 170 <= pt < 200:
        weight = 0.902615
    elif 200 <= pt < 230:
        weight = 0.898827
    elif 230 <= pt < 260:
        weight = 0.959081
    elif 260 <= pt < 290:
        weight = 0.891248
    elif 290 <= pt < 320:
        weight = 0.860188
    elif 320 <= pt < 350:
        weight = 0.884811
    elif 350 <= pt < 390:
        weight = 0.868131
    elif 390 <= pt < 430:
        weight = 0.848655
    elif 430 <= pt < 470:
        weight = 0.806186
    elif 470 <= pt < 510:
        weight = 0.848507
    elif 510 <= pt < 550:
        weight = 0.83763
    elif 550 <= pt < 590:
        weight = 0.792152
    elif 590 <= pt < 640:
        weight = 0.730731
    elif 640 <= pt < 690:
        weight = 0.778061
    elif 690 <= pt < 740:
        weight = 0.771811
    elif 740 <= pt < 790:
        weight = 0.795004
    elif 790 <= pt < 840:
        weight = 0.757859
    elif 840 <= pt < 900:
        weight = 0.709571
    elif 900 <= pt < 960:
        weight = 0.702751
    elif 960 <= pt < 1020:
        weight = 0.657821
    elif 1020 <= pt < 1090:
        weight = 0.762559
    elif 1090 <= pt < 1160:
        weight = 0.845925
    elif pt >= 1160:
        weight = 0.674034

    return weight

def get2016QCDW(pt):
    weight = 1
    if pt < 170:
        weight = 1.43896
    elif 170 <= pt < 200:
        weight = 1.45307
    elif 200 <= pt < 230:
        weight = 1.41551
    elif 230 <= pt < 260:
        weight = 1.42199
    elif 260 <= pt < 290:
        weight = 1.3477
    elif 290 <= pt < 320:
        weight = 1.35302
    elif 320 <= pt < 350:
        weight = 1.34289
    elif 350 <= pt < 390:
        weight = 1.32474
    elif 390 <= pt < 430:
        weight = 1.23267
    elif 430 <= pt < 470:
        weight = 1.22641
    elif 470 <= pt < 510:
        weight = 1.23149
    elif 510 <= pt < 550:
        weight = 1.21593
    elif 550 <= pt < 590:
        weight = 1.16506
    elif 590 <= pt < 640:
        weight = 1.01718
    elif 640 <= pt < 690:
        weight = 1.01575
    elif 690 <= pt < 740:
        weight = 1.05425
    elif 740 <= pt < 790:
        weight = 1.05992
    elif 790 <= pt < 840:
        weight = 1.01503
    elif 840 <= pt < 900:
        weight = 1.01761
    elif 900 <= pt < 960:
        weight = 0.947194
    elif 960 <= pt < 1020:
        weight = 0.932754
    elif 1020 <= pt < 1090:
        weight = 1.00849
    elif 1090 <= pt < 1160:
        weight = 0.94805
    elif pt >= 1160:
        weight = 0.86956

    return weight

def get2016RenUpW(pt):
    weight = 1
    if pt < 170:
        weight = 1.33939
    elif 170 <= pt < 200:
        weight = 1.34649
    elif 200 <= pt < 230:
        weight = 1.3115
    elif 230 <= pt < 260:
        weight = 1.30782
    elif 260 <= pt < 290:
        weight = 1.24619
    elif 290 <= pt < 320:
        weight = 1.24531
    elif 320 <= pt < 350:
        weight = 1.23472
    elif 350 <= pt < 390:
        weight = 1.20709
    elif 390 <= pt < 430:
        weight = 1.13368
    elif 430 <= pt < 470:
        weight = 1.12227
    elif 470 <= pt < 510:
        weight = 1.12612
    elif 510 <= pt < 550:
        weight = 1.1118
    elif 550 <= pt < 590:
        weight = 1.06549
    elif 590 <= pt < 640:
        weight = 0.931838
    elif 640 <= pt < 690:
        weight = 0.929282
    elif 690 <= pt < 740:
        weight = 0.959553
    elif 740 <= pt < 790:
        weight = 0.955823
    elif 790 <= pt < 840:
        weight = 0.920614
    elif 840 <= pt < 900:
        weight = 0.917243
    elif 900 <= pt < 960:
        weight = 0.855649
    elif 960 <= pt < 1020:
        weight = 0.84587
    elif 1020 <= pt < 1090:
        weight = 0.906862
    elif 1090 <= pt < 1160:
        weight = 0.858763
    elif pt >= 1160:
        weight = 0.794909

    return weight

def get2016RenDownW(pt):
    weight = 1
    if pt < 170:
        weight = 1.52924
    elif 170 <= pt < 200:
        weight = 1.55189
    elif 200 <= pt < 230:
        weight = 1.5104
    elif 230 <= pt < 260:
        weight = 1.53402
    elif 260 <= pt < 290:
        weight = 1.44328
    elif 290 <= pt < 320:
        weight = 1.45859
    elif 320 <= pt < 350:
        weight = 1.4486
    elif 350 <= pt < 390:
        weight = 1.44818
    elif 390 <= pt < 430:
        weight = 1.33058
    elif 430 <= pt < 470:
        weight = 1.33239
    elif 470 <= pt < 510:
        weight = 1.33976
    elif 510 <= pt < 550:
        weight = 1.32336
    elif 550 <= pt < 590:
        weight = 1.26854
    elif 590 <= pt < 640:
        weight = 1.10489
    elif 640 <= pt < 690:
        weight = 1.1055
    elif 690 <= pt < 740:
        weight = 1.15559
    elif 740 <= pt < 790:
        weight = 1.17682
    elif 790 <= pt < 840:
        weight = 1.118
    elif 840 <= pt < 900:
        weight = 1.13097
    elif 900 <= pt < 960:
        weight = 1.04988
    elif 960 <= pt < 1020:
        weight = 1.02796
    elif 1020 <= pt < 1090:
        weight = 1.12438
    elif 1090 <= pt < 1160:
        weight = 1.04704
    elif pt >= 1160:
        weight = 0.94791

    return weight

def get2016FacUpW(pt):
    weight = 1
    if pt < 170:
        weight = 1.44488
    elif 170 <= pt < 200:
        weight = 1.45847
    elif 200 <= pt < 230:
        weight = 1.41004
    elif 230 <= pt < 260:
        weight = 1.39864
    elif 260 <= pt < 290:
        weight = 1.3391
    elif 290 <= pt < 320:
        weight = 1.33663
    elif 320 <= pt < 350:
        weight = 1.32073
    elif 350 <= pt < 390:
        weight = 1.3076
    elif 390 <= pt < 430:
        weight = 1.20904
    elif 430 <= pt < 470:
        weight = 1.20066
    elif 470 <= pt < 510:
        weight = 1.20462
    elif 510 <= pt < 550:
        weight = 1.18286
    elif 550 <= pt < 590:
        weight = 1.12586
    elif 590 <= pt < 640:
        weight = 0.990615
    elif 640 <= pt < 690:
        weight = 0.984473
    elif 690 <= pt < 740:
        weight = 1.0171
    elif 740 <= pt < 790:
        weight = 1.01706
    elif 790 <= pt < 840:
        weight = 0.986107
    elif 840 <= pt < 900:
        weight = 0.972452
    elif 900 <= pt < 960:
        weight = 0.910183
    elif 960 <= pt < 1020:
        weight = 0.885284
    elif 1020 <= pt < 1090:
        weight = 0.950662
    elif 1090 <= pt < 1160:
        weight = 0.89605
    elif pt >= 1160:
        weight = 0.823212

    return weight

def get2016FacDownW(pt):
    weight = 1
    if pt < 170:
        weight = 1.43694
    elif 170 <= pt < 200:
        weight = 1.45181
    elif 200 <= pt < 230:
        weight = 1.42649
    elif 230 <= pt < 260:
        weight = 1.4513
    elif 260 <= pt < 290:
        weight = 1.3624
    elif 290 <= pt < 320:
        weight = 1.37472
    elif 320 <= pt < 350:
        weight = 1.3718
    elif 350 <= pt < 390:
        weight = 1.34637
    elif 390 <= pt < 430:
        weight = 1.26276
    elif 430 <= pt < 470:
        weight = 1.25762
    elif 470 <= pt < 510:
        weight = 1.26389
    elif 510 <= pt < 550:
        weight = 1.25504
    elif 550 <= pt < 590:
        weight = 1.2116
    elif 590 <= pt < 640:
        weight = 1.04846
    elif 640 <= pt < 690:
        weight = 1.052
    elif 690 <= pt < 740:
        weight = 1.09788
    elif 740 <= pt < 790:
        weight = 1.10967
    elif 790 <= pt < 840:
        weight = 1.046
    elif 840 <= pt < 900:
        weight = 1.06961
    elif 900 <= pt < 960:
        weight = 0.989001
    elif 960 <= pt < 1020:
        weight = 0.987722
    elif 1020 <= pt < 1090:
        weight = 1.07552
    elif 1090 <= pt < 1160:
        weight = 1.00802
    elif pt >= 1160:
        weight = 0.922172

    return weight

###########################################################################################

def get2016EWKZ(pt):
    weight = 1
    if pt < 170:
        weight = 0.970592
    elif 170 <= pt < 200:
        weight = 0.964424
    elif 200 <= pt < 230:
        weight = 0.956695
    elif 230 <= pt < 260:
        weight = 0.948747
    elif 260 <= pt < 290:
        weight = 0.941761
    elif 290 <= pt < 320:
        weight = 0.934246
    elif 320 <= pt < 350:
        weight = 0.927089
    elif 350 <= pt < 390:
        weight = 0.919181
    elif 390 <= pt < 430:
        weight = 0.909926
    elif 430 <= pt < 470:
        weight = 0.900911
    elif 470 <= pt < 510:
        weight = 0.892561
    elif 510 <= pt < 550:
        weight = 0.884353
    elif 550 <= pt < 590:
        weight = 0.8761
    elif 590 <= pt < 640:
        weight = 0.867687
    elif 640 <= pt < 690:
        weight = 0.858047
    elif 690 <= pt < 740:
        weight = 0.849014
    elif 740 <= pt < 790:
        weight = 0.840317
    elif 790 <= pt < 840:
        weight = 0.832017
    elif 840 <= pt < 900:
        weight = 0.823545
    elif 900 <= pt < 960:
        weight = 0.814596
    elif 960 <= pt < 1020:
        weight = 0.806229
    elif 1020 <= pt < 1090:
        weight = 0.798038
    elif 1090 <= pt < 1160:
        weight = 0.789694
    elif pt >= 1160:
        weight = 0.781163

    return weight


def get2016QCDZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.47528
    elif 170 <= pt < 200:
        weight = 1.5428
    elif 200 <= pt < 230:
        weight = 1.49376
    elif 230 <= pt < 260:
        weight = 1.39119
    elif 260 <= pt < 290:
        weight = 1.40538
    elif 290 <= pt < 320:
        weight = 1.44661
    elif 320 <= pt < 350:
        weight = 1.38176
    elif 350 <= pt < 390:
        weight = 1.37381
    elif 390 <= pt < 430:
        weight = 1.29145
    elif 430 <= pt < 470:
        weight = 1.33452
    elif 470 <= pt < 510:
        weight = 1.25765
    elif 510 <= pt < 550:
        weight = 1.24265
    elif 550 <= pt < 590:
        weight = 1.24331
    elif 590 <= pt < 640:
        weight = 1.16187
    elif 640 <= pt < 690:
        weight = 1.07349
    elif 690 <= pt < 740:
        weight = 1.10748
    elif 740 <= pt < 790:
        weight = 1.06617
    elif 790 <= pt < 840:
        weight = 1.05616
    elif 840 <= pt < 900:
        weight = 1.1149
    elif 900 <= pt < 960:
        weight = 1.03164
    elif 960 <= pt < 1020:
        weight = 1.06872
    elif 1020 <= pt < 1090:
        weight = 0.981645
    elif 1090 <= pt < 1160:
        weight = 0.81729
    elif pt >= 1160:
        weight = 0.924246

    return weight


def get2016RenUpZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.38864
    elif 170 <= pt < 200:
        weight = 1.43444
    elif 200 <= pt < 230:
        weight = 1.39055
    elif 230 <= pt < 260:
        weight = 1.2934
    elif 260 <= pt < 290:
        weight = 1.30618
    elif 290 <= pt < 320:
        weight = 1.33453
    elif 320 <= pt < 350:
        weight = 1.27239
    elif 350 <= pt < 390:
        weight = 1.25981
    elif 390 <= pt < 430:
        weight = 1.18738
    elif 430 <= pt < 470:
        weight = 1.21686
    elif 470 <= pt < 510:
        weight = 1.14351
    elif 510 <= pt < 550:
        weight = 1.13305
    elif 550 <= pt < 590:
        weight = 1.12346
    elif 590 <= pt < 640:
        weight = 1.05153
    elif 640 <= pt < 690:
        weight = 0.98056
    elif 690 <= pt < 740:
        weight = 1.00651
    elif 740 <= pt < 790:
        weight = 0.96642
    elif 790 <= pt < 840:
        weight = 0.956161
    elif 840 <= pt < 900:
        weight = 0.998079
    elif 900 <= pt < 960:
        weight = 0.927421
    elif 960 <= pt < 1020:
        weight = 0.954094
    elif 1020 <= pt < 1090:
        weight = 0.883954
    elif 1090 <= pt < 1160:
        weight = 0.729348
    elif pt >= 1160:
        weight = 0.833531

    return weight


def get2016RenDownZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.55029
    elif 170 <= pt < 200:
        weight = 1.64816
    elif 200 <= pt < 230:
        weight = 1.59095
    elif 230 <= pt < 260:
        weight = 1.48394
    elif 260 <= pt < 290:
        weight = 1.49876
    elif 290 <= pt < 320:
        weight = 1.55923
    elif 320 <= pt < 350:
        weight = 1.49283
    elif 350 <= pt < 390:
        weight = 1.49279
    elif 390 <= pt < 430:
        weight = 1.39829
    elif 430 <= pt < 470:
        weight = 1.46136
    elif 470 <= pt < 510:
        weight = 1.38252
    elif 510 <= pt < 550:
        weight = 1.36081
    elif 550 <= pt < 590:
        weight = 1.37848
    elif 590 <= pt < 640:
        weight = 1.28575
    elif 640 <= pt < 690:
        weight = 1.17335
    elif 690 <= pt < 740:
        weight = 1.21841
    elif 740 <= pt < 790:
        weight = 1.17707
    elif 790 <= pt < 840:
        weight = 1.16846
    elif 840 <= pt < 900:
        weight = 1.25225
    elif 900 <= pt < 960:
        weight = 1.15246
    elif 960 <= pt < 1020:
        weight = 1.20461
    elif 1020 <= pt < 1090:
        weight = 1.09389
    elif 1090 <= pt < 1160:
        weight = 0.921666
    elif pt >= 1160:
        weight = 1.02897

    return weight


def get2016FacUpZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.48202
    elif 170 <= pt < 200:
        weight = 1.54551
    elif 200 <= pt < 230:
        weight = 1.49889
    elif 230 <= pt < 260:
        weight = 1.37781
    elif 260 <= pt < 290:
        weight = 1.39565
    elif 290 <= pt < 320:
        weight = 1.42918
    elif 320 <= pt < 350:
        weight = 1.35817
    elif 350 <= pt < 390:
        weight = 1.34966
    elif 390 <= pt < 430:
        weight = 1.26377
    elif 430 <= pt < 470:
        weight = 1.30245
    elif 470 <= pt < 510:
        weight = 1.22426
    elif 510 <= pt < 550:
        weight = 1.20589
    elif 550 <= pt < 590:
        weight = 1.20726
    elif 590 <= pt < 640:
        weight = 1.12468
    elif 640 <= pt < 690:
        weight = 1.04086
    elif 690 <= pt < 740:
        weight = 1.06898
    elif 740 <= pt < 790:
        weight = 1.02577
    elif 790 <= pt < 840:
        weight = 1.01591
    elif 840 <= pt < 900:
        weight = 1.06486
    elif 900 <= pt < 960:
        weight = 0.986133
    elif 960 <= pt < 1020:
        weight = 1.01224
    elif 1020 <= pt < 1090:
        weight = 0.937427
    elif 1090 <= pt < 1160:
        weight = 0.77774
    elif pt >= 1160:
        weight = 0.87435

    return weight


def get2016FacDownZ(pt):
    weight = 1
    if pt < 170:
        weight = 1.47257
    elif 170 <= pt < 200:
        weight = 1.54471
    elif 200 <= pt < 230:
        weight = 1.49357
    elif 230 <= pt < 260:
        weight = 1.41154
    elif 260 <= pt < 290:
        weight = 1.42113
    elif 290 <= pt < 320:
        weight = 1.46933
    elif 320 <= pt < 350:
        weight = 1.41145
    elif 350 <= pt < 390:
        weight = 1.40379
    elif 390 <= pt < 430:
        weight = 1.32597
    elif 430 <= pt < 470:
        weight = 1.37278
    elif 470 <= pt < 510:
        weight = 1.2973
    elif 510 <= pt < 550:
        weight = 1.28653
    elif 550 <= pt < 590:
        weight = 1.28624
    elif 590 <= pt < 640:
        weight = 1.20542
    elif 640 <= pt < 690:
        weight = 1.11142
    elif 690 <= pt < 740:
        weight = 1.1518
    elif 740 <= pt < 790:
        weight = 1.11242
    elif 790 <= pt < 840:
        weight = 1.1022
    elif 840 <= pt < 900:
        weight = 1.17247
    elif 900 <= pt < 960:
        weight = 1.08313
    elif 960 <= pt < 1020:
        weight = 1.13457
    elif 1020 <= pt < 1090:
        weight = 1.03169
    elif 1090 <= pt < 1160:
        weight = 0.862166
    elif pt >= 1160:
        weight = 0.981254

    return weight

def ensureTFile(filename, option="READ"):
    """Open TFile, checking if the file in the given path exists."""
    if not os.path.isfile(filename):
        print(
            (
                '>>> ERROR! ScaleFactorTool.ensureTFile: File in path "%s" does not exist!'
                % (filename)
            )
        )
        exit(1)
    file = TFile(filename, option)
    if not file or file.IsZombie():
        print(
            (
                '>>> ERROR! ScaleFactorTool.ensureTFile Could not open file by name "%s"'
                % (filename)
            )
        )
        exit(1)
    return file


def ensureFile(*paths, **kwargs):
    """Ensure file exists."""
    filepath = os.path.join(*paths)
    stop = kwargs.get("stop", True)
    if "*" in filepath or "?" in filepath:
        exists = len(glob.glob(filepath)) > 0
    else:
        exists = os.path.isfile(filepath)
    if not exists and stop:
        raise OSError('File "%s" does not exist' % (filepath))
    return filepath


def extractTH1(file, histname):
    """Get histogram from a file, and do SetDirectory(0)."""
    if not file or file.IsZombie():
        print(">>> ERROR! ScaleFactorTool.extractTH1 Could not open file!")
        exit(1)
    hist = file.Get(histname)
    if not hist:
        print(
            (
                '>>> ERROR! ScaleFactorTool.extractTH1: Did not find histogtam "%s" in file %s!'
                % (histname, file.GetName())
            )
        )
        exit(1)
    if isinstance(hist, TH1):
        hist.SetDirectory(0)
    return hist


def warning(string, **kwargs):
    """Print warning with color."""
    pre = kwargs.get("pre", "") + "\033[1m\033[93mWarning!\033[0m \033[93m"
    title = kwargs.get("title", "")
    if title:
        pre = "%s%s: " % (pre, title)
    string = "%s%s\033[0m" % (pre, string)
    print((string.replace("\n", "\n" + " " * (len(pre) - 18))))


class KFactorTool:
    def __init__(self, year="2024"):
        self.year = year
        self.f_nlo_qcd = ensureTFile("ZorW_NLO_corrections.root")

        self.nlo_ewk_WJets = self.f_nlo_qcd.Get("EWK_W_nominal")
        self.nlo_qcd_WJets = self.f_nlo_qcd.Get("QCD_W")
        self.nlo_ewk_ZJets = self.f_nlo_qcd.Get("EWK_Z_nominal")
        self.nlo_qcd_ZJets = self.f_nlo_qcd.Get("QCD_Z")

    def getEWKW(self, pt):
        ##return get2016EWKW(pt)
        # if self.year == 2016:
        #     return get2016EWKW(pt)

        bin = self.nlo_ewk_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_ewk_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_ewk_WJets.GetBinContent(bin)
        return sf

    def getDeborahEWKW(self, pt):
        return get2016EWKW(pt)

    def getEWKZ(self, pt):
        bin = self.nlo_ewk_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_ewk_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_ewk_ZJets.GetBinContent(bin)
        return sf

    def getDeborahEWKZ(self, pt):
        return get2016EWKZ(pt)

    def getQCDW(self, pt):
        bin = self.nlo_qcd_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_WJets.GetBinContent(bin)
        return sf

    def getQCDZTo2Nu(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        return sf

    def getQCDZTo2L(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        return sf

    def getRenUpW(self, pt):
        bin = self.nlo_qcd_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_WJets.GetBinContent(bin)
        ratio = get2016RenUpW(pt) / get2016QCDW(pt)
        return sf * ratio

    def getRenDownW(self, pt):
        bin = self.nlo_qcd_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_WJets.GetBinContent(bin)
        ratio = get2016RenDownW(pt) / get2016QCDW(pt)
        return sf * ratio

    def getFacUpW(self, pt):
        bin = self.nlo_qcd_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_WJets.GetBinContent(bin)
        ratio = get2016FacUpW(pt) / get2016QCDW(pt)
        return sf * ratio

    def getFacDownW(self, pt):
        bin = self.nlo_qcd_WJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_WJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_WJets.GetBinContent(bin)
        ratio = get2016FacDownW(pt) / get2016QCDW(pt)
        return sf * ratio

    def getRenUpZTo2Nu(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016RenUpZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getRenDownZTo2Nu(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016RenDownZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getFacUpZTo2Nu(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016FacUpZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getFacDownZTo2Nu(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016FacDownZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getRenUpZTo2L(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016RenUpZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getRenDownZTo2L(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016RenDownZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getFacUpZTo2L(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016FacUpZ(pt) / get2016QCDZ(pt)
        return sf * ratio

    def getFacDownZTo2L(self, pt):
        bin = self.nlo_qcd_ZJets.GetXaxis().FindBin(pt)
        if bin == 0:
            bin = 1
        elif bin > self.nlo_qcd_ZJets.GetXaxis().GetNbins():
            bin -= 1
        sf = self.nlo_qcd_ZJets.GetBinContent(bin)
        ratio = get2016FacDownZ(pt) / get2016QCDZ(pt)
        return sf * ratio


class wandzgenptweight(Module):
    def __init__(self, filename, year="2024", doSysVar=False):
        print("ACTIAVTE: wandzgenptweight")
        self.filename = filename
        self.year = year
        self.kFactorTool = KFactorTool(year=self.year)
        self.doSysVar = doSysVar


    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        # Define output branches to check if preselection is working
        self.eventCount = 0
        self.out = wrappedOutputTree
        self.out.branch("ewkWWeight", "F")
        self.out.branch("ewkZWeight", "F")
        self.out.branch("ewkWDeborahsWeight", "F")
        self.out.branch("ewkZDeborahsWeight", "F")
        self.out.branch("qcdWWeight", "F")
        self.out.branch("qcdZTo2LWeight", "F")
        self.out.branch("GenV_pt", "F")
        self.out.branch("WnnloWeight", "F")
        self.out.branch("ZnnloWeight", "F")
        self.out.branch("combinedWZgenPtWeight", "F")
        self.out.branch("combinedWZgenPtDeborahWeight", "F")
        self.out.branch("combinedWZgenPtDylanWeight", "F")

        if self.doSysVar:
            self.out.branch("qcdWWeightRenUp", "F")
            self.out.branch("qcdWWeightRenDown", "F")
            self.out.branch("qcdWWeightFacUp", "F")
            self.out.branch("qcdWWeightFacDown", "F")
            self.out.branch("qcdZTo2LWeightRenUp", "F")
            self.out.branch("qcdZTo2LWeightRenDown", "F")
            self.out.branch("qcdZTo2LWeightFacUp", "F")
            self.out.branch("qcdZTo2LWeightFacDown", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        ewkWWeight = ewkZWeight = qcdWWeight = qcdZTo2LWeight = ewkWDeborahsWeight = ewkZDeborahsWeight = combinedWZgenPtDylanWeight = WnnloWeight = ZnnloWeight = combinedWZgenPtWeight = combinedWZgenPtDeborahWeight = 1
        qcdWWeightRenUp = qcdWWeightRenDown = qcdWWeightFacUp = qcdWWeightFacDown = 1
        qcdZTo2LWeightRenUp = qcdZTo2LWeightRenDown = qcdZTo2LWeightFacUp = (qcdZTo2LWeightFacDown) = 1
        GenV_pt = -9
        GenV = []
        sample = os.path.basename(self.filename)
        if search(r"^DYto2L", sample, IGNORECASE):
            genParticles = Collection(event, "GenPart")
            GenV = list(
                [
                    gen
                    for gen in genParticles
                    if (gen.pdgId == 23) and (gen.status == 22)
                ]
            )
            if len(GenV) > 0:
                GenV_pt = GenV[0].pt
                ewkZWeight *= self.kFactorTool.getEWKZ(GenV_pt)
                ewkZDeborahsWeight *= self.kFactorTool.getDeborahEWKZ(GenV_pt)
                qcdZTo2LWeight *= self.kFactorTool.getQCDZTo2L(GenV_pt)
                qcdZTo2LWeightRenUp *= self.kFactorTool.getRenUpZTo2L(GenV_pt)
                qcdZTo2LWeightRenDown *= self.kFactorTool.getRenDownZTo2L(GenV_pt)
                qcdZTo2LWeightFacUp *= self.kFactorTool.getFacUpZTo2L(GenV_pt)
                qcdZTo2LWeightFacDown *= self.kFactorTool.getFacDownZTo2L(GenV_pt)
                ZnnloWeight = 0.934

                # if self.year == "2016APV":
                #     combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight
                #     combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight
                # elif self.year == "2016":
                #     combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight
                #     combinedWZgenPtDeborahWeight = qcdZTo2LWeight * ewkZDeborahsWeight
                # elif self.year == "2017":
                #     combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight * 0.934
                #     combinedWZgenPtDeborahWeight = (qcdZTo2LWeight * ewkZDeborahsWeight * 0.934)
                if self.year == "2024":
                    combinedWZgenPtWeight = qcdZTo2LWeight * ewkZWeight * 0.934
                    combinedWZgenPtDeborahWeight = (qcdZTo2LWeight * ewkZDeborahsWeight * 0.934)

            combinedWZgenPtDylanWeight = 1.23
        # combinedWZgenPtWeight = qcdZTo2LWeight*ewkZWeight*0.934
        # combinedWZgenPtDeborahWeight = qcdZTo2LWeight*ewkZDeborahsWeight*0.934

        elif search(r"^WtoLNu", sample, IGNORECASE):
            genParticles = Collection(event, "GenPart")
            GenV = list(
                [
                    gen
                    for gen in genParticles
                    if (abs(gen.pdgId) == 24) and gen.status == 22
                ]
            )
            if len(GenV) > 0:
                GenV_pt = GenV[0].pt
                ewkWWeight *= self.kFactorTool.getEWKW(GenV_pt)
                ewkWDeborahsWeight *= self.kFactorTool.getDeborahEWKW(GenV_pt)
                qcdWWeight *= self.kFactorTool.getQCDW(GenV_pt)
                qcdWWeightRenUp *= self.kFactorTool.getRenUpW(GenV_pt)
                qcdWWeightRenDown *= self.kFactorTool.getRenDownW(GenV_pt)
                qcdWWeightFacUp *= self.kFactorTool.getFacUpW(GenV_pt)
                qcdWWeightFacDown *= self.kFactorTool.getFacDownW(GenV_pt)
                WnnloWeight = 0.9135
                # if self.year == "2016APV":
                #     combinedWZgenPtWeight = qcdWWeight * ewkWWeight
                #     combinedWZgenPtDeborahWeight = (
                #         qcdWWeight * ewkWDeborahsWeight
                #     )  # VERIFY again..Deborah said QCD to be taken from python for 2016 but we have all same tunes unlike victor
                # elif self.year == "2016":
                #     combinedWZgenPtWeight = qcdWWeight * ewkWWeight
                #     combinedWZgenPtDeborahWeight = (
                #         qcdWWeight * ewkWDeborahsWeight
                #     )  # VERIFY again..Deborah said QCD to be taken from python for 2016 but we have all same tunes unlike victor
                # elif self.year == "2017":
                #     combinedWZgenPtWeight = qcdWWeight * ewkWWeight * 0.9135
                #     combinedWZgenPtDeborahWeight = (
                #         qcdWWeight * ewkWDeborahsWeight * 0.9135
                #     )
                if self.year == "2024":
                    combinedWZgenPtWeight = qcdWWeight * ewkWWeight * 0.9135
                    combinedWZgenPtDeborahWeight = (
                        qcdWWeight * ewkWDeborahsWeight * 0.9135
                    )

            combinedWZgenPtDylanWeight = 1.21
            # combinedWZgenPtWeight = qcdWWeight*ewkWWeight*0.9135
            # combinedWZgenPtDeborahWeight = qcdWWeight*ewkWDeborahsWeight*0.9135

        self.out.fillBranch("ewkWWeight", ewkWWeight)
        self.out.fillBranch("ewkZWeight", ewkZWeight)
        self.out.fillBranch("ewkWDeborahsWeight", ewkWDeborahsWeight)
        self.out.fillBranch("ewkZDeborahsWeight", ewkZDeborahsWeight)
        self.out.fillBranch("qcdWWeight", qcdWWeight)
        self.out.fillBranch("qcdZTo2LWeight", qcdZTo2LWeight)
        self.out.fillBranch("GenV_pt", GenV_pt)
        self.out.fillBranch("WnnloWeight", WnnloWeight)
        self.out.fillBranch("ZnnloWeight", ZnnloWeight)
        self.out.fillBranch("combinedWZgenPtWeight", combinedWZgenPtWeight)
        self.out.fillBranch(
            "combinedWZgenPtDeborahWeight", combinedWZgenPtDeborahWeight
        )
        self.out.fillBranch("combinedWZgenPtDylanWeight", combinedWZgenPtDylanWeight)

        if self.doSysVar:
            # Systematics - QCD Scale Factors
            self.out.fillBranch("qcdWWeightRenUp", qcdWWeightRenUp)
            self.out.fillBranch("qcdWWeightRenDown", qcdWWeightRenDown)
            self.out.fillBranch("qcdWWeightFacUp", qcdWWeightFacUp)
            self.out.fillBranch("qcdWWeightFacDown", qcdWWeightFacDown)
            self.out.fillBranch("qcdZTo2LWeightRenUp", qcdZTo2LWeightRenUp)
            self.out.fillBranch("qcdZTo2LWeightRenDown", qcdZTo2LWeightRenDown)
            self.out.fillBranch("qcdZTo2LWeightFacUp", qcdZTo2LWeightFacUp)
            self.out.fillBranch("qcdZTo2LWeightFacDown", qcdZTo2LWeightFacDown)

        return True