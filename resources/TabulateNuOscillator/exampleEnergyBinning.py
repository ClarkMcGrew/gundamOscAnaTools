#!/usr/bin/python3
#
# An example script to build an 1D energy binning histogram.  This
# specific example is based around the DUNE LBL oscillation analysis,
# but is intended to provide an example starting point for any fix
# oscillation length analysis.
#
# The x axis of the TH1 defines the energy binning in GeV.

import sys
import os.path
import math
import argparse
import ROOT
import numpy

# Fill a list with energies for the bin lower edges. The last entry is
# the upper edge of the final bin.
def makeInverseEnergy(eBins,eMin,eMax,eRes):
    energies = []
    # ADJUST THIS COMMENT IF THE EXPECTED BINNING CHANGES: The fit is
    # expected to have about 20 energy bins roughly uniform in log(E), and
    # there should be at least three or four energy steps calculated per bin.
    # This leads to 80 energy grid points.
    minFraction = math.exp(-(math.log(eMax)-math.log(eMin))/80)
    # When the total number of energy samples is small, the fraction must
    # be bigger.  The fraction is picked to need slightly fewer limited steps
    # than energy grid points.
    maxFraction = math.exp(-1.5*(math.log(eMax)-math.log(eMin))/eBins)
    print("Minimum fraction is ", minFraction)
    targetFraction = 1-eRes
    fractionLimit = targetFraction
    fractionLimit = max(fractionLimit,minFraction)
    fractionLimit = min(fractionLimit,maxFraction)
    if (eRes < 1E-8): fractionLimit = 0.0
    print("Energy step",
          "target:",targetFraction,
          "min:",minFraction,
          "max:",maxFraction,
          "used:",fractionLimit)
    step = (1.0/eMin - 1.0/eMax)/eBins-1
    lastInvE = 1/eMax
    energies.append(1.0/lastInvE)
    while len(energies) <= eBins:
        invE = lastInvE + step
        if 1/invE < fractionLimit/lastInvE:
            invE = lastInvE/fractionLimit
            step = (1.0/eMin - invE)/(eBins-len(energies))
        energies.append(1.0/invE)
        lastInvE = invE
    energies.sort()
    return energies

parser = argparse.ArgumentParser()
parser.add_argument("file",help="Name of the output file")
parser.add_argument("-e","--energy-bins",type=int,default=1000,
                    help="Number of energy bins")
parser.add_argument("-m","--min-energy",type=float,default=0.1,
                    help="Minimum energy (GeV)")
parser.add_argument("-M","--max-energy",type=float,default=100,
                    help="Maximum energy (GeV)")
parser.add_argument("--energy-step",default="inverse",
                    help="Energy step (inverse or logarithmic)")
parser.add_argument("--energy-resolution",type=float,default=0.1,
                    help="The target limit on ratio between energy steps")
args = parser.parse_args()

output = ROOT.TFile(args.file,"recreate")

if "inv" in args.energy_step:
    energies = makeInverseEnergy(eBins = args.energy_bins,
                                 eMin = args.min_energy,
                                 eMax = args.max_energy,
                                 eRes = args.energy_resolution);
else:
    print("Invalid step type")
    sys.exit(1)

hist = ROOT.TH1D("energyBinning","Energy bins",
                 len(energies)-1,numpy.array(energies))

print("Energy bins",hist.GetXaxis().GetNbins(),
      "  Zenith bins", hist.GetYaxis().GetNbins(),
      "  Height bins", hist.GetZaxis().GetNbins())

output.Write()
sys.exit(0)
