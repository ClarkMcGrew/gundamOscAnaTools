# Gundam Oscillation Analysis Toolkit

This is a set of tools for neutrino oscillation analyses using [GUNDAM (Generalized Unified Neutrino Data Analysis Method)](https://github.com/gundam-organization/gundam.git).  The primary (and initial) target for the tool kit is the DUNE experiment, but it aims to be applicable to any experiment.  The initial concentration are PNMS analyses.  

# Using the tool kit

The tool kit depends on the GUNDAM software found at (https://github.com/gundam-organization/gundam.git), and provides the applications specific tools needed for an oscillation analysis.

## Resources to calculate oscillation probabilities

Oscillation probabilities are provided to GUNDAM using the `Tabulated` and `Kriged` dial types.  The main implementation is the [TabulatedNuOscillator](./resources/TabulateNuOscillator/README.md) directory which uses NuOscillator by Dan Barrow (https://github.com/dbarrow257/NuOscillator.git) to calculate oscillation probabilities for both LBL and ATM neutrinos (the difference is that LBL oscillations use a single path-length, and density, while ATM oscillations depend on the zenith angle, and handle the PREM earth density model).
