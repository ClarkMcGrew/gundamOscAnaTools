// Fill a table of oscillation weights for the GUNDAM Tabulated dial type
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <cmath>

#include "TabulatedNuOscillator.hh"

// Add a "local" logging facility
#define LIB_NAME "TabulatedNuOscillator"
#ifndef LOUD_AND_PROUD
#define LOUD_AND_PROUD true
#endif

#ifndef LIB_COUT
#define LIB_COUT if (true) std::cout << LIB_NAME << " -- "
#endif
#ifndef LIB_CERR
#define LIB_CERR std::cerr << "ERROR: " << LIB_NAME << " -- "
#endif

// Define the global lookup tables.  These are normally hidden from outside
// code since they only exist in the shared library (the symbol isn't usually
// loaded).
TabulatedNuOscillator::ConfigLookup TabulatedNuOscillator::configLookup;
TabulatedNuOscillator::GlobalLookup TabulatedNuOscillator::globalLookup;

void TabulatedNuOscillator::FillInverseEnergyArray(
    std::vector<FLOAT_T>& energies, double eMin, double eMax, double eRes) {
    // ADJUST THIS COMMENT IF THE EXPECTED BINNING CHANGES: The fit is
    // expected to have about 20 energy bins roughly uniform in log(E), and
    // there should be at least three or four energy steps calculated per bin.
    // This leads to 80 energy grid points.
    double minFraction = std::exp(-(std::log(eMax)-std::log(eMin))/80);
    // When the total number of energy samples is to small, the fraction must
    // be bigger.  The fraction is picked to need slightly fewer limited steps
    // than energy grid points.
    double maxFraction = std::exp(-1.5*(std::log(eMax)-std::log(eMin))/energies.size());
    // The user wants steps closer than a particular energy resolution.
    double targetFraction = 1.0-eRes;
    // Choose the actual limiting fraciton.
    double fractionLimit = targetFraction;
    fractionLimit = std::max(fractionLimit,minFraction);
    fractionLimit = std::min(fractionLimit,maxFraction);
    if (eRes < 1E-8) fractionLimit = 0.0; // zero turns off limit
    LIB_COUT << "Energy step "
             << " target: " << targetFraction
             << " min: " << minFraction
             << " max: " << maxFraction
             << " used: " << fractionLimit
             << std::endl;
    double step = (1.0/eMin - 1.0/eMax)/(energies.size()-1);
    std::size_t bin = 0;
    double lastInvE = 1.0/eMax;
    energies[bin++] = 1.0/lastInvE;
    while (bin < energies.size()) {
        double invE = lastInvE + step;
        if (1.0/invE < fractionLimit/lastInvE) {
            invE = lastInvE/fractionLimit;
            if (bin+1 < energies.size()) {
                step = (1.0/eMin - invE)/(energies.size() - bin);
            }
        }
        energies[bin++] = 1.0/invE;
        lastInvE = invE;
    };
    if (energies[bin-1]>eMin) energies[bin-1] = eMin;
}

void TabulatedNuOscillator::FillLogarithmicEnergyArray(
    std::vector<FLOAT_T>& energies, double eMin, double eMax) {
    double eMaxLog = std::log(eMax);
    double eMinLog = std::log(eMin);
    double step=(eMaxLog - eMinLog)/(energies.size()-1);
    std::size_t bin = 0;
    do {
        double v = eMinLog + step*bin;
        energies[bin++] = std::exp(v);
    } while (bin < energies.size());
}

void TabulatedNuOscillator::FillEnergyArray(
    std::vector<FLOAT_T>& energies, const std::string& type,
    double eMin, double eMax, double eRes) {
    if (type.find("inv") != std::string::npos) {
        FillInverseEnergyArray(energies,eMin,eMax,0.10);
        // FillInverseEnergyArray(energies,eMin,eMax,0.05);
    }
    else {
        LIB_COUT << "WARNING -- Logarithmic energy step loses precision"
                 << std::endl;
        LIB_COUT << "WARNING -- Use inverse, instead logarithmic energy step"
                 << std::endl;
        FillLogarithmicEnergyArray(energies,eMin,eMax);
    }

    // NuOscillator needs the energies in increasing order.
    std::sort(energies.begin(), energies.end());
    energies[0] = eMin;
    energies[energies.size()-1] = eMax;

#ifdef DUMP_ENERGY
    for (int i = 0; i < energies.size(); ++i) {
        std::cout << "E " << i
                  << " " << energies[i]
                  << " " << std::log10(energies[i])
                  << std::endl;
    }
#endif

}

void TabulatedNuOscillator::FillZenithArray(std::vector<FLOAT_T>& zenith) {
    if (zenith.size() < 1) return;
    if (zenith.size() < 2) {
        zenith[0] = -1.0;
        return;
    }
    double minCos = -1.0;
    double maxCos = 1.0;
#ifdef DEBUG_LIMIT_ZENITH
#warning Limit zenith angle range
    maxCos = -0.5;
#endif
    double minPath = RoughZenithPath(maxCos);
    double maxPath = RoughZenithPath(minCos);
    double step = (maxPath - minPath)/(zenith.size() - 1);
    double maxCosZStep = 0.05;  // 40 bins, but could be a parameter.
    double path = minPath;
    double lastC = maxCos;
    int bin = 0;
    zenith[bin++] = lastC;
    for (double c = zenith[0]; bin < zenith.size(); c -= 1E-8) {
        double p = RoughZenithPath(c);
        if (lastC - c < maxCosZStep &&  p - path < step) continue;
        else if (p-path < step) {
            step = (maxPath - p) / (zenith.size() - bin - 1);
        }
        zenith[bin++] = c;
        lastC = c;
        path = p;
        if (bin >= zenith.size()) break;
    }
    std::sort(zenith.begin(), zenith.end());
    zenith[0] = minCos;
    zenith[zenith.size()-1] = maxCos;
}

double TabulatedNuOscillator::RoughZenithPath(double cosz) {
    const double Rd{6371}; //Average Earth Radius in km (average)
    const double Rp{Rd + 20.0}; // Very rough production elevation.
    if (cosz > 1.0) cosz = 1.0;
    if (cosz < -1.0) cosz = -1.0;
    double L = std::sqrt(Rd*Rd*(cosz*cosz-1.0) + Rp*Rp) - Rd*cosz;
    return L;
}

bool TabulatedNuOscillator::AlmostEqual(double a, double b) {
    const double diff = std::abs(a-b);
    const double avg = std::abs(a) + std::abs(b);
    const double delta = 1E-10;
    if (diff > delta*avg+delta) return false;
    return true;
}

void TabulatedNuOscillator::ConfigureNuOscillator(const TableGlobals& globals) {
    if (globals.nuOscillatorConfig.empty()) return;
    if (globals.oscEnergies.empty()) return;

    ConfigLookup::iterator config
        = configLookup.find(globals.nuOscillatorConfig);
    if (config != configLookup.end()) {
        // Already initialized, check it's the same.
        if (globals.oscEnergies.size() != config->second.energies.size()) {
            LIB_CERR << "Only one energy binning is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for (std::size_t i = 0; i < globals.oscEnergies.size(); ++i) {
            if (AlmostEqual(globals.oscEnergies[i],
                            config->second.energies[i])) {
                continue;
            }
            LIB_CERR << "Only one energy binning is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (globals.oscZenith.size() != config->second.zenith.size()) {
            LIB_CERR << "Only one zenith binning is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        for (std::size_t i = 0; i < globals.oscZenith.size(); ++i) {
            if (AlmostEqual(globals.oscZenith[i],
                            config->second.zenith[i])) {
                continue;
            }
            LIB_CERR << "Only one zenith binning is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (not AlmostEqual(globals.oscPath,
                            config->second.oscPath)) {
            LIB_CERR << "Only one path length is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (not AlmostEqual(globals.oscDensity,
                            config->second.oscDensity)) {
            LIB_CERR << "Only one density is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (not AlmostEqual(globals.oscElectronDensity,
                            config->second.oscElectronDensity)) {
            LIB_CERR << "Only one electron density is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (not AlmostEqual(globals.oscProdHeight,
                            config->second.oscProdHeight)) {
            LIB_CERR << "Only one production height is allowed for "
                     << globals.nuOscillatorConfig
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }

        return;
    }

    NuOscillatorConfig& newConfig
        = configLookup[globals.nuOscillatorConfig];
    newConfig.name = globals.nuOscillatorConfig;
    newConfig.energies = globals.oscEnergies;
    newConfig.zenith = globals.oscZenith;
    newConfig.oscParIndex = globals.oscParIndex;

    std::unique_ptr<OscillatorFactory> factory
        = std::make_unique<OscillatorFactory>();
#ifdef TABULATED_NUOSCILLATOR_DECONSTRUCTABLE
    newConfig.oscillator.reset(factory->CreateOscillator(newConfig.name));
#else
    newConfig.oscillator = factory->CreateOscillator(newConfig.name);
#endif

    if (newConfig.oscillator->ReturnImplementationName().find("Unbinned_")
        == std::string::npos) {
        LIB_CERR << "NuOscillator must use an unbinned configuration"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }

    newConfig.oscillator->SetEnergyArrayInCalcer(newConfig.energies);
    if (not newConfig.oscillator->CosineZIgnored()) {
        if (newConfig.zenith.empty()) {
            LIB_CERR << "Zenith angle bins not filled for atmospheric oscillator"
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        newConfig.oscillator->SetCosineZArrayInCalcer(newConfig.zenith);
    }
    else if (not newConfig.zenith.empty()) {
        LIB_CERR << "Zenith angle filled for LBL oscillator"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }
    newConfig.oscPath = globals.oscPath;
    newConfig.oscDensity = globals.oscDensity;
    newConfig.oscElectronDensity = globals.oscElectronDensity;
    newConfig.oscProdHeight = globals.oscProdHeight;

    newConfig.oscillator->Setup();
    newConfig.oscParams.resize(newConfig.oscillator->ReturnNOscParams());

    LIB_COUT << "Configured: " << newConfig.name << std::endl;
}

// Provide the initializetable entry point required by the GUNDAM tabulated
// dials.  This requires string arguments:
//
// CONFIG <file-name>
// FLUX_FLAVOR [anti-]{electron,muon,tau}
// INTERACTION_FLAVOR [anti-]{electron,muon,tau}
// PARAMETERS <list of SS12, SS23, SS13, DM21, DM32, DCP>
// ENERGY_BINS <integer> -- Number of energy bins for each neutrino type
// MIN_ENERGY <double> -- Minimum neutrino energy in GeV
// MAX_ENERGY <double> -- Maximum neutrino energy in GeV
// ENERGY_STEP {inverse,logarithmic} -- The energy binning to use (def: inverse)
// ENERGY_SMOOTH <double> -- The 1/E (1/GeV) smoothing (dev 0.1, limits bins considered).
// ENERGY_RESOLUTION <double> -- Fractional energy resolution to smooth over (def: 0.02)
// ZENITH_BINS <integer>  -- Number of zenith cosine bins (def: 0)
// ZENITH_SMOOTH <integer> -- The pathlength (km) smoothing (def: 100, limits bins considered).
// ZENITH_RESOLUTION <double> -- Angle (radian) to smooth over at horizon (def: 0.05).
// DENSITY <double>    -- Density in gm/cc
// ELECTRON_DENSITY <double> -- Almost always 0.5
// PATH <double>       -- Path length in kilometers (for LBL)
// PRODUCTION_HEIGHT <double>  -- Neutrino production height in km (for ATM)
//
extern "C"
int initializeTable(const char* name, int argc, const char* argv[],
                    int suggestedBins) {
    LIB_COUT << "Initialize: " << name << std::endl;
    TabulatedNuOscillator::TableGlobals& globals
        = TabulatedNuOscillator::globalLookup[name];

    // Set default values.
    globals.name = name;
    globals.oscEnergyBins = 1000;
    globals.oscZenithBins = 0;
    globals.oscPath = 1300.0; // km
    globals.oscProdHeight = 17.0; // km
    globals.oscMinEnergy = 0.050; // GeV
    globals.oscDensity = 2.6; // gm/cc
    globals.oscElectronDensity = 0.5;
    globals.oscEnergyStep = "inverse";
    globals.oscEnergySmooth = 0.4; // 1/GeV
    globals.oscEnergyResol = 0.05; // relative
    globals.oscZenithSmooth = 800.0; // km
    globals.oscZenithResol = 0.0; // radian

    for (int i = 0; i < argc; ++i) {
        LIB_COUT << "Argument: " << argv[i] << std::endl;
        globals.arguments.emplace_back(argv[i]);
    }

    // Get the configuration file name.  GUNDAM will have already expanded
    // any environment variables.
    for (std::string arg: globals.arguments) {
        if (arg.find("CONFIG") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.nuOscillatorConfig;
    }

    // Get the flavor of the flux (i.e. the parent neutrino)
    for (std::string arg: globals.arguments) {
        if (arg.find("FLUX_FLAVOR") != 0) continue;
        std::istringstream tmp(arg);
        std::string flavor;
        tmp >> arg >> flavor;
        int nuType = 1;
        if (flavor.find("anti") != std::string::npos) {
            nuType = -1;
            flavor = flavor.substr(flavor.find("anti-")+5);
        }
        globals.oscInitialFlavor = nuType*NeutrinoFlavour_StrToInt(flavor);
        break;
    }

    // Get the flavor of the interaction type (i.e. the final neutrino)
    for (std::string arg: globals.arguments) {
        if (arg.find("INTERACTION_FLAVOR") != 0) continue;
        std::istringstream tmp(arg);
        std::string flavor;
        tmp >> arg >> flavor;
        int nuType = 1;
        if (flavor.find("anti") != std::string::npos) {
            nuType = -1;
            flavor = flavor.substr(flavor.find("anti-")+5);
        }
        globals.oscFinalFlavor = nuType*NeutrinoFlavour_StrToInt(flavor);
        break;
    }

    // Get order of the oscillation parameters from the fit.
    for (std::string arg: globals.arguments) {
        if (arg.find("PARAMETERS") != 0) continue;
        arg.erase(arg.find("PARAMETERS"),10);
        arg.erase(std::remove_if(
                      arg.begin(), arg.end(),
                      [](unsigned char c){return std::isspace(c);}),
                  arg.end());
        std::istringstream tmp(arg);
        std::string param;
        int index = 0;
        while (std::getline(tmp, param, ',')) {
            if      (param == "SS12") {globals.oscParIndex.ss12 = index++;}
            else if (param == "SS23") {globals.oscParIndex.ss23 = index++;}
            else if (param == "SS13") {globals.oscParIndex.ss13 = index++;}
            else if (param == "DM21") {globals.oscParIndex.dm21 = index++;}
            else if (param == "DM32") {globals.oscParIndex.dm32 = index++;}
            else if (param == "DCP")  {globals.oscParIndex.dcp = index++;}
            else {
                LIB_CERR << "Unknown name parameter: " << param << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
        break;
    }

    // Get the number of oscillation bins per neutrino
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_BINS") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergyBins;
        break;
    }

    // Get the minimum energy for neutrinos
    for (std::string arg: globals.arguments) {
        if (arg.find("MIN_ENERGY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscMinEnergy;
        break;
    }

    // Get the maximum energy for neutrinos
    for (std::string arg: globals.arguments) {
        if (arg.find("MAX_ENERGY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscMaxEnergy;
        break;
    }

    // Get the type of step for the energies
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_STEP") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergyStep;
        break;
    }


    // Get the number of bins to smooth across for energies.
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_SMOOTH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergySmooth;
        break;
    }

    // Get the fractional energy resolution for energy to smooth over.
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_RESOLUTION") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergyResol;
        break;
    }

    // Get the number of oscillation bins per neutrino
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_BINS") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithBins;
        break;
    }

    // Get the number of bins to smooth across for zenith angle
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_SMOOTH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithSmooth;
        break;
    }

    // Get the angle (radians) to some over zenith angle.
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_RESOLUTION") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithResol;
        break;
    }

    // The material density along the path length.  This only works for LBL
    // when the beam doesn't go deep into the earth.  DUNE goes 30 km deep, so
    // assume the density is constant.
    for (std::string arg: globals.arguments) {
        if (arg.find("DENSITY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscDensity;
        break;
    }

    // The material electron density along the path length.  This only works
    // for LBL when the beam doesn't go deep into the earth.  DUNE goes 30 km
    // deep, so assume the density is constant.
    for (std::string arg: globals.arguments) {
        if (arg.find("ELECTRON_DENSITY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscElectronDensity;
        break;
    }

    // The path length.
    for (std::string arg: globals.arguments) {
        if (arg.find("PATH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscPath;
        break;
    }

    // The path length.
    for (std::string arg: globals.arguments) {
        if (arg.find("PRODUCTION_HEIGHT") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscProdHeight;
        break;
    }

    // It's brutal, but stop if the user forgot to configure the oscillator.
    if (globals.nuOscillatorConfig.empty()) {
        LIB_COUT << "The NuOscillator config file must be provided"
                 << std::endl;
        LIB_CERR << "The NuOscillator config file must be provided"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Set the upper limit for the L over E.
    globals.oscMaxLoverE = globals.oscPath/globals.oscMinEnergy;

    if (globals.oscEnergyBins < 2) {
        LIB_COUT << "The NuOscillator must have at least two energy bins"
                 << std::endl;
        LIB_CERR << "The NuOscillator must have at least two energy bins"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }
    globals.oscEnergies.resize(globals.oscEnergyBins);
    globals.oscZenith.resize(globals.oscZenithBins);

    // Make sure the upper energy limit is sane.
    if (globals.oscMaxEnergy < globals.oscMinEnergy) {
        double minLoverE = globals.oscMaxLoverE/globals.oscEnergies.size();
        globals.oscMaxEnergy = globals.oscPath/minLoverE;
    }

    TabulatedNuOscillator::FillEnergyArray(globals.oscEnergies,
                                           globals.oscEnergyStep,
                                           globals.oscMinEnergy,
                                           globals.oscMaxEnergy,
                                           globals.oscEnergyResol);

    TabulatedNuOscillator::FillZenithArray(globals.oscZenith);

#ifdef DEBUG_ENERGY_BINNING
    // Print the energy binning
    for (std::size_t bin = 0; bin < globals.oscEnergies.size(); ++bin) {
        double roughDMS = 2.5E-3;
        double LoverE = globals.oscPath/globals.oscEnergies[bin];
        const int defPrec = std::cout.precision();
        LIB_COUT << bin << " Approx Phase: "
                 << std::setprecision(3) << std::fixed
                 << 1.27*roughDMS*LoverE/3.14 << "*pi"
                 << std::setprecision(1)
                 << " L/E: " << LoverE
                 << std::setprecision(defPrec) << std::defaultfloat
                 << " E: " << globals.oscEnergies[bin]
                 << std::endl;
    }
#endif

    ConfigureNuOscillator(globals);
    TabulatedNuOscillator::NuOscillatorConfig& config = TabulatedNuOscillator::configLookup[globals.nuOscillatorConfig];

    int zSmooth = 0;
    int eSmooth = 0;

    int zenithIndex = 0;
    // The zenith loops are done like this since oscZenith.size() might be
    // zero.
    do {
        int iz = std::max(0,zenithIndex-zSmooth);
        do {
            double zenith = -999.0;
            if (iz < globals.oscZenith.size()) {
                zenith = globals.oscZenith[iz];
            }
            for (int energyIndex = 0;
                 energyIndex < (int) globals.oscEnergies.size();
                 ++energyIndex) {
                std::size_t bin
                    = zenithIndex*globals.oscEnergies.size() + energyIndex;
                for (int ie = std::max(0, energyIndex - eSmooth);
                     ie < std::min((int) globals.oscEnergies.size(),
                                   energyIndex + eSmooth + 1);
                     ++ie) {
                    double energy = globals.oscEnergies[energyIndex];
                    const FLOAT_T* address
                        = config.oscillator->ReturnWeightPointer(
                            globals.oscInitialFlavor,globals.oscFinalFlavor,
                            energy, zenith);
                    globals.weightAddress.emplace_back(
                        TabulatedNuOscillator::TableGlobals::OscWeight(
                            {bin, address, 1.0}));
                }
            }
        } while (++iz < std::min((int) globals.oscZenith.size(),
                                 zenithIndex+zSmooth+1));
    } while (++zenithIndex < globals.oscZenith.size());

    // Find the maximum bin in the table
    std::size_t bins = 0;
    for (TabulatedNuOscillator::TableGlobals::OscWeight &weight
             : globals.weightAddress) {
        bins = std::max(bins,weight.index+1);
    }

    // Find the sum of the weights for a particular bin.
    std::vector<double> work(bins);
    for (TabulatedNuOscillator::TableGlobals::OscWeight &weight
             : globals.weightAddress) {
        work[weight.index] += weight.weight;
    }

    // Rescale the weights so they sum to one
    for (TabulatedNuOscillator::TableGlobals::OscWeight &weight
             : globals.weightAddress) {
        if (work[weight.index] > 0.0) weight.weight /= work[weight.index];
    }

    LIB_COUT << "Oscillation size: " << globals.weightAddress.size()
             << " Table size: " << bins
             << " (suggested size: " << suggestedBins << ")" << std::endl;

    return bins;
}

// Calculate the approximate "delta" along the energy axis.  The table spacing
// is approximately 1/E while the bins are labeled by E.  This is the
// difference in 1/E.  This returns the absolute value of the change.
double energyBinDelta(double e2, double e1) {
    double v = 1/e2 - 1/e1;
    return std::abs(v);
}

// Calculate the approximate "delta" along the zenith angle axis.  The table
// spacing is approximately by path length while the bins are labeled in
// cos(zenithAngle).  This is the approximate difference in path length. This
// returns the absolute value of the change.
double zenithBinDelta(double c2, double c1) {
    double v = TabulatedNuOscillator::RoughZenithPath(c2)
        - TabulatedNuOscillator::RoughZenithPath(c1);
    return std::abs(v);
}

#define IDX(e,z) ((z)*globals.oscEnergies.size() + (e))

// Add a weight to the output list.  The bin in the table is ie,iz (energy,
// zenith).  The weight is calculated from smooth*window*area, and it is added
// if smooth*window more than zero.
int AddWeight(TabulatedNuOscillator::TableGlobals& globals,
              int ie, int iz,
              double smooth, double window, double area,
              int entry, int index[], double weights[]) {
    if (0 <= ie and ie < globals.oscEnergies.size()
        and 0 <= iz and iz < globals.oscZenith.size()) {
        index[entry] = IDX(ie, iz);
        weights[entry] = smooth;
        weights[entry] *= window;
        if (weights[entry] > 0.0) {
            weights[entry] *= area;
            ++entry;
        }
    }
    return entry;
}


// Provide the weightTable entry point required by the GUNDAM tabulated dials.
// The `index[]` and `weights[]` arrays must have at least `entries` elements
// allocated.  The function returns the number of entries filled in the
// `index[]` and `weights[]` arrays.
extern "C"
int weightTable(const char* name, int bins,
                int varc, double varv[],
                int entries, int index[], double weights[]) {
    double energyValue = varv[0];
    double zenithValue = -1.0;
    if (varc>1) {
        zenithValue = varv[1];
        double zenithPath
            = TabulatedNuOscillator::RoughZenithPath(zenithValue);
    }
    TabulatedNuOscillator::TableGlobals& globals
        = TabulatedNuOscillator::globalLookup[name];

    // Find the index of the entry oscZenith table that is greater than or
    // equal to the zenithValue. The index will be zero if the zenith angle
    // table doesn't exist.
    int zenithIndex = 0;
    if (globals.oscZenith.size() > 1) {
        // double value = zenithValue - 1.0/globals.oscZenith.size();
        double value = zenithValue;
        auto zenithIt = std::lower_bound(globals.oscZenith.begin(),
                                         globals.oscZenith.end(),
                                         value);
        if (value < globals.oscZenith.front()) {
            zenithIt = globals.oscZenith.begin();
        }
        else if ( value >= globals.oscZenith.back()) {
            zenithIt = globals.oscZenith.end()-1;
        }
        zenithIndex = std::distance(globals.oscZenith.begin(), zenithIt);
    }
    if (0 < zenithIndex) {
        if (zenithValue < globals.oscZenith[zenithIndex-1]
            or globals.oscZenith[zenithIndex] < zenithValue) {
#ifdef DUMP_Z_PROBLEM
            std::cout << "ZenithIndex " << zenithIndex
                      << " " << globals.oscZenith[zenithIndex-1]
                      << " " << zenithValue
                      << " " << globals.oscZenith[zenithIndex]
                      << std::endl;
#endif
        }
    }

    // Find the index for the energy in the oscEnergies table that is greater
    // or equal to the energyValue.
    auto energyIt = std::lower_bound(globals.oscEnergies.begin(),
                                     globals.oscEnergies.end(),
                                     energyValue);
    if (energyValue < (globals.oscEnergies.front())) {
        energyIt = globals.oscEnergies.begin();
    }
    else if (energyValue >= globals.oscEnergies.back()) {
        energyIt = globals.oscEnergies.end()-1;
    }
    int energyIndex = std::distance(globals.oscEnergies.begin(), (energyIt));
    if (0<energyIndex) {
        if (energyValue < globals.oscEnergies[energyIndex-1]
            or globals.oscEnergies[energyIndex] < energyValue) {
#ifdef DUMP_E_PROBLEM
            std::cout << "EnergyIndex " << energyIndex
                      << " " << globals.oscEnergies[energyIndex-1]
                      << " " << energyValue
                      << " " << globals.oscEnergies[energyIndex]
                      << std::endl;
#endif
        }
    }

    // Check that the table will be the right size.
    {
        int requiredBins = std::max(std::size_t(1),globals.oscZenith.size());
        requiredBins *= globals.oscEnergies.size();
        if (requiredBins != bins) {
            LIB_CERR << "Table is the wrong size: " << bins
                     << " != " << requiredBins << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Check that there is enough space.
    if (entries < 4) {
        LIB_CERR << "Not enough entries in weight arrays" << std::endl;
        std::exit(EXIT_FAILURE);
    }
    int entry = 0;

    // Find the local differential (and bin spacing) based on the size of the
    // "central" bin
    double energyDelta = 1.0;
    double energyDifferential = 1.0;
    {
        double lower;
        double higher;
        if (0 < energyIndex and energyIndex < globals.oscEnergies.size()) {
            lower = globals.oscEnergies[energyIndex-1];
            higher = globals.oscEnergies[energyIndex];
        }
        else if (energyIndex<1) {
            lower = globals.oscEnergies[0];
            higher = globals.oscEnergies[1];
        }
        if (globals.oscEnergies.size() <= energyIndex) {
            lower = globals.oscEnergies[globals.oscEnergies.size()-2];
            higher = globals.oscEnergies[globals.oscEnergies.size()-1];
        }
        energyDelta = higher - lower;
        energyDifferential = energyBinDelta(higher,lower);
    }

    double zenithDelta = 1.0;
    double zenithDifferential = 1.0;
    {
        double lower;
        double higher;
        if (0 < zenithIndex and zenithIndex < globals.oscZenith.size()) {
            lower = globals.oscZenith[zenithIndex-1];
            higher = globals.oscZenith[zenithIndex];
        }
        else if (zenithIndex<1) {
            lower = globals.oscZenith[0];
            higher = globals.oscZenith[1];
        }
        if (globals.oscZenith.size() <= zenithIndex) {
            lower = globals.oscZenith[globals.oscZenith.size()-2];
            higher = globals.oscZenith[globals.oscZenith.size()-1];
        }
        zenithDelta = higher - lower;
        zenithDifferential = zenithBinDelta(higher,lower);
    }

    // Average over a window
    const double energyBinSigma = globals.oscEnergySmooth;
    const double pathSigma = globals.oscZenithSmooth;

    // Smooth over the resolution.  The energy resolution is relative, and the
    // zenith resolution is in radians.
    const double sigmaE = globals.oscEnergyResol;
    const double sigmaZ = globals.oscZenithResol;

    const double windowThres = 0.1; // about +/- 2 sigma around center
    for (int ie = 0; ie < globals.oscEnergies.size(); ++ie) {
        // The lowE is the index of lower energy bin.  The location of the
        // bin will always be less than or equal to energyValue.  The lowE
        // value may be less than zero (which means the index should be
        // ignored).
        int lowE = energyIndex - ie - 1;

        // The highE is the index of the upper zenith bin.  The location of
        // the bin will always be greater than or equal to the energyValue The
        // highE value may be greater than or equal to the size of the
        // oscZenith vector (which means the index should be ignored).
        int highE = energyIndex + ie;

        // Find the weight for the averaging window.  This will be ignored for
        // the central box.
        double lowerBinEnergy = 0.0;
        if (0 <= lowE) {
            if (ie == 0) {
                if (highE < globals.oscEnergies.size()) {
                    double delta = energyBinDelta(energyValue,
                                                  globals.oscEnergies[lowE]);
                    double diff = energyBinDelta(globals.oscEnergies[highE],
                                                 globals.oscEnergies[lowE]);
                    lowerBinEnergy = 1.0 - delta/diff;
                }
                else lowerBinEnergy = 1.0;
            }
            else if (energyBinSigma > 2.0*energyDifferential) {
                double deltaInvE = energyBinDelta(globals.oscEnergies[lowE],
                                                  energyValue);
                deltaInvE /= energyBinSigma;
                lowerBinEnergy = std::exp(-0.5*deltaInvE*deltaInvE);
            }
        }

        double upperBinEnergy = 0.0;
        if (highE < globals.oscEnergies.size()) {
            if (ie == 0) {
                if (0 <= lowE) {
                    double delta = energyBinDelta(globals.oscEnergies[highE],
                                                  energyValue);
                    double diff = energyBinDelta(globals.oscEnergies[highE],
                                                 globals.oscEnergies[lowE]);
                    upperBinEnergy = 1.0 - delta/diff;
                }
                else upperBinEnergy = 1.0;
            }
            else if (energyBinSigma > 2.0*energyDifferential) {
                double deltaInvE = energyBinDelta(globals.oscEnergies[highE]
                                                  ,energyValue);
                deltaInvE /= energyBinSigma;
                upperBinEnergy = std::exp(-0.5*deltaInvE*deltaInvE);
            }
        }

        // Find the weight for the energy resolution smoothing.  This will be
        // ignored for the central box.
        double lowerEnergy = 0.0;
        if (ie == 0) lowerEnergy = 1.0;
        else if (0 <= lowE
                 and sigmaE*energyValue > energyDelta) {
            double binValue = globals.oscEnergies[lowE];
            // Smooth of the neutrino energy resolution
            double deltaE = binValue - energyValue;
            deltaE /= energyValue*sigmaE;
            lowerEnergy = std::exp(-0.5*deltaE*deltaE);
        }

        double upperEnergy = 0.0;
        if (ie == 0) upperEnergy = 1.0;
        else if (highE < globals.oscEnergies.size()
                 and sigmaE*energyValue > energyDelta) {
            double binValue = globals.oscEnergies[highE];
            double deltaE = binValue - energyValue;
            deltaE /= energyValue*sigmaE;
            upperEnergy = std::exp(-0.5*deltaE*deltaE);
        }

        int iz = 0;
        do {

            // The lowZ is the index of lower zenith bin.  The location of the
            // bin will always be less than or equal to zenithValue.  The lowZ
            // value may be less than zero (which means the index should be
            // ignored).
            int lowZ = zenithIndex-iz - 1;

            // The highZ is the index of the upper zenith bin.  The location
            // of the bin will always be greater than or equal to the
            // zenithValue The highZ value may be greater than or equal to the
            // size of the oscZenith vector (which means the index should be
            // ignored).
            int highZ = zenithIndex + iz;

            double lowerBinZenith = 0.0;
            if (0 <= lowZ) {
                if (iz == 0) {
                    if (highZ < globals.oscZenith.size()) {
                        double delta = zenithBinDelta(globals.oscZenith[lowZ],
                                                      zenithValue);
                        double diff = zenithBinDelta(globals.oscZenith[highZ],
                                                     globals.oscZenith[lowZ]);
                        lowerBinZenith = 1.0 - delta/diff;
                    }
                    else lowerBinZenith = 1.0;
                }
                else if (pathSigma > 2.0*zenithDifferential) {
                    double deltaPath = zenithBinDelta(globals.oscZenith[lowZ],
                                                      zenithValue);
                    deltaPath /= pathSigma;
                    lowerBinZenith = std::exp(-0.5*deltaPath*deltaPath);
                }
            }

            double upperBinZenith = 0.0;
            if (highZ < globals.oscZenith.size()) {
                if (iz == 0) {
                    if (0 <= lowZ) {
                        double delta = zenithBinDelta(globals.oscZenith[highZ],
                                                      zenithValue);
                        double diff = zenithBinDelta(globals.oscZenith[highZ],
                                                     globals.oscZenith[lowZ]);
                        upperBinZenith = 1.0 - delta/diff;
                    }
                    else upperBinZenith = 1.0;
                }
                else if (pathSigma > 2.0*zenithDifferential) {
                    double deltaPath = zenithBinDelta(globals.oscZenith[highZ],
                                                      zenithValue);
                    deltaPath /= pathSigma;
                    upperBinZenith = std::exp(-0.5*deltaPath*deltaPath);
                }
            }
            if (upperBinZenith<windowThres && lowerBinZenith<windowThres) break;

            double lowerZenith = 0.0;
            if (iz == 0) lowerZenith = 1.0;
            else if (0 <= lowZ and sigmaZ > zenithDelta) {
                double binValue = globals.oscZenith[lowZ];
#ifdef ANGLE_SMOOTHING
                // Smooth over angle since it's the direction that is
                // uncertain
                double deltaZ = std::acos(binValue) - std::acos(zenithValue);
#else
                // Smooth over angle since it's the direction that is
                // uncertain
                double deltaZ = binValue - zenithValue;
#endif
                deltaZ /= sigmaZ;
                lowerZenith = std::exp(-0.5*deltaZ*deltaZ);
            }

            double upperZenith = 0.0;
            if (iz == 0) upperZenith = 1.0;
            else if (highZ < globals.oscZenith.size() and sigmaZ > zenithDelta) {
                double binValue = globals.oscZenith[highZ];
#ifdef ANGLE_SMOOTHING
                // Smooth over angle since it's the direction that is
                // uncertain
                double deltaZ = std::acos(binValue) - std::acos(zenithValue);
#else
                // Smooth over angle since it's the direction that is
                // uncertain
                double deltaZ = binValue - zenithValue;
#endif
                deltaZ /= sigmaZ;
                upperZenith = std::exp(-0.5*deltaZ*deltaZ);
            }
            if (upperZenith<windowThres and lowerZenith<windowThres) break;

            // Calculate the area correction around the high and low points.
            double dUpperE = 0.0;
            if (0 <= highE and highE < globals.oscEnergies.size()-1) {
                dUpperE += 0.5*energyBinDelta(globals.oscEnergies[highE],
                                              globals.oscEnergies[highE+1]);
            }
            if (0 < highE and highE < globals.oscEnergies.size()) {
                dUpperE += 0.5*energyBinDelta(globals.oscEnergies[highE-1],
                                              globals.oscEnergies[highE]);
            }

            double dLowerE = 0.0;
            if (0 <= lowE and lowE < globals.oscEnergies.size()-1) {
                dLowerE += 0.5*energyBinDelta(globals.oscEnergies[lowE],
                                              globals.oscEnergies[lowE+1]);
            }
            if (0 < lowE and lowE < globals.oscEnergies.size()) {
                dLowerE += 0.5*energyBinDelta(globals.oscEnergies[lowE-1],
                                              globals.oscEnergies[lowE]);
            }

            // Check if the area is for the central bin and override the
            // area correction if it is.
            if (ie == 0) {
                if (0 <= lowE and highE < globals.oscEnergies.size()) {
                    dLowerE = energyBinDelta(globals.oscEnergies[highE],
                                             globals.oscEnergies[lowE]);
                }
                else if (lowE < 0) {
                    dLowerE = energyBinDelta(globals.oscEnergies[1],
                                             globals.oscEnergies[0]);
                }
                else if (globals.oscEnergies.size() <= highE) {
                    dLowerE = energyBinDelta(
                        globals.oscEnergies[globals.oscEnergies.size()-1],
                        globals.oscEnergies[globals.oscEnergies.size()-2]);
                }
                else dLowerE = 0.0;
                dLowerE *= 4.0;
                dUpperE = dLowerE;

            }

            double dUpperZ = 0.0;
            if (0 <= highZ and highZ < globals.oscZenith.size()-1) {
                dUpperZ += 0.5*zenithBinDelta(globals.oscZenith[highZ],
                                              globals.oscZenith[highZ+1]);
            }
            if (0 < highZ and highZ < globals.oscZenith.size()) {
                dUpperZ += 0.5*zenithBinDelta(globals.oscZenith[highZ-1],
                                              globals.oscZenith[highZ]);
            }

            double dLowerZ = 0.0;
            if (0 <= lowZ and lowZ < globals.oscZenith.size()-1) {
                dLowerZ += 0.5*zenithBinDelta(globals.oscZenith[lowZ],
                                              globals.oscZenith[lowZ+1]);
            }
            if (0 < lowZ and lowZ < globals.oscZenith.size()) {
                dLowerZ += 0.5*zenithBinDelta(globals.oscZenith[lowZ-1],
                                              globals.oscZenith[lowZ]);
            }


            // Check if the area is for the central bin and override the
            // area correction if it is.
            if (iz == 0) {
                if (0 <= lowZ and highZ < globals.oscZenith.size()) {
                    dLowerZ = zenithBinDelta(globals.oscZenith[highZ],
                                             globals.oscZenith[lowZ]);
                }
                else if (lowZ < 0) {
                    dLowerZ = zenithBinDelta(globals.oscZenith[1],
                                             globals.oscZenith[0]);
                }
                else if (globals.oscZenith.size() <= highZ) {
                    dLowerZ = zenithBinDelta(
                        globals.oscZenith[globals.oscZenith.size()-1],
                        globals.oscZenith[globals.oscZenith.size()-2]);
                }
                else dLowerZ = 0.0;
                dLowerZ *= 4.0;
                dUpperZ = dLowerZ;
            }

            if (dUpperE < 0 or dLowerE<0 or dUpperZ < 0 or dLowerZ < 0) {
                LIB_COUT << "smoothing area"
                         << " dUpperE " << dUpperE
                         << " dLowerE " << dLowerE
                         << " dUpperZ " << dUpperZ
                         << " dLowerZ " << dLowerZ
                         << " high Z " << highZ
                         << " low Z " << lowZ
                         << std::endl;
                LIB_CERR << "smoothing area wrong"
                         << std::endl;
                std::exit(1);
            }

            int lastEntry = entry;
            if (entries < entry+4) {
                LIB_CERR << "Weight table too small: " << entries
                        << std::endl;
                break; // Not enough space for more points
            }

            entry = AddWeight(globals,highE,highZ,
                              upperEnergy*upperZenith,
                              upperBinEnergy*upperBinZenith,
                              dUpperE*dUpperZ,
                              entry,index,weights);

            entry = AddWeight(globals,lowE,lowZ,
                              lowerEnergy*lowerZenith,
                              lowerBinEnergy*lowerBinZenith,
                              dLowerE*dLowerZ,
                              entry,index,weights);

            entry = AddWeight(globals,highE,lowZ,
                              upperEnergy*lowerZenith,
                              upperBinEnergy*lowerBinZenith,
                              dUpperE*dLowerZ,
                              entry,index,weights);

            entry = AddWeight(globals,lowE,highZ,
                              lowerEnergy*upperZenith,
                              lowerBinEnergy*upperBinZenith,
                              dLowerE*dUpperZ,
                              entry,index,weights);

            if (lastEntry == entry) break;
        } while (++iz < globals.oscZenith.size());
#ifdef INTERPOLATE_ONLY
#warning Skip smoothing in energy
        break;
#endif
        if (upperBinEnergy<windowThres and lowerBinEnergy<windowThres) break;
        if (upperEnergy<windowThres and lowerEnergy<windowThres) break;
    }

    // Find the total weight, and the weight of the maximum entry.
    double maxWeight = 0.0;
    double totalWeight = 0.0;
    for (int i = 0; i < entry; ++i) {
        maxWeight = std::max(maxWeight, weights[i]);
        totalWeight += weights[i];
    }

    double weightThres = 0.01*totalWeight;
    if (maxWeight > 0.3*totalWeight) {
        weightThres = std::max(weightThres,0.05*maxWeight);
    }
    weightThres = 0.05*maxWeight;

#ifdef SORTED_WEIGHTS
    // Fill work area with weights in decreasing order.
    static std::vector<double> workArea;
    workArea.clear();
    for (int i=0; i<entry; ++i) workArea.emplace_back(weights[i]);
    std::sort(workArea.begin(), workArea.end(),std::greater{}); // reverse sort

    // Find the weight threshold.
    {
        double s = 0.0;
        for (int i = 0; i < workArea.size(); ++i) {
            weightThres = workArea[i];
            s += workArea[i];
            if (i > 4 and s > 0.5*totalWeight) break;
        }
    }
#endif

    // Copy weights equal to or above the threshold.
    int iEntry = 0;
    for (int i = 0; i < entry; ++i) {
        if (weights[i] < weightThres) continue;
        weights[iEntry] = weights[i];
        index[iEntry] = index[i];
        ++iEntry;
    }

    // Fix the normalization.
    double sum = 0;
    for (int i = 0; i < iEntry; ++i) sum += weights[i];
    for (int i = 0; i < iEntry; ++i) weights[i] /= sum;

#ifdef DEBUG_WEIGHT_TABLE
    std::cout << "reweight "
              << energyValue << " " << zenithValue
              << " with " << iEntry << "/" << entry << " entries"
              << ", sum " << sum
              << ", and norm " << norm
              << std::endl;
#endif

    return iEntry;
}

// Provide the binTable entry point required by the GUNDAM tabulated dials.
//
// Find the bin in the table.  The table will always be for a single neutrino
// type and has the order {(Z0,E0) to (Z0,En); (Z1,E0) to (Z1,En); ...}
extern "C"
double binTable(const char* name,
                int varc, double varv[],
                int bins) {
    double energyValue = varv[0];
    double zenithValue = (varc<2) ? -1.0: varv[1];

    TabulatedNuOscillator::TableGlobals& globals
        = TabulatedNuOscillator::globalLookup[name];

    int expectedBins = std::max((std::size_t) 1, globals.oscZenith.size());
    expectedBins = expectedBins * globals.oscEnergies.size();

    if (bins != expectedBins) {
        LIB_CERR << "Table " << name << "has wrong number of bins."
                 << " (Bins: " << bins
                 << " Expected: " << expectedBins << ")"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Find the index for the zenith cosine in the oscZenith table.  The index
    // will be zero if the zenith angle table doesn't exist.  This rounds to
    // the nearest value.
    int zenithIndex = 0;
    if (globals.oscZenith.size() > 1) {
        double value = zenithValue - 1.0/globals.oscZenith.size();
        auto zenithIt = std::lower_bound(globals.oscZenith.begin(),
                                         globals.oscZenith.end(),
                                         value);
        if (zenithIt == globals.oscZenith.end()) {
            zenithIt = globals.oscZenith.end()-1;
        }
        zenithIndex = std::distance(globals.oscZenith.begin(), zenithIt);
    }

    // Find the index for the energy in the oscEnergies table.
    auto energyIt = std::upper_bound(globals.oscEnergies.begin(),
                                     globals.oscEnergies.end(),
                                     energyValue);
    if (energyValue < *(globals.oscEnergies.begin()+1)) {
        energyIt = globals.oscEnergies.begin()+1;
    }
    else if (energyIt == globals.oscEnergies.end()) {
        energyIt = globals.oscEnergies.end()-1;
    }

    int energyIndex = std::distance(globals.oscEnergies.begin(), (energyIt-1));

    double energyBase = zenithIndex*globals.oscEnergies.size() + energyIndex;
    double energyFrac = (energyValue-*(energyIt-1))/(*(energyIt)-*(energyIt-1));

    energyFrac = std::max(0.0,std::min(1.0,energyFrac));

#ifdef DEBUG_EVENT_BINNING
    LIB_COUT << "values: " << energyValue << "," << zenithValue
             << " lo E: " << *(energyIt-1)
             << " hi E: " << *(energyIt)
             << " Z: " << zenithIndex
             << " bin: " << energyBase + energyFrac
             << std::endl;
#endif

    return energyBase + energyFrac;
}

// Provide the updateTable entry point required by the GUNDAM tabulated dials.
extern "C"
int updateTable(const char* name,
                double table[], int bins,
                const double par[], int npar) {

#ifdef DEBUG_UPDATE_TABLE
    LIB_COUT << "Fill table " << name
             << " @ " << (void*) table
             << " bins: " << bins << std::endl;

    LIB_COUT << "    PAR --";
    for (int i = 0; i<npar; ++i) {
        std::cout << " " << i << ": " << par[i];
    }
    std::cout << std::endl;
#endif

    // Sanity check on all of the parameter values;
    for (int i = 0; i<npar; ++i) {
        if (not std::isnan(par[i])) continue;
        LIB_CERR << ": " << name
                 << " WITH NAN PARAMETER VALUE: " << i
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }

    TabulatedNuOscillator::TableGlobals& globals
        = TabulatedNuOscillator::globalLookup[name];
    std::string configName = globals.nuOscillatorConfig;
    TabulatedNuOscillator::NuOscillatorConfig& config
        = TabulatedNuOscillator::configLookup[configName];

    bool oscParamsFilled = false;
    ////////////////////////////////////////////////////////////////////////
    // NuOscillator reverses the convention on the label index order for
    // delta-mass-squared.  NuOscillator kDM23 is (m3^2 - m2^2) which is the
    // PDG value for DM32.
    ////////////////////////////////////////////////////////////////////////
#ifdef UseNuFASTLinear
    if (config.oscillator->ReturnImplementationName()
        .find("Unbinned_NuFASTLinear") != std::string::npos) {
        oscParamsFilled = true;
        // This one only works for LBL neutrino oscillations
        using Calcer = OscProbCalcerNuFASTLinear;
        if (Calcer::kNOscParams != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        config.oscParams[Calcer::kTH12]
            = std::max(par[config.oscParIndex.ss12],1E-12);
        config.oscParams[Calcer::kTH13]
            = std::max(par[config.oscParIndex.ss13],1E-12);
        config.oscParams[Calcer::kTH23]
            = std::max(par[config.oscParIndex.ss23],1E-12);
        // NuOscillator reverses the index order on delta-mass-squared
        config.oscParams[Calcer::kDM12] = par[config.oscParIndex.dm21];
        config.oscParams[Calcer::kDM23] = par[config.oscParIndex.dm32];
        config.oscParams[Calcer::kDCP] = par[config.oscParIndex.dcp];
        config.oscParams[Calcer::kPATHL] = config.oscPath;
        config.oscParams[Calcer::kDENS] = config.oscDensity;
        config.oscParams[Calcer::kELECDENS] = config.oscElectronDensity;
#ifdef DEBUG_NUFAST_PARAMS
        LIB_COUT << "kTH12 " << config.oscParams[Calcer::kTH12] << std::endl;
        LIB_COUT << "kTH23 " << config.oscParams[Calcer::kTH23] << std::endl;
        LIB_COUT << "kTH13 " << config.oscParams[Calcer::kTH13] << std::endl;
        LIB_COUT << "kDM12 " << config.oscParams[Calcer::kDM12] << std::endl;
        LIB_COUT << "kDM23 " << config.oscParams[Calcer::kDM23] << std::endl;
        LIB_COUT << "kPATHL " << config.oscParams[Calcer::kPATHL] << std::endl;
        LIB_COUT << "kDENS " << config.oscParams[Calcer::kDENS] << std::endl;
        LIB_COUT << "kELECDENS " << config.oscParams[Calcer::kELECDENS]
                 << std::endl;
#endif
    }
#endif
#ifdef UseProb3ppLinear
    if (config.oscillator->ReturnImplementationName()
        .find("Unbinned_Prob3ppLinear") != std::string::npos) {
        oscParamsFilled = true;
        // This one only works for LBL neutrino oscillations
        using Calcer = OscProbCalcerProb3ppLinear;
        if (Calcer::kNOscParams != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        config.oscParams[Calcer::kTH12] = par[config.oscParIndex.ss12];
        config.oscParams[Calcer::kTH13] = par[config.oscParIndex.ss13];
        config.oscParams[Calcer::kTH23] = par[config.oscParIndex.ss23];
        config.oscParams[Calcer::kDM12] = par[config.oscParIndex.dm21];
        config.oscParams[Calcer::kDM23] = par[config.oscParIndex.dm32];
        config.oscParams[Calcer::kDCP] = par[config.oscParIndex.dcp];
        config.oscParams[Calcer::kPATHL] = config.oscPath;
        config.oscParams[Calcer::kDENS] = config.oscDensity;
    }
#endif
#ifdef UseOscProb
    if (config.oscillator->ReturnImplementationName()
        .find("Unbinned_OscProb") != std::string::npos) {
        oscParamsFilled = true;
        // This one only works for atmospheric neutrino oscillations
        using Calcer = OscProbCalcerOscProb;
        if (Calcer::kNOscParams+1 == config.oscillator->ReturnNOscParams()) {
            config.oscParams[Calcer::kNOscParams] = config.oscProdHeight;
        }
        else if (Calcer::kNOscParams+3
                 == config.oscillator->ReturnNOscParams()) {
            config.oscParams[Calcer::kNOscParams] = config.oscPath;
            config.oscParams[Calcer::kNOscParams+1] = config.oscDensity;
            config.oscParams[Calcer::kNOscParams+2] = 0.5;
        }
        else if (Calcer::kNOscParams
                 != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << Calcer::kNOscParams
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        config.oscParams[Calcer::kTH12] = par[config.oscParIndex.ss12];
        config.oscParams[Calcer::kTH13] = par[config.oscParIndex.ss13];
        config.oscParams[Calcer::kTH23] = par[config.oscParIndex.ss23];
        config.oscParams[Calcer::kDM12] = par[config.oscParIndex.dm21];
        config.oscParams[Calcer::kDM23] = par[config.oscParIndex.dm32];
        config.oscParams[Calcer::kDCP] = par[config.oscParIndex.dcp];
    }
#endif
#ifdef UseCUDAProb3
    if (config.oscillator->ReturnImplementationName()
        .find("Unbinned_CUDAProb3") != std::string::npos) {
        oscParamsFilled = true;
        // This one only works for atmospheric neutrino oscillations
        using Calcer = OscProbCalcerCUDAProb3;
        if (7 != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << 7
                     << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: "
                     << config.oscillator->ReturnNOscParams()
                     << " Needed: " << 7
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        config.oscParams[0] = par[config.oscParIndex.ss12];
        config.oscParams[2] = par[config.oscParIndex.ss13];
        config.oscParams[1] = par[config.oscParIndex.ss23];
        config.oscParams[3] = par[config.oscParIndex.dm21];
        config.oscParams[4] = par[config.oscParIndex.dm32];
        config.oscParams[5] = par[config.oscParIndex.dcp];
        config.oscParams[6] = config.oscProdHeight;
    }
#endif

    if (not oscParamsFilled) {
        LIB_COUT << "Incorrect oscillation configuration" << std::endl;
        LIB_CERR << "Incorrect oscillation configuration" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    try {
        // See if the table needs to be recalculated.  The oscillator is
        // clever and will only recalculate if the parameters have changed.
        config.oscillator->CalculateProbabilities(config.oscParams);
    } catch (...) {
        LIB_CERR << "Invalid probability or other throw from NuOscillator"
                 << std::endl;
        LIB_CERR << "Filling table " << name
                 << " @ " << (void*) table
                 << " bins: " << bins << std::endl;
        for (int i = 0; i<npar; ++i) {
            LIB_CERR << "     Parameter: " << i
                     << " is " << par[i] << std::endl;
        }
        std::exit(EXIT_FAILURE);
    }

    for (int i = 0; i<bins; ++i) table[i] = 0.0;
    for (TabulatedNuOscillator::TableGlobals::OscWeight &weight
             : globals.weightAddress) {
        const std::size_t i = weight.index;
        const double v = *weight.address;
        const double w = weight.weight;
#ifdef ERROR_CHECKING
        if (i < 0 or bins <= i or w < 0.0 or w > 1.0) {
            LIB_CERR << "Error filling " << name << std::endl;
            LIB_CERR << "    Expecting bin: 0 <= " << i << " < " << bins
                     << std::endl;
            LIB_CERR << "    Expecting weight 0.0 < " << w << " < " << 1.0
                     << std::endl;
            std::exit(EXIT_FAILURE);
        }
        if (not std::isfinite(v) or v < 0.0 or v > 1.0) {
            LIB_CERR << "Error filling " << name << std::endl;
            for (int j = 0; j < npar; ++j) {
                LIB_CERR << "   Parameter " << j
                         << " is " << par[j] << std::endl;
            }
            LIB_CERR << "Table bin is " << i << std::endl;
            LIB_CERR << "Oscillation weight is " << v << std::endl;
            LIB_CERR << "Smoothing weight is " << w << std::endl;
            std::exit(EXIT_FAILURE);
        }
#endif
        table[i] += w*v;
    }

#ifdef DEBUG_UPDATE_TABLE
    std::cout << "Table: " << name << std::endl;
    for (int i=0; i<bins; ++i) {
        std::cout << "  bin" << i
                  << " value: " << table[i]
                  << " energy: " << globals.oscEnergies[i]
                  << std::endl;
    }
#endif

    return 0;
}
