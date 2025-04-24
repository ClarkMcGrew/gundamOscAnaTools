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
    std::vector<FLOAT_T>& energies, double eMin, double eMax) {
    double step = (1.0/eMin - 1.0/eMax)/(energies.size()-1);
    for (std::size_t bin = 0; bin < energies.size(); ++bin) {
        double v = 1.0/eMax + step*bin;
        energies[bin] = 1.0/v;
    }
}

void TabulatedNuOscillator::FillLogarithmicEnergyArray(
    std::vector<FLOAT_T>& energies, double eMin, double eMax) {
    double eMaxLog = std::log(eMax);
    double eMinLog = std::log(eMin);
    double step=(eMaxLog - eMinLog)/(energies.size()-1);
    for (std::size_t bin = 0; bin < energies.size(); ++bin) {
        double v = eMinLog + step*bin;
        energies[bin] = std::exp(v);
    }
}

void TabulatedNuOscillator::FillEnergyArray(
    std::vector<FLOAT_T>& energies, const std::string& type,
    double eMin, double eMax) {
    if (type.find("inv") != std::string::npos) {
        FillInverseEnergyArray(energies,eMin,eMax);
    }
    else {
        FillLogarithmicEnergyArray(energies,eMin,eMax);
    }

    // NuOscillator needs the energies in increasing order.
    std::sort(energies.begin(), energies.end());
    energies[0] = eMin;
    energies[energies.size()-1] = eMax;
}

void TabulatedNuOscillator::FillZenithArray(std::vector<FLOAT_T>& zenith) {
    if (zenith.size() < 1) return;
    if (zenith.size() < 2) {
        zenith[0] = -1.0;
        return;
    }
    double minPath = RoughZenithPath(1.0);
    double maxPath = RoughZenithPath(-1.0);
    double step = (maxPath - minPath)/(zenith.size()-1);
    double maxCosZStep = 0.5/std::sqrt(zenith.size());
    double path = minPath;
    double lastC = 1.0;
    int bin = 0;
    zenith[bin++] = lastC;
    for (double c = zenith[0]; c > -1.0; c -= 1E-8) {
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
    zenith[0] = -1.0;
    zenith[zenith.size()-1] = 1.0;
}

double TabulatedNuOscillator::RoughZenithPath(double cosz) {
    const double Rd{6371}; //Average Earth Radius in km (average)
    const double Rp{Rd + 30.0};
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
// ENERGY_SMOOTH <integer> -- Number of energy bins to smooth over (def: 0)
// ZENITH_BINS <integer>  -- Number of zenith cosine bins (def: 0)
// ZENITH_SMOOTH <integer> -- Number of zenith bins to smooth over (def: 0).
// DENSITY <double>    -- Density in gm/cc
// ELECTRON_DENSITY <double> -- Almost always 0.5
// PATH <double>       -- Path length in kilometers (for LBL)
// PRODUCTION_HEIGHT <double>  -- Neutrino production height in km (for ATM)
//
extern "C"
int initializeTable(const char* name, int argc, const char* argv[],
                    int suggestedBins) {
    LIB_COUT << "Initialize: " << name << std::endl;
    TabulatedNuOscillator::TableGlobals& globals = TabulatedNuOscillator::globalLookup[name];

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
    globals.oscEnergySmooth = 0;
    globals.oscZenithSmooth = 0;

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


    // Get the smoothing for the osc weights in energies
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_SMOOTH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergySmooth;
        break;
    }

    // Get the number of oscillation bins per neutrino
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_BINS") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithBins;
        break;
    }

    // Get the smoothing for the osc weights in energies
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_SMOOTH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithSmooth;
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
                                           globals.oscMaxEnergy);

    TabulatedNuOscillator::FillZenithArray(globals.oscZenith);

#ifdef DUMP_ENERGY_BINNING
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

    int zSmooth = globals.oscZenithSmooth;
    int eSmooth = globals.oscEnergySmooth;
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

    TabulatedNuOscillator::TableGlobals& globals = TabulatedNuOscillator::globalLookup[name];

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

    TabulatedNuOscillator::TableGlobals& globals = TabulatedNuOscillator::globalLookup[name];
    std::string configName = globals.nuOscillatorConfig;
    TabulatedNuOscillator::NuOscillatorConfig& config = TabulatedNuOscillator::configLookup[configName];

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
            LIB_COUT << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
            std::exit(EXIT_FAILURE);
        }
        config.oscParams[Calcer::kTH12] = std::max(par[config.oscParIndex.ss12],1E-12);
        config.oscParams[Calcer::kTH13] = std::max(par[config.oscParIndex.ss13],1E-12);
        config.oscParams[Calcer::kTH23] = std::max(par[config.oscParIndex.ss23],1E-12);
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
        LIB_COUT << "kELECDENS " << config.oscParams[Calcer::kELECDENS] << std::endl;
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
            LIB_COUT << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
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
        else if (Calcer::kNOscParams+3 == config.oscillator->ReturnNOscParams()) {
            config.oscParams[Calcer::kNOscParams] = config.oscPath;
            config.oscParams[Calcer::kNOscParams+1] = config.oscDensity;
            config.oscParams[Calcer::kNOscParams+2] = 0.5;
        }
        else if (Calcer::kNOscParams != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << Calcer::kNOscParams  << std::endl;
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
    //std::cout<<"MyDebug: "<<config.oscillator->ReturnImplementationName()<<std::endl;
    if (config.oscillator->ReturnImplementationName()
        .find("Unbinned_CUDAProb3") != std::string::npos) {
        oscParamsFilled = true;
        // This one only works for atmospheric neutrino oscillations
        using Calcer = OscProbCalcerCUDAProb3;
        if (7 != config.oscillator->ReturnNOscParams()) {
            LIB_COUT << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << 7  << std::endl;
            LIB_CERR << "Wrong number of parameters.  Provided: " << config.oscillator->ReturnNOscParams() << " Needed: " << 7  << std::endl;
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
        // See if the table needs to be recalculated.  The oscillator is clever
        // and will only recalculate if the parameters have changed.
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

#ifdef DUMP_UPDATE_TABLE
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
