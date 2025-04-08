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

#include <Oscillator/OscillatorFactory.h>
#ifdef UseNuFASTLinear
#warning Including NuFASTLinear in the build
#include <OscProbCalcer/OscProbCalcer_NuFASTLinear.h>
#endif
#ifdef UseProb3ppLinear
#warning Including Prob3PlusPlus in the build
#include <OscProbCalcer/OscProbCalcer_Prob3ppLinear.h>
#endif
#ifdef UseOscProb
#warning Including OscProb in the build
#include <OscProbCalcer/OscProbCalcer_OscProb.h>
#endif
//#undef UseCUDAProb3     // (as of 25/01) BUG in OscProbCalcer_CUDAProb3.h
#ifdef UseCUDAProb3
#warning Including CUDAProb3 in the build
#include <OscProbCalcer/OscProbCalcer_CUDAProb3.h>
#endif

#define LIB_NAME "TabulatedNuOscillator"
#ifndef LOUD_AND_PROUD
#define LOUD_AND_PROUD true
#endif

#define LIB_COUT if (LOUD_AND_PROUD) std::cout << LIB_NAME << " -- "
#define LIB_CERR std::cerr << "ERROR: " << LIB_NAME << " -- "

namespace {
    struct OscillationParameters {
        int ss12;
        int ss13;
        int ss23;
        int dm21;
        int dm32;
        int dcp;
    };

    struct NuOscillatorConfig {
        std::string name;      // The configuration file to use.
        OscillationParameters oscParIndex;
        double oscPath;
        double oscDensity;
        double oscElectronDensity;
        double oscProdHeight;
        // NuOscillator interface values:
        //    -- FLOAT_T is defined in OscillatorConstants.h (no namespace).
        std::unique_ptr<OscillatorBase> oscillator;
        std::vector<FLOAT_T> energies; // The energies for each bin
        std::vector<FLOAT_T> zenith;   // The cosines for each bin
        std::vector<FLOAT_T> oscParams;
    };
    using ConfigLookup = std::map<std::string, NuOscillatorConfig>;
    ConfigLookup configLookup;

    struct TableGlobals {
        std::string name;           // The table name
        std::vector<std::string> arguments; // initialization arguments
        std::string nuOscillatorConfig;     // The configuration file to use.
        double oscMinEnergy;       // Minimum neutrino energy (GeV)
        double oscMaxEnergy;       // Maximum neutrino energy (GeV)
        double oscMaxLoverE;       // Osc table upper limit
        std::string oscEnergyStep; // The type of step (inverse or logarithmic)
        int oscEnergyBins;         // Number of energy bins per neutrino type.
        int oscZenithBins;         // Number of angle bins per neutrino type.
        OscillationParameters oscParIndex;
        // NuOscillator interface values:
        //    -- FLOAT_T is defined in OscillatorConstants.h (no namespace).
        //
        // The flavors are defined by NuOscillators with values of
        // NuOscillator::kElectron==1, NuOscillator::kMuon==2, and
        // NuOscillator::kTau==3.  Anti-neutrinos are specified using a
        // negative value.  The oscInitialFlavor and oscFinalFlavor must have
        // the same sign to be valid.
        std::unique_ptr<OscillatorBase> oscillator;
        int oscInitialFlavor;       // Flavor of the parent (neg. for anti)
        int oscFinalFlavor;         // Flaver of the interacting
        FLOAT_T oscDensity;         // The density along the path (gm/cc)
        FLOAT_T oscElectronDensity; // The electron density (usually 0.5).
        FLOAT_T oscPath;            // The path length for the table (km).
        FLOAT_T oscProdHeight;      // The production height for the table (km).
        std::vector<FLOAT_T> oscEnergies; // Energies for bins
        std::vector<FLOAT_T> oscZenith; // Zenith cosines for bins (optional)
        std::vector<const FLOAT_T*> weightAddress; // NuOscillator to Tabulate map
    };
    std::map<std::string,TableGlobals> globalLookup;

    void FillInverseEnergyArray(std::vector<FLOAT_T>& energies,
                                double eMin, double eMax) {
        double step = (1.0/eMin - 1.0/eMax)/(energies.size()-1);
        for (std::size_t bin = 0; bin < energies.size(); ++bin) {
            double v = 1.0/eMax + step*bin;
            energies[bin] = 1.0/v;
        }
    }

    void FillLogarithmicEnergyArray(std::vector<FLOAT_T>& energies,
                                    double eMin, double eMax) {
        double step=(std::log(eMax)-std::log(eMin))/(energies.size()-1);
        for (std::size_t bin = 0; bin < energies.size(); ++bin) {
            double v = std::log(eMin) + step*bin;
            energies[bin] = std::exp(v);
        }
    }

    void FillEnergyArray(std::vector<FLOAT_T>& energies,
                         const std::string& type,
                         double eMin, double eMax) {
        if (type.find("inv") != std::string::npos) {
            FillInverseEnergyArray(energies,eMin,eMax);
        }
        else {
            FillLogarithmicEnergyArray(energies,eMin,eMax);
        }

        // NuOscillator needs the energies in increasing order.
        std::sort(energies.begin(), energies.end());
    }

    void FillZenithArray(std::vector<FLOAT_T>& zenith) {
        if (zenith.size() < 1) return;
        if (zenith.size() < 2) {
            zenith[0] = -1.0;
            return;
        }
        double step = 2.0/(zenith.size()-1);
        for (std::size_t bin = 0; bin < zenith.size(); ++bin) {
            zenith[bin] = -1.0 + step*bin;
        }
    }

    bool AlmostEqual(double a, double b) {
        const double diff = std::abs(a-b);
        const double avg = std::abs(a) + std::abs(b);
        const double delta = 1E-10;
        if (diff > delta*avg+delta) return false;
        return true;
    }

    void ConfigureNuOscillator(const TableGlobals& globals) {
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
        newConfig.oscillator.reset(factory->CreateOscillator(newConfig.name));

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
    }
};

// This requires string arguments
//
// CONFIG <file-name>
// FLUX_FLAVOR [anti-]{electron,muon,tau}
// INTERACTION_FLAVOR [anti-]{electron,muon,tau}
// PARAMETERS <list of SS12, SS23, SS13, DM21, DM32, DCP>
// ENERGY_BINS <integer> -- Number of energy bins for each neutrino type
// MIN_ENERGY <double> -- Minimum neutrino energy in GeV
// MAX_ENERGY <double> -- Maximum neutrino energy in GeV
// ENERGY_STEP {inverse,logarithmic} -- The binning to use for energy
// ZENITH_BINS <integer>  -- Number of zenith cosine bins (def: 0)
// ENERGY_STEP inverse or logarithmic -- The type of energy step.
// DENSITY <double>    -- Density in gm/cc
// ELECTRON_DENSITY <double> -- Almost always 0.5
// PATH <double>       -- Path length in kilometers (for LBL)
// PRODUCTION_HEIGHT <double>  -- Neutrino production height in km (for ATM)
//
extern "C"
int initializeTable(const char* name, int argc, const char* argv[],
                    int bins) {
    LIB_COUT << "Initialize" << name << std::endl;
    TableGlobals& globals = globalLookup[name];

    // Set default values.
    globals.oscEnergyBins = 1000;
    globals.oscZenithBins = 0;
    globals.oscPath = 1300.0; // km
    globals.oscProdHeight = 17.0; // km
    globals.oscMinEnergy = 0.050; // GeV
    globals.oscDensity = 2.6; // gm/cc
    globals.oscElectronDensity = 0.5;
    globals.oscEnergyStep = "inverse";

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

    // Get the type of step for the energy array.
    for (std::string arg: globals.arguments) {
        if (arg.find("ENERGY_STEP") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscEnergyStep;
        break;
    }

    // Get the number of oscillation bins per neutrino
    for (std::string arg: globals.arguments) {
        if (arg.find("ZENITH_BINS") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscZenithBins;
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

    FillEnergyArray(globals.oscEnergies, globals.oscEnergyStep,
                    globals.oscMinEnergy,globals.oscMaxEnergy);

    FillZenithArray(globals.oscZenith);

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
    NuOscillatorConfig& config = configLookup[globals.nuOscillatorConfig];

    int zenithIndex = 0;
    double zenith = -999.0;
    while (true) {
        if (zenithIndex < globals.oscZenith.size()) {
            zenith = globals.oscZenith[zenithIndex];
        }
        for (double energy : globals.oscEnergies) {
            const FLOAT_T* address
                = config.oscillator->ReturnWeightPointer(
                    globals.oscInitialFlavor,globals.oscFinalFlavor,
                    energy, zenith);
            globals.weightAddress.emplace_back(address);
        }
        ++zenithIndex;
        if (zenithIndex < globals.oscZenith.size()) continue;
        break;
    }

    LIB_COUT << "Table size: " << globals.weightAddress.size() << std::endl;

    return globals.weightAddress.size();
}

// Find the bin in the table.  The table will always be for a single neutrino
// type and has the order {(Z0,E0) to (Z0,En); (Z1,E0) to (Z1,En); ...}
extern "C"
double binTable(const char* name,
                int varc, double varv[],
                int bins) {
    double energyValue = varv[0];
    double zenithValue = (varc<2) ? -1.0: varv[1];

    TableGlobals& globals = globalLookup[name];

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

    TableGlobals& globals = globalLookup[name];
    std::string configName = globals.nuOscillatorConfig;
    NuOscillatorConfig& config = configLookup[configName];

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

    // Copy the table.
    if (bins != globals.weightAddress.size()) {
        LIB_CERR << "Mismatched table size"
                 << std::endl;
        std::exit(EXIT_FAILURE);
    }

    for (int i=0; i<bins; ++i) {
        double v = *globals.weightAddress[i];
        if (not std::isfinite(v) or v < 0.0 or v > 1.0) {
            LIB_CERR << "Error filling " << name << std::endl;
            for (int j = 0; j < npar; ++j) {
                LIB_CERR << "   Parameter " << j << " is " << par[j] << std::endl;
            }
            LIB_CERR << "Weight is " << v << std::endl;
            std::exit(EXIT_FAILURE);
        }
        table[i] = v;
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
