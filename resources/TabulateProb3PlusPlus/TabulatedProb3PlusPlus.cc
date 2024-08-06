// Fill a table of oscillation weights for the GUNDAM Tabulated dial type
// using the Prob3PlusPlus library found at
//
// https://github.com/rogerwendell/Prob3plusplus.git
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <vector>
#include <map>
#include <algorithm>
#include <tuple>
#include <cmath>

#include <BargerPropagator.h>

#define LIB_NAME "TabulatedProb3PlusPlus -- "
#undef LOUD_AND_PROUD

struct OscillationParameters {
    int index_ss12;
    int index_ss13;
    int index_ss23;
    int index_dm21;
    int index_dm32;
    int index_dcp;
};

// A tuple holding the flux information.  <0> is the energy, and <1> is the
// flux (arbitrary units since this is only used in ratios.
using FluxElement = std::tuple<double,double>;

// The fluxes for each flavor of neutrinos.  The vector of flux elements is
// sorted by energy (element <0>).
using FluxTable = std::map<int,std::vector<FluxElement>>;

// A tuple holding the flux weighting intformation: <0> is the energy; <1> is
// the weight; and <2> is the parameter index.
using WeightElement = std::tuple<double,double,int>;

// The weights for each flavor of neutrino.  The vector is sorted by the
// neutrino energy (element <0>).
using WeightTable = std::map<int, std::vector<WeightElement>>;

struct TableGlobals {
    std::string name;           // The table name
    std::vector<std::string> arguments; // initialization arguments
    double oscPath;             // The path length for the table (km).
    double oscDensity;          // The density along the pat.
    double oscMinEnergy;        // Minimum neutrino energy (GeV)
    int oscNeutrinos;           // Number of neutrino types.
    int oscBins;                // Number of bins per neutrino type.
    double oscMaxLoverE;        // Osc table upper limit
    OscillationParameters oscPar; // The MNS parameters
    FluxTable fluxTable;          // Flux for the neutrino vs energy
    WeightTable weightTable;      // Flux weight for the neutrino vs energy
    std::map<int,int> tableIndex;  // Index in the osc table for the neutrino
    std::unique_ptr<BargerPropagator> propagator;
};
std::map<std::string,TableGlobals> globalLookup;

struct FileLine {
    // The field name which is the first string-like entry on the line.
    std::string field;
    // The numerical values come aftger the name.  Depends on doubles
    // being exact for all reasonable integer values.
    std::vector<double> values;
};

// Yes, this is returning a copy by value.
std::vector<FileLine> parseFile(std::string fileName) {
    std::ifstream input(fileName.c_str());

    if (not input.is_open()) {
        std::cout << LIB_NAME << "File not found: " << fileName << std::endl;
        return std::vector<FileLine>();
    }

    std::vector<FileLine> lines;
    std::string line;
    while (std::getline(input,line)) {
        // truncate at the comment character
        std::size_t end = line.find("#");
        line = line.substr(0,end);
        if (line.empty()) continue;
        // replace delimiters with space
        std::replace(line.begin(), line.end(),',',' ');  // comman separated
        std::replace(line.begin(), line.end(),'\t',' '); // tab separated
        std::istringstream parser(line);
        FileLine entry;
        parser >> entry.field;
        while (true) {
            double value;
            parser >> value;
            if (not parser.good()) break;
            entry.values.emplace_back(value);
        }
        lines.emplace_back(entry);
    }

    return lines;
}

// Translate between the MC standard neutrino type and the index in the barger
// code.  This is mostly for documentation.
int mcCodeToBarger(int mcCode) {
    if (mcCode == 12) return 1;
    if (mcCode == 14) return 2;
    if (mcCode == 16) return 3;
    if (mcCode == -12) return -1;
    if (mcCode == -14) return -2;
    if (mcCode == -16) return -3;
    throw std::runtime_error("Unknown MC Code for neutrino");
}

// Translate between the index in the barger code and the MC standard neutrino
// types.  This is mostly for documentation.
int BargerToMcCode(int barger) {
    if (barger == 1) return 12;
    if (barger == 2) return 14;
    if (barger == 3) return 16;
    if (barger == -1) return -12;
    if (barger == -2) return -14;
    if (barger == -3) return -16;
    throw std::runtime_error("Unknown Barger index for neutrino");
}

// Do the interpolation in a table to find the value for the energy.
template <typename Table>
double interpolateTable(Table table, double energy) {
    if (table.size() < 1) return 0.0;
    typename Table::iterator elem
        = std::lower_bound(table.begin(), table.end(), energy,
                           [](const typename Table::iterator::value_type& lhs,
                              const double& rhs){
                               return std::get<0>(lhs) < rhs;
                           });
    if (std::next(elem) == table.end()) return std::get<1>(*elem);
    if (energy < std::get<0>(*elem)) return std::get<1>(*elem);
    const double n = energy - std::get<0>(*elem);
    const double d = std::get<0>(*std::next(elem)) - std::get<0>(*elem);
    const double f = n/d;
    double v = f*std::get<1>(*elem) + (1.0-f)*std::get<1>(*std::next(elem));
    if (not std::isfinite(v)) std::runtime_error("Invalid interpolation");
    return v;
}

// Read a flux file and add it to the flux table.
bool readFluxTable(std::string fileName, FluxTable& fluxTable) {
    std::vector<FileLine> lines = parseFile(fileName);
    if (lines.size() < 1) return false;
    // Get the flux table ("FLUX nutype energy flux")
    for (auto &line: lines) {
        if (line.field != "FLUX") continue;
        if (line.values.size() != 3) {
            std::cout << LIB_NAME
                      << " Invalid flux line: " << line.field;
            for(double v: line.values) std::cout << " " << v;
            std::cout << std::endl;
            std::exit(EXIT_FAILURE);
        }
        int type = line.values[0];
        double energy = line.values[1];
        double flux = line.values[2];
        fluxTable[type].emplace_back(energy,flux);
    }
    // Sort the fluxes by neutrino energy.
    for (auto& flux : fluxTable) {
        std::sort(flux.second.begin(), flux.second.end());
    }
    return true;
}

// Read a parameter definition file.  There can be only one parameter
// definition file per TabulatedIterator.
bool readParameterDefinitions(std::string fileName,
                              OscillationParameters& oscPar,
                              WeightTable& weightTable) {
    // Get the parameters
    //
    // SS12, SS13, SS23, DM21, DM32, DCP, WEIGHT
    //
    // If the line is a WEIGHT, it takes three numbers
    // WEIGHT nutype energy default
    oscPar.index_ss12 = -1;
    oscPar.index_ss13 = -1;
    oscPar.index_ss23 = -1;
    oscPar.index_dm21 = -1;
    oscPar.index_dm32 = -1;
    oscPar.index_dcp = -1;
    // Read the parameters!
    int index = 0;
    std::vector<FileLine> lines = parseFile(fileName);
    for (auto& line : lines) {
        if (line.field == "SS12") oscPar.index_ss12 = index++;
        else if (line.field == "SS13") oscPar.index_ss13 = index++;
        else if (line.field == "SS23") oscPar.index_ss23 = index++;
        else if (line.field == "DM21") oscPar.index_dm21 = index++;
        else if (line.field == "DM32") oscPar.index_dm32 = index++;
        else if (line.field == "DCP") oscPar.index_dcp = index++;
        else if (line.field == "WEIGHT" and line.values.size() >= 2) {
            int index_wght = index++;
            int type = line.values[0];
            double energy = line.values[1];
            double weight = 1.0;
            if (line.values.size() > 2) weight = line.values[2];
            weightTable[type].emplace_back(energy,weight,index_wght);
        }
        else {
            std::cout << LIB_NAME
                      << " Invalid parameter line: " << line.field;
            for(double v: line.values) std::cout << " " << v;
            std::cout << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    // Sort the weights by neutrino energy.
    for (auto& weight : weightTable) {
        std::sort(weight.second.begin(), weight.second.end());
    }

    return true;
}

// This requires string arguments
//
// BINS <integer>
// PATH <double>       -- Path length in kilometers
// DENSITY <double>
// MIN_ENERGY <double> -- Minimum neutrino energy in GeV
// FLUX <file-name>
// PARAMETERS <file-name>
//
extern "C"
int initializeTable(const char* name, int argc, const char* argv[],
                    int bins) {
    std::cout << LIB_NAME << "initialize" << std::endl;
    TableGlobals& globals = globalLookup[name];

    // Make sure that the weight and flux tables exist (but are empty)
    globals.weightTable[12]; globals.fluxTable[12];
    globals.weightTable[14]; globals.fluxTable[14];
    globals.weightTable[16]; globals.fluxTable[16];
    globals.weightTable[-12]; globals.fluxTable[-12];
    globals.weightTable[-14]; globals.fluxTable[-14];
    globals.weightTable[-16]; globals.fluxTable[-16];

    for (int i = 0; i < argc; ++i) {
        std::cout << LIB_NAME << "Argument: " << argv[i] << std::endl;
        globals.arguments.emplace_back(argv[i]);
    }

    // Get the number of oscillation bins per neutrino
    for (std::string arg: globals.arguments) {
        if (arg.find("BINS") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscBins;
        break;
    }
    // Get the minimum energy for neutrinos
    for (std::string arg: globals.arguments) {
        if (arg.find("MIN_ENERGY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscMinEnergy;
        break;
    }
    // Get the path length for neutrinos
    for (std::string arg: globals.arguments) {
        if (arg.find("PATH") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscPath;
        break;
    }
    // The material density along the path length
    for (std::string arg: globals.arguments) {
        if (arg.find("DENSITY") != 0) continue;
        std::istringstream tmp(arg);
        tmp >> arg >> globals.oscDensity;
        break;
    }
    // Get the parameter file and read it.
    for (std::string arg: globals.arguments) {
        if (arg.find("PARAMETERS") != 0) continue;
        std::istringstream tmp(arg);
        std::string fileName;
        tmp >> arg >> fileName;
        readParameterDefinitions(fileName, globals.oscPar, globals.weightTable);
    }
    // Get the flux files and read them.
    for (std::string arg: globals.arguments) {
        if (arg.find("FLUX") != 0) continue;
        std::istringstream tmp(arg);
        std::string fileName;
        tmp >> arg >> fileName;
        readFluxTable(fileName,globals.fluxTable);
    }

    // Record where each neutrino type starts in the table.
    int nextIndex = 0;
    for (auto& flux: globals.fluxTable) {
        int neutrinoType = flux.first;
        globals.tableIndex[neutrinoType] = nextIndex;
        nextIndex += globals.oscBins;
    }

    // Set the maximum L/E for the table (the minimum is always zero).
    globals.oscMaxLoverE = globals.oscPath/globals.oscMinEnergy;

    // Count the number of neutrino types base on the number of flux files.
    globals.oscNeutrinos = globals.fluxTable.size();

    // Create the oscillation propagator.
    globals.propagator = std::make_unique<BargerPropagator>();
    globals.propagator->UseMassEigenstates(false);

    return globals.oscNeutrinos*globals.oscBins;
}

extern "C"
double binTable(const char* name,
                int varc, double varv[],
                int bins) {
    TableGlobals& globals = globalLookup[name];
    if (varc < 3) {
        throw std::runtime_error("must provide type, energy and pathlength");
    }

    int neutrinoType = varv[0];
    double neutrinoEnergy = varv[1];
    double neutrinoPath = varv[2];
    double LE = neutrinoPath/neutrinoEnergy;
    if (LE < 0.0) return -1;

    LE = LE/globals.oscMaxLoverE;
    if (LE >= globals.oscBins) LE = globals.oscBins - 1E-8;

    // Neutrino outside of range of oscillations
    std::map<int,int>::iterator elem = globals.tableIndex.find(neutrinoType);
    if (elem == globals.tableIndex.end()) return -1;

    double value = elem->second;
    value += LE;

    return value;
}

extern "C"
int updateTable(const char* name,
                double table[], int bins,
                const double par[], int npar) {
    TableGlobals& globals = globalLookup[name];

    // Sanity check on all of the parameter values;
    for (int i = 0; i<npar; ++i) {
        if (not std::isnan(par[i])) continue;
        std::cerr << LIB_NAME << ": " << name
                  << " NAN PARAMETER VALUE: " << i
                  << std::endl;
        throw std::runtime_error("NAN parameter value");
    }

    if (bins < globals.oscBins*globals.oscNeutrinos) {
        std::runtime_error("Invalid input table (not large enough)");
    }

    // Get the oscillation parameters.
#define CHECK(n) if ((n)<0||(n)>=npar) throw std::runtime_error("missing index")
    CHECK(globals.oscPar.index_ss12);
    CHECK(globals.oscPar.index_ss13);
    CHECK(globals.oscPar.index_ss23);
    CHECK(globals.oscPar.index_dm21);
    CHECK(globals.oscPar.index_dm32);
    CHECK(globals.oscPar.index_dcp);
#undef CHECK

    double ss12 = par[globals.oscPar.index_ss12];
    double ss13 = par[globals.oscPar.index_ss13];
    double ss23 = par[globals.oscPar.index_ss23];
    double dm21 = par[globals.oscPar.index_dm21];
    double dm32 = par[globals.oscPar.index_dm32];
    double dcp = par[globals.oscPar.index_dcp];

    // Fill the weights from the parameters
    for (auto& weight : globals.weightTable) {
        for (auto& entry: weight.second) {
            int index = std::get<2>(entry);
            if (index < 0) continue;
            if (index >= npar) std::runtime_error("Invalid weight");
            std::get<1>(entry) = par[index];
        }
    }

    // The flux for neutrino type 1, 2, and 3 (i.e. 12, 14 and 16 -- meaning
    // e, mu, and tau).
    double currentFlux[3];

    // NEUTRINOS: Loop over each energy, calculate the oscillation weights,
    // calculate the fluxes, and fill the table.
    for (int i = 0; i < globals.oscBins; ++i) { // From 1 since 0 is always 1.0
        if (i < 1) {
            table[globals.tableIndex[12]] = 1.0;
            table[globals.tableIndex[14]] = 1.0;
            table[globals.tableIndex[16]] = 1.0;
            continue;
        }
        double currentLoE = i*globals.oscMaxLoverE/globals.oscBins;
        double currentEnu = globals.oscPath/currentLoE;
        currentFlux[0] = interpolateTable(globals.weightTable[12], currentEnu);
        currentFlux[0] *= interpolateTable(globals.fluxTable[12], currentEnu);
        currentFlux[1] = interpolateTable(globals.weightTable[14], currentEnu);
        currentFlux[1] *= interpolateTable(globals.fluxTable[14], currentEnu);
        currentFlux[2] = interpolateTable(globals.weightTable[16], currentEnu);
        currentFlux[2] *= interpolateTable(globals.fluxTable[16], currentEnu);
        globals.propagator->SetMNS(ss12, ss13, ss23, dm21, dm32, dcp,
                                   currentEnu, true);
        globals.propagator->propagateLinear(+1, // neutrinos
                                            globals.oscPath,
                                            globals.oscDensity);
        for (int i = 1; i<=3; ++i) {
            double weight = 0;
            for (int j = 1; j<=3; ++j) {
                weight += globals.propagator->GetProb(i,j)*currentFlux[j-1];
            }
            weight /= currentFlux[i-1];
            table[globals.tableIndex[BargerToMcCode(i)]] = weight;
        }
    }

    // ANTINEUTRINOS: Loop over each energy, calculate the oscillation
    // weights, calculate the fluxes, and fill the table.
    for (int i = 0; i < globals.oscBins; ++i) { // From 1 since 0 is always 1.0
        if (i < 1) {
            table[globals.tableIndex[-12]] = 1.0;
            table[globals.tableIndex[-14]] = 1.0;
            table[globals.tableIndex[-16]] = 1.0;
            continue;
        }
        double currentLoE = i*globals.oscMaxLoverE/globals.oscBins;
        double currentEnu = globals.oscPath/currentLoE;
        currentFlux[0] = interpolateTable(globals.weightTable[-12], currentEnu);
        currentFlux[0] *= interpolateTable(globals.fluxTable[-12], currentEnu);
        currentFlux[1] = interpolateTable(globals.weightTable[-14], currentEnu);
        currentFlux[1] *= interpolateTable(globals.fluxTable[-14], currentEnu);
        currentFlux[2] = interpolateTable(globals.weightTable[-16], currentEnu);
        currentFlux[2] *= interpolateTable(globals.fluxTable[-16], currentEnu);
        globals.propagator->SetMNS(ss12, ss13, ss23, dm21, dm32, dcp,
                                   currentEnu, true);
        globals.propagator->propagateLinear(-1, // antineutrinos
                                            globals.oscPath,
                                            globals.oscDensity);
        for (int i = 1; i<=3; ++i) {
            double weight = 0;
            for (int j = 1; j<=3; ++j) {
                weight += globals.propagator->GetProb(-i,-j)*currentFlux[j-1];
            }
            weight /= currentFlux[i-1];
            table[globals.tableIndex[BargerToMcCode(-i)]] = weight;
        }
    }

    return 0;
}
