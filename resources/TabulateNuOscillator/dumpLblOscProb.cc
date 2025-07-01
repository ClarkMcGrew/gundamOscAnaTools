#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>

#include <dlfcn.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile2D.h>
#include <TPad.h>
#include <TStyle.h>

#include "TabulatedNuOscillator.hh"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

struct TableEntry {
    std::string name;
    std::string config;
    std::string flux;
    std::string interaction;

    std::vector<double> table;

    struct KrigWeight {
        float eMin;
        float eMax;
        float zMin;
        float zMax;
        float energy;
        int stepE;
        int stepZ;
        int ie;
        int iz;
        std::vector<int> index;
        std::vector<float> weights;
    };
    std::vector<KrigWeight> krigging;

    void *library;
    int (*initFunc)(const char* name,
                    int argc, const char* argv[],
                    int bins);
    int (*updateFunc)(const char* name,
                      double table[], int bins,
                      const double par[], int npar);
    double (*binningFunc)(const char* name,
                          int varc, double varv[],
                          int bins);
    int (*weightingFunc)(const char* name, int bins,
                         int varc, double varv[],
                         int entries, int index[], double weights[]);
    TabulatedNuOscillator::GlobalLookup *globals;
    TabulatedNuOscillator::ConfigLookup *configs;
};
std::vector<TableEntry> gOscTables;

void AddTable(std::string config,
              std::string flux,
              std::string interaction,
              std::string binFile,
              std::string binName,
              std::string energySmooth,
              std::string energyResol,
              std::string zenithSmooth,
              std::string zenithResol) {
    gOscTables.emplace_back();
    std::cout << "Add table " << gOscTables.size() << std::endl;

    TableEntry& tableEntry = gOscTables.back();
    tableEntry.name = config + "/" + flux + ":" + interaction ;
    tableEntry.config = config;
    tableEntry.flux = flux;
    tableEntry.interaction = interaction;

    tableEntry.library = dlopen("libTabulatedNuOscillator.so", RTLD_LAZY );
    if( tableEntry.library == nullptr ){
        std::cout << "Cannot load library: " << dlerror() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    void* symbol = dlsym(tableEntry.library, "initializeTable");
    if( symbol == nullptr ){
        std::cout << "Initialization function symbol not found: "
                  << "initializeTable"
                  << std::endl;
        std::exit(EXIT_FAILURE); // Exit, not throw!
    }

    std::cout << "LIBRARY:  " << (void*) tableEntry.library << std::endl;

    tableEntry.initFunc
        = reinterpret_cast<
            int(*)(const char* name,int argc, const char* argv[], int bins)
        >(symbol);

    std::cout << "initializationTable address: " << (void*) tableEntry.initFunc
              << std::endl;
    std::vector<std::string> initFunc_arguments;
    initFunc_arguments.push_back("FLUX_FLAVOR "+tableEntry.flux);
    initFunc_arguments.push_back("INTERACTION_FLAVOR "+tableEntry.interaction);
    initFunc_arguments.push_back("PARAMETERS SS12,SS23,SS13,DM21,DM32,DCP");
    initFunc_arguments.push_back("BINNING_FILE "+binFile);
    initFunc_arguments.push_back("BINNING_HIST "+binName);
    initFunc_arguments.push_back("ZENITH_SMOOTH "+zenithSmooth);
    initFunc_arguments.push_back("ZENITH_RESOLUTION "+zenithResol);
    initFunc_arguments.push_back("ENERGY_SMOOTH "+energySmooth);
    initFunc_arguments.push_back("ENERGY_RESOLUTION "+energyResol);
    initFunc_arguments.push_back("DENSITY 2.6");
    initFunc_arguments.push_back("PATH 1300.0");
    initFunc_arguments.push_back("CONFIG "+tableEntry.config);

    std::vector<const char*> argv;
    for (std::string& arg : initFunc_arguments)
        argv.push_back(arg.c_str());

    int bins
        = (*tableEntry.initFunc)(tableEntry.name.c_str(),
                                 argv.size(),
                                 argv.data(),
                                 -1);
    tableEntry.table.resize(bins);

    std::cout << "Table size: " << tableEntry.table.size()
              << " Returned " << bins
              << std::endl;

    // Get the update function
    symbol = dlsym(tableEntry.library, "updateTable");
    if( symbol == nullptr ){
        std::cout << "Update function symbol not found: "
                  << "updateTable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    tableEntry.updateFunc
        = reinterpret_cast<
            int(*)(const char* name, double table[],
                   int bins, const double par[], int npar)>(symbol);

    std::cout << "updateTable address: " << (void*) tableEntry.updateFunc
              << std::endl;

    // Get the binning function
    symbol = dlsym(tableEntry.library, "binTable");
    if( symbol == nullptr ){
        std::cout << "Binning function symbol not found: "
                  << "binTable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    tableEntry.binningFunc
        = reinterpret_cast<
            double(*)(const char* name,
                   int varc, double varv[], int bins)>(symbol);

    std::cout << "binTable address: " << (void*) tableEntry.binningFunc
              << std::endl;

    // Get some of the TabulateNuOscillator internal symbols for debugging.
    symbol = dlsym(tableEntry.library, "globalLookup");
    if( symbol == nullptr ){
        std::cout << "Global lookup table symbol not found "
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    tableEntry.globals = reinterpret_cast<
        TabulatedNuOscillator::GlobalLookup*>(symbol);

    symbol = dlsym(tableEntry.library, "configLookup");
    if( symbol == nullptr ){
        std::cout << "Config lookup table symbol not found "
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    tableEntry.configs = reinterpret_cast<
        TabulatedNuOscillator::ConfigLookup*>(symbol);

}

double TableLookup(int enr,
                   std::vector<double>& table,
                   std::vector<FLOAT_T> energies) {
    if (enr < 0 or energies.size() <= enr) {
        std::cout << "Energy bin out of bounds: " << enr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    int index = enr;
    if (index < 0 or table.size() <= index) {
        std::cout << "Table index out of bounds: " << index << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return table.at(index);
}

void PlotProbabilities(std::string name, std::vector<double> table,
                       OscillatorBase* oscillator,
                       std::string flux, int oscInitialFlavor,
                       std::string inter, int oscFinalFlavor,
                       std::vector<FLOAT_T> energies,
                       std::vector<FLOAT_T> params) {
    std::cout << "    Plot " << name << std::endl;
    gStyle->SetCanvasDefH(1000);
    gStyle->SetCanvasDefW(1400);

    std::ostringstream title;
    title << "Oscillation probability for " << flux << " to " << inter;

    int iLow = 0;
    int steps = 100;
    while (iLow < energies.size()) {
        TGraph energyPlot;

        int ePlot = 0;
        for (int i = iLow;
             i < std::min((std::size_t) iLow+steps,energies.size());
             ++i) {
            double e = energies[i];
            double p = TableLookup(i,table,energies);
            energyPlot.SetPoint(ePlot++, (e),p);
        }

        energyPlot.Draw("AL*");
        gPad->Update();
        {
            std::ostringstream fName;
            fName << "energyPlot-" << flux
                  << "-" << inter
                  << "_" << iLow
                  << ".png";
            gPad->Print(fName.str().c_str());
        }
        iLow += steps;
    }

}

int main(int argc, char** argv) {
    std::string enrSmt{"0.4"};
    std::string enrRes{"0.05"};
    std::string zenSmt{"800"};
    std::string zenRes{"0.0"};
    std::string oscer{"nufast"};

    if (argc > 1) enrRes = argv[1];
    if (argc > 2) zenRes = argv[2];
    if (argc > 3) enrSmt = argv[3];
    if (argc > 4) zenSmt = argv[4];

#define TestNuFASTLinear
#ifdef  TestNuFASTLinear
#warning "Test NuFASTLinear"
    if (oscer == "nufast") {
        std::cout << "Testing NuFASTLinear" << std::endl;
        AddTable("./Configs/GUNDAM_NuFASTLinear","muon","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_NuFASTLinear","muon","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_NuFASTLinear","muon","tau",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_NuFASTLinear","electron","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_NuFASTLinear","electron","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_NuFASTLinear","electron","tau",
                 "./Configs/exampleEnergyBinning.root","energyBinning",
                 enrSmt,enrRes,zenSmt,zenRes);
    }
#else
#warning "Not testing NuFASTLinear"
#endif

    // The PDG values for oscillation parameters
    std::vector<double> pdgPar = {3.07e-1,  // ss12
                                  5.28e-1,  // ss23
                                  2.18e-2,  // ss13
                                  7.53e-5,  // dms21
                                  2.509e-3, // dms32
                                  -1.601};  // dcp

    std::vector<double> nullOsc = {0.0,  // ss12
                                   0.0,  // ss23
                                   0.0,  // ss13
                                   0.0,  // dms21
                                   0.0, // dms32
                                   0.0};  // dcp

    std::cout << "TABLES INITIALIZED" << std::endl;

    // Check the update function.
#ifdef DEBUG_NULL_OSCILLATIONS
    std::vector<double> par = nullOsc;
#else
    std::vector<double> par = pdgPar;
#endif
    for (TableEntry& t : gOscTables) {
        std::cout << "Test " << t.name
                  << " update " << t.table.size()
                  << std::endl;
        t.updateFunc(t.name.c_str(),
                     t.table.data(), t.table.size(),
                     par.data(), par.size());
        TabulatedNuOscillator::TableGlobals& global = (*t.globals)[t.name];
        std::cout << "  global name: " << global.name << std::endl;
        TabulatedNuOscillator::NuOscillatorConfig& config
            = (*t.configs)[global.nuOscillatorConfig];
        PlotProbabilities(global.name, t.table, config.oscillator,
                          t.flux, global.oscInitialFlavor,
                          t.interaction, global.oscFinalFlavor,
                          config.energies, config.oscParams);
    }

    std::exit(EXIT_SUCCESS);
}
