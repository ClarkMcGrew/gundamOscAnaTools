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

    int bins;
    std::vector<double> table;

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
    TabulatedNuOscillator::GlobalLookup *globals;
    TabulatedNuOscillator::ConfigLookup *configs;
};
std::vector<TableEntry> gOscTables;

void AddTable(std::string config,
              std::string flux,
              std::string interaction,
              std::string energyBins, std::string zenithBins) {
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
    initFunc_arguments.push_back("ENERGY_BINS "+energyBins);
    initFunc_arguments.push_back("ENERGY_STEP inverse");
    initFunc_arguments.push_back("MIN_ENERGY 0.10");
    initFunc_arguments.push_back("MAX_ENERGY 100.0");
    if (not zenithBins.empty()) {
        initFunc_arguments.push_back("ZENITH_BINS "+zenithBins);
    }
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

    std::cout << "Table size: " << tableEntry.table.size() << std::endl;

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

    for (auto [table, global] : *tableEntry.globals) {
        std::cout << "table: " << table << " " << global.name << std::endl;
        std::cout << " energies: " << global.oscEnergies.size() << std::endl;
        std::cout << " zeniths: " << global.oscZenith.size() << std::endl;
        TabulatedNuOscillator::NuOscillatorConfig& config
            = (*tableEntry.configs)[global.nuOscillatorConfig];
        std::cout << " config energies: " << config.energies.size() << std::endl;
    }
}

double RoughZenithPath(double cosz) {
    const double Rd{6371}; //Average Earth Radius in km (average)
    const double Rp{Rd + 30.0};
    double L = std::sqrt(Rd*Rd*(cosz*cosz-1.0) + Rp*Rp) - Rd*cosz;
    return L;
}


void PlotProbabilities(std::string name, OscillatorBase* oscillator,
                       std::string flux, int oscInitialFlavor,
                       std::string inter, int oscFinalFlavor,
                       std::vector<FLOAT_T> energies,
                       std::vector<FLOAT_T> zenith,
                       std::vector<FLOAT_T> params) {
    std::cout << "    Plot " << name << std::endl;
    try {
        std::ostringstream title;
        title << "Oscillation probability for " << flux << " to " << inter;

        oscillator->CalculateProbabilities(params);
        TGraph2D energyZenithPlot;
        energyZenithPlot.SetNpx(std::min((std::size_t)500,energies.size()));
        energyZenithPlot.SetNpy(std::min((std::size_t)500,zenith.size()));
        TGraph2D energyCosZPlot;
        energyCosZPlot.SetNpx(std::min((std::size_t)500,energies.size()));
        energyCosZPlot.SetNpy(std::min((std::size_t)500,zenith.size()));
        int i = 0;
        int j = 0;
        for (double e : energies) {
            for (double z: zenith) {
                double p = oscillator->ReturnOscillationProbability(
                    oscInitialFlavor, oscFinalFlavor, e, z);
                energyZenithPlot.SetPoint(i++, std::log10(e),
                                          RoughZenithPath(z),
                                          p);
                energyCosZPlot.SetPoint(j++, std::log10(e),
                                        z,
                                        p);
            }
        }
        energyZenithPlot.Draw("colz");
        energyZenithPlot.SetTitle(title.str().c_str());
        energyZenithPlot.GetXaxis()->SetTitle("Energy (log10(GeV))");
        energyZenithPlot.GetYaxis()->SetTitle("Rough Path length (km)");
        gPad->Update();
        // energyZenithPlot.GetYaxis()->SetTitle("Zenith Angle (cos(theta))");
        {
            std::ostringstream fName;
            fName << "OscProbLength-" << flux
                  << "-" << inter
                  << ".png";
            gPad->Print(fName.str().c_str());
        }

        energyCosZPlot.Draw("colz");
        energyCosZPlot.SetTitle(title.str().c_str());
        energyCosZPlot.GetXaxis()->SetTitle("Energy (log10(GeV))");
        energyCosZPlot.GetYaxis()->SetTitle("Rough Path length (km)");
        gPad->Update();
        // energyCosZPlot.GetYaxis()->SetTitle("Zenith Angle (cos(theta))");
        {
            std::ostringstream fName;
            fName << "OscProbCosZ-" << flux
                  << "-" << inter
                  << ".png";
            gPad->Print(fName.str().c_str());
        }

        TProfile2D oscProfile("OscProbProfile",title.str().c_str(),
                              50, -1.0, 2.01, 20, -1.0, 1.0);
        for (double e = 0.1; e < 100.0; e *=1.01) {
            for (double z = -1.0; z < 1.0; z += 0.01) {
                double p = energyCosZPlot.Interpolate(std::log10(e),z);
                oscProfile.Fill(std::log10(e),z,p);
            }
        }

        oscProfile.Draw("colz");
        oscProfile.SetTitle(title.str().c_str());
        oscProfile.GetXaxis()->SetTitle("Energy (log10(GeV))");
        oscProfile.GetYaxis()->SetTitle("Zenith Angle (cosine)");
        gPad->Update();
        {
            std::ostringstream fName;
            fName << "OscProfile-" << flux
                  << "-" << inter
                  << ".png";
            gPad->Print(fName.str().c_str());
        }

        TGraph zenithPlot;
        zenithPlot.SetTitle(title.str().c_str());
        zenithPlot.GetYaxis()->SetTitle("Probability at 100 MeV");
        zenithPlot.GetXaxis()->SetTitle("Path Length (km)");
        i = 0;
        for (double z: zenith) {
            double p = oscillator->ReturnOscillationProbability(
                oscInitialFlavor, oscFinalFlavor, 0.1, z);
            double l = RoughZenithPath(z);
            if (l < 12000) continue;
            zenithPlot.SetPoint(i++, l, p);
        }
        zenithPlot.Draw("AC*");
        gPad->Update();
        {
            std::ostringstream fName;
            fName << "ZenithProb-" << flux
                  << "-" << inter
                  << ".png";
            gPad->Print(fName.str().c_str());
        }

        TGraph energyPlot;
        energyPlot.SetTitle(title.str().c_str());
        energyPlot.GetYaxis()->SetTitle("Probabilty at cosZ of -1");
        energyPlot.GetXaxis()->SetTitle("Energy (GeV)");
        i = 0;
        for (double e: energies) {
            if (e > 0.120) continue;
            double p = oscillator->ReturnOscillationProbability(
                oscInitialFlavor, oscFinalFlavor, e, -1.0);
            energyPlot.SetPoint(i++, e, p);
        }
        energyPlot.Draw("AC*");
        gPad->Update();
        {
            std::ostringstream fName;
            fName << "EnergyProb-" << flux
                  << "-" << inter
                  << ".png";
            gPad->Print(fName.str().c_str());
        }
    }
    catch (...) {
        std::cout << "Problem with oscillation calculation" << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv) {
    std::string enr{"100"};
    std::string zen{"100"};
    if (argc > 1) enr = argv[1];
    if (argc > 2) zen = argv[2];

    std::cout << "Generating a " << enr << " x " << zen << " grid" <<std::endl;

#ifdef NOTTestOscProb
#warning "Test OscProb"
    AddTable("./Configs/GUNDAM_OscProb","muon","muon",enr,zen);
    AddTable("./Configs/GUNDAM_OscProb","muon","electron",enr,zen);
    AddTable("./Configs/GUNDAM_OscProb","muon","tau",enr,zen);
    AddTable("./Configs/GUNDAM_OscProb","electron","electron",enr,zen);
    AddTable("./Configs/GUNDAM_OscProb","electron","muon",enr,zen);
    AddTable("./Configs/GUNDAM_OscProb","electron","tau",enr,zen);
#else
#warning "Not testing OscProb"
#endif

#ifdef TestCUDAProb3
#warning "Test CUDAProb3"
    AddTable("./Configs/GUNDAM_CUDAProb3","muon","muon",enr,zen);
    AddTable("./Configs/GUNDAM_CUDAProb3","muon","electron",enr,zen);
    AddTable("./Configs/GUNDAM_CUDAProb3","muon","tau",enr,zen);
    AddTable("./Configs/GUNDAM_CUDAProb3","electron","electron",enr,zen);
    AddTable("./Configs/GUNDAM_CUDAProb3","electron","muon",enr,zen);
    AddTable("./Configs/GUNDAM_CUDAProb3","electron","tau",enr,zen);
#else
#warning "Not testing CUDAProb3"
#endif

    std::random_device random_seed;
    std::default_random_engine engine{random_seed()};
    std::uniform_real_distribution<double> uniform{0.0,1.0};
    std::normal_distribution<double> normal{0.0,1.0};
    std::lognormal_distribution<double> lognormal{0.0, 0.75};

    // The PDG values for oscillation parameters
    std::vector<double> pdgPar = {3.07e-1,  // ss12
                                  5.28e-1,  // ss23
                                  2.18e-2,  // ss13
                                  7.53e-5,  // dms21
                                  2.509e-3, // dms32
                                  -1.601};  // dcp

    // Check the update function.
    std::vector<double> par = pdgPar;
    for (TableEntry& t : gOscTables) {
        std::cout << "Test " << t.name
                  << " update " << t.table.size()
                  << std::endl;
        t.updateFunc(t.name.c_str(),
                     t.table.data(), t.table.size(),
                     par.data(), par.size());
        TabulatedNuOscillator::TableGlobals& global = (*t.globals)[t.name];
        std::cout << "  global name: " << global.name << std::endl;
        TabulatedNuOscillator::NuOscillatorConfig& config = (*t.configs)[global.nuOscillatorConfig];
        PlotProbabilities(global.name, config.oscillator,
                          t.flux, global.oscInitialFlavor,
                          t.interaction, global.oscFinalFlavor,
                          config.energies, config.zenith, config.oscParams);
    }

    std::exit(EXIT_SUCCESS);
}
