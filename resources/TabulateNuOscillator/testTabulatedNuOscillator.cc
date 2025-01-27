#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <random>
#include <chrono>

#include <dlfcn.h>

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
    initFunc_arguments.push_back("MIN_ENERGY 0.50");
    initFunc_arguments.push_back("MAX_ENERGY 50.0");
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
}

void PrintHist(const std::map<int,int>& hist, double binning) {
    double maxBin = 1.0;
    for (auto [x,y] : hist) maxBin = std::max(maxBin,(double)y);
    for (auto [x,y] : hist) {
        std::cout << std::setw(6)
                  << x*binning
                  << ' ' << std::string(30.0*y/maxBin, '*') << std::endl;
    }
}

void FillHist(std::map<int,int>& hist, double binning, double value) {
    ++hist[std::round(value/binning)];
}

double CreateHist(std::map<int,int>& hist, double maxBin, int bins) {
    double binning = maxBin/bins;
    for (int i=0; i<bins+1; ++i) hist[i] = 0;
    return binning;
}

int main(int argc, char** argv) {

#define TestNUFAST
#ifdef TestNUFAST
    AddTable("./Configs/GUNDAM_NuFASTLinear","muon","muon","1000","");
    AddTable("./Configs/GUNDAM_NuFASTLinear","muon","electron","1000","");
    AddTable("./Configs/GUNDAM_NuFASTLinear","electron","muon","1000","");
    AddTable("./Configs/GUNDAM_NuFASTLinear","anti-muon","anti-muon","1000","");
    AddTable("./Configs/GUNDAM_NuFASTLinear","anti-muon","anti-electron","1000","");
    AddTable("./Configs/GUNDAM_NuFASTLinear","anti-electron","anti-muon","1000","");
#endif

#ifdef TestOscProb
    AddTable("./Configs/GUNDAM_OscProb","muon","muon","100","20");
    AddTable("./Configs/GUNDAM_OscProb","muon","electron","100","20");
    AddTable("./Configs/GUNDAM_OscProb","muon","tau","100","20");
    AddTable("./Configs/GUNDAM_OscProb","electron","electron","100","20");
    AddTable("./Configs/GUNDAM_OscProb","electron","muon","100","20");
    AddTable("./Configs/GUNDAM_OscProb","electron","tau","100","20");
#endif

    std::random_device random_seed;
    std::default_random_engine engine{random_seed()};
    std::uniform_real_distribution<double> uniform{0.0,1.0};
    std::normal_distribution<double> normal{0.0,1.0};
    std::lognormal_distribution<double> lognormal{0.0, 0.75};

    auto makeEnergy = [&](){
        const double energyMedian = 5.0;
        return energyMedian*lognormal(engine);
    };

    // Check binning
    double energyMax = 20.0;
    std::map<int,int> energyHist;
    double energyBin = CreateHist(energyHist, energyMax, 40);
    for (int i=0; i<10000; ++i) {
        double v = makeEnergy();
        if (v > energyMax) continue;
        FillHist(energyHist, energyBin, v);
    }
    PrintHist(energyHist,energyBin);

    // Check the binning function.
    for (TableEntry& t : gOscTables) {
        std::cout << "Test " << t.name
                  << " bins " << t.table.size()
                  << std::endl;
        for (int i=0; i<10; ++i) {
            double energy = makeEnergy()-1.0;
            double cosZenith = 2.0*uniform(engine) - 1.0;
            std::vector<double> argv = {energy, cosZenith};
            double bin = t.binningFunc(t.name.c_str(),
                                       argv.size(),argv.data(),
                                       t.table.size());
            if (bin < 0.0 or bin > t.table.size()-1) {
                std::cout << "Bad bin: Must satisfity 0.0 <= " << bin
                          << " <= " << t.table.size()-1
                          << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }

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
    }

    // Dump the oscillation probabilities for NuFAST
    auto default_precision{std::cout.precision()};
    std::cout << std::fixed
              << std::setprecision(6);
    for (int i = 0; i<9999; ++i) {
        std::cout << "Entry " << i << " ";
        for (TableEntry& t : gOscTables) {
            if (t.config.find("NuFASTLinear") == std::string::npos) continue;
            if (t.table.size() <= i) {
                std::cout << "END";
                i = 9999;
                break;
            }
            std::cout  << std::setw(10) << t.table[i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::defaultfloat
              << std::setprecision(default_precision)
              << std::setw(0);

    // Iterate toward no oscillations
    int iter = 0;
    for (double ss23 = 0.60; ss23 > 1e-107; ss23 *= 0.7) {
        par = pdgPar;
        par[0] = ss23;
        par[1] = ss23;
        par[2] = ss23;
        std::cout << "SS23 " << ss23;
        for (TableEntry& t : gOscTables) {
            if (t.config.find("NuFASTLinear") == std::string::npos) continue;
            t.updateFunc(t.name.c_str(),
                         t.table.data(), t.table.size(),
                         par.data(), par.size());
            std::cout  << std::setw(15)
                       << t.table[800];
        }
        std::cout << std::endl;
    }

    int iterations = 10000;
    std::cout << "Time " << iterations << " NuFASTLinear iterations"
              << " (takes several seconds)" << std::endl;
    // Time the calls
    auto t1 = high_resolution_clock::now();

    for (int i=0; i<iterations; ++i) {
        par[4] = 1.0E-4*normal(engine) + 2.5E-3;
        for (TableEntry& t : gOscTables) {
            if (t.config.find("NuFASTLinear") == std::string::npos) continue;
            t.updateFunc(t.name.c_str(),
                         t.table.data(), t.table.size(),
                         par.data(), par.size());
        }
    }
    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> elapsed = t2 - t1;

    std::cout << "NuFASTLinear Elapsed time: " << elapsed.count() << " ms total"
              << " " << 1000*elapsed.count()/iterations << " us per iteration"
              << std::endl;

    iterations = 100;
    std::cout << "Time " << iterations << " OscProb iterations"
              << " (takes several seconds)" << std::endl;
    // Time the calls
    t1 = high_resolution_clock::now();

    for (int i=0; i<iterations; ++i) {
        par[4] = 1.0E-4*normal(engine) + 2.5E-3;
        for (TableEntry& t : gOscTables) {
            if (t.config.find("OscProb") == std::string::npos) continue;
            t.updateFunc(t.name.c_str(),
                         t.table.data(), t.table.size(),
                         par.data(), par.size());
        }
    }
    t2 = high_resolution_clock::now();

    elapsed = t2 - t1;

    std::cout << "OscProb Elapsed time: " << elapsed.count() << " ms total"
              << " " << elapsed.count()/iterations << " ms per iteration"
              << std::endl;

    std::exit(EXIT_SUCCESS);
}
