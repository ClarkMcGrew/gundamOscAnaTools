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
              std::string binFile,
              std::string binName) {
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

void UpdateAllTables(std::vector<double>& par) {
    for (TableEntry& t : gOscTables) {
        t.updateFunc(t.name.c_str(),
                     t.table.data(), t.table.size(),
                     par.data(), par.size());
    }
}

TableEntry& FindTable(const std::string& config,
                      const std::string& flux,
                      const std::string& interaction) {
    for (TableEntry& t : gOscTables) {
        if (t.config.find(config) == std::string::npos) continue;
        if (t.flux.find(flux) == std::string::npos) continue;
        if (t.interaction.find(interaction) == std::string::npos) continue;
        return t;
    }
    std::cout << "Could not find " << config
              << ":" << flux
              << ":" << interaction
              << std::endl;
    std::exit(1);
}

std::tuple<double,double,double>
CheckTables(const std::string& c1, const std::string& c2,
            const std::string& flux, const std::string& interaction) {
    TableEntry& t1 = FindTable(c1,flux,interaction);
    TableEntry& t2 = FindTable(c2,flux,interaction);
    if (t1.table.size() != t2.table.size()) {
        std::cout << "Wrong table sizes"
                  << std::endl;
        std::exit(1);
    }
    double maxDiff = 0.0;
    double avgDiff = 0.0;
    double rmsDiff = 0.0;
    for (int i=0; i<t1.table.size(); ++i) {
        double delta = t1.table[i] - t2.table[i];
        double sum = 1.0;
#define RELATIVE_DIFFERENCE
#ifdef RELATIVE_DIFFERENCE
        sum =  0.5*(t1.table[i] + t2.table[i]);
#endif
        maxDiff = std::max(maxDiff, std::abs(delta));
        avgDiff += delta/sum;
        rmsDiff += delta*delta/sum/sum;
    }
    avgDiff /= t1.table.size();
    rmsDiff /= t1.table.size();
#ifdef VERBOSE
    std::cout << "Check " << c1
              << " & " << c2
              << " " << flux << " -> " << interaction
              << " max diff: " << maxDiff
              << std::endl;
#endif

    return {maxDiff,avgDiff,std::sqrt(rmsDiff)};
}

std::tuple<double,double,double>
CheckAllTables(const std::string& c1, const std::string c2) {
    double maxTmp;
    double avgTmp;
    double rmsTmp;
    double maxDiff = 0.0;
    double avgDiff = 0.0;
    double avgSum = 0.0;
    double rmsDiff = 0.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"muon","muon");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"electron","electron");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"muon","electron");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"electron", "muon");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"anti-muon","anti-muon");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"anti-electron","anti-electron");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"anti-muon","anti-electron");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    std::tie(maxTmp, avgTmp, rmsTmp) = CheckTables(c1,c2,"anti-electron", "anti-muon");
    maxDiff = std::max(maxDiff,std::abs(maxTmp));
    avgDiff += avgTmp;
    rmsDiff += rmsTmp*rmsTmp;
    avgSum += 1.0;

    return {maxDiff, avgDiff/avgSum, std::sqrt(rmsDiff/avgSum)};
}

std::vector<double> MakeOscPar(const std::vector<double>& base,
                               std::default_random_engine& engine) {
    std::normal_distribution<double> normal{0.0,1.0};
    std::vector<double> par = base;
    do {par[0] = base[0]*(1.0+0.7*normal(engine));}
    while (par[0] < 0.1*base[0] or 0.7 < par[0]);
    do {par[1] = base[1]*(1.0+0.7*normal(engine));}
    while (par[1] < 0.1*base[1] or 0.7 < par[1]);
    do {par[2] = base[2]*(1.0+0.7*normal(engine));}
    while (par[2] < 0.1*base[2] or 0.7 < par[2]);
    do {par[3] = base[3]*(1.0+0.7*normal(engine));}
    while (par[3] < 0.1*base[3]);
    do {par[4] = base[4]*(1.0+0.7*normal(engine));}
    while (par[4] < 0.1*base[4]);
    par[5] = base[4] + normal(engine);
    return par;
}

int main(int argc, char** argv) {

    std::random_device random_seed;
    std::default_random_engine engine{random_seed()};
    std::uniform_real_distribution<double> uniform{0.0,1.0};
    std::normal_distribution<double> normal{0.0,1.0};
    std::lognormal_distribution<double> lognormal{0.0, 0.75};

    auto makeEnergy = [&](){
        const double energyMedian = 5.0;
        return energyMedian*lognormal(engine);
    };

    std::vector<std::string> configs;
#ifdef TestNuFASTLinear
#warning "Test NuFAST"
    {
        configs.push_back("NuFASTLinear");
        std::string enr{"1000"};
        std::string zen{""};
        AddTable("./Configs/GUNDAM_NuFASTLinear","muon","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","muon","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","electron","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","electron","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","anti-muon","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","anti-muon","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","anti-electron","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_NuFASTLinear","anti-electron","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
    }
#else
#warning "Not testing NUFAST"
#endif

#ifdef TestProbGPULinear
#warning "Test ProbGPULinear"
    {
        configs.push_back("ProbGPULinear");
        std::string enr{"1000"};
        std::string zen{""};
        AddTable("./Configs/GUNDAM_ProbGPULinear","muon","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","muon","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","electron","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","electron","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","anti-muon","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","anti-muon","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","anti-electron","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_ProbGPULinear","anti-electron","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
    }
#else
#warning "Not testing ProbGPULinear"
#endif

#ifdef TestProb3ppLinear
#warning "Test Prob3ppLinear"
    {
        configs.push_back("Prob3ppLinear");
        std::string enr{"1000"};
        std::string zen{""};
        AddTable("./Configs/GUNDAM_Prob3ppLinear","muon","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","muon","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","electron","electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","electron","muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","anti-muon","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","anti-muon","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","anti-electron","anti-electron",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
        AddTable("./Configs/GUNDAM_Prob3ppLinear","anti-electron","anti-muon",
                 "./Configs/exampleEnergyBinning.root","energyBinning");
    }
#else
#warning "Not testing Prob3ppLinear"
#endif

    // The PDG values for oscillation parameters
    std::vector<double> pdgPar = {3.07e-1,  // ss12
                                  5.28e-1,  // ss23
                                  2.18e-2,  // ss13
                                  7.53e-5,  // dms21
                                  2.509e-3, // dms32
                                  -1.601};  // dcp

    for (int conf1 = 0; conf1 < configs.size(); ++conf1) {
        for (int conf2 = conf1; conf2 < configs.size(); ++conf2) {
            std::cout << "Check " << configs[conf1]
                      << " and " << configs[conf2]
                      << std::endl;
            UpdateAllTables(pdgPar);
            double maxDiff;
            double avgDiff;
            double rmsDiff;
            std::tie(maxDiff, avgDiff, rmsDiff)
                = CheckAllTables(configs[conf1],configs[conf2]);
            rmsDiff *= rmsDiff;
#define TRIALS
#ifdef TRIALS
            double trials = 1;
            for (int i = 0; i< 100; ++i) {
                std::vector<double> par = MakeOscPar(pdgPar,engine);
                // for (double p : par) std::cout << " " << p;
                // std::cout << std::endl;
                UpdateAllTables(par);
                const auto [m,a,r] = CheckAllTables(configs[conf1],configs[conf2]);
                maxDiff = std::max(maxDiff,m);
                avgDiff += a;
                rmsDiff += r*r;
                trials += 1;
            }
            avgDiff /= trials;
            rmsDiff /= trials;
            std::cout << "Average Difference : " << avgDiff
                      << " for " << trials << " trials"
                      << std::endl;
            std::cout << "RMS Difference : " << std::sqrt(rmsDiff)
                      << std::endl;

#endif
            std::cout << "Maximum Difference : " << maxDiff << std::endl;
        }
    }
    std::exit(EXIT_SUCCESS);
}
