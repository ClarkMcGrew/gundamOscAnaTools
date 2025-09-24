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

std::vector<double> eLikeEnergyBins{
    0.0, 0.07, 0.10, 0.15, 0.22, 0.33, 0.50, 0.70, 1.1, 1.6,
    2.3, 3.4, 5.0, 7.6, 11.0, 17.0, 25.0, 37.0, 54.0, 100.0,
};

std::vector<double> muLikeEnergyBins{
    0.0, 0.12, 0.20, 0.33, 0.54, 0.90, 1.5, 2.5,
    4.0, 7.0, 11.0, 19.0, 32.0, 100.0,
};

std::vector<double> eLikeZenithBins{
    -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
};

std::vector<double> muLikeZenithBins{
    -1.0, -0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0,
};

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
        float zenith;
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
    initFunc_arguments.push_back("PARAMETERS SS12,SS23,SS13,DM21,DM32,DCP,SIGN32");
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

    // Get the Krigging weight function
    symbol = dlsym(tableEntry.library, "weightTable");
    if( symbol == nullptr ){
        std::cout << "Weighting function symbol not found: "
                  << "weightTable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    tableEntry.weightingFunc
        = reinterpret_cast<
            int(*)(const char* name, int bins,
                   int varc, double varv[],
                   int entries, int index[], double weights[])>(symbol);

    std::cout << "weightTable address: " << (void*) tableEntry.weightingFunc
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

    {
        TabulatedNuOscillator::TableGlobals& global
            = (*tableEntry.globals)[tableEntry.name];
        TabulatedNuOscillator::NuOscillatorConfig& config
            = (*tableEntry.configs)[global.nuOscillatorConfig];
        std::cout << "Nu Oscillator Config for Krigging "
                  << global.nuOscillatorConfig << std::endl;
        int brake = 0;
        int stepE = 1000;
        int stepZ = 1000;
        std::cout << "Krig energy range " << config.energies.front()
                  << " " << config.energies.back()
                  << std::endl;
        double eMin = std::log10(config.energies.front());
        double eMax = std::log10(config.energies.back());
        std::cout << "Krig zenith range " << config.zenith.front()
                  << " " << config.zenith.back()
                  << std::endl;
        double zMin = config.zenith.front();
        double zMax = config.zenith.back();
        int totalEntries = 0;
        const int tableSize = 100000;
        int index[tableSize];
        double weights[tableSize];
        for (int ie = 0; ie<stepE; ++ie) {
            double e = eMin + ie*(eMax-eMin)/(stepE);
            e = std::exp(e*std::log(10.0));
            for (int iz = 0; iz < stepZ; ++iz) {

                double z = zMin + iz*(zMax-zMin)/(stepZ);
                double varv[]{e,z};
                int entries = tableEntry.weightingFunc(tableEntry.name.c_str(),
                                                       tableEntry.table.size(),
                                                       2, varv,
                                                       tableSize,index,weights);
                totalEntries += entries;
                TableEntry::KrigWeight w;
                w.eMin = eMin;
                w.eMax = eMax;
                w.zMin = zMin;
                w.zMax = zMax;
                w.stepE = stepE;
                w.stepZ = stepZ;
                w.ie = ie;
                w.iz = iz;
                w.energy = e;
                w.zenith = z;
                for (int i = 0; i<entries; ++i) {
                    w.index.emplace_back(index[i]);
                    w.weights.emplace_back(weights[i]);
                }
#ifdef DUMP_KRIG
                if (brake++ < 2000) {
                    std::cout << "KRIG " << w.energy << " " << w.zenith
                              << " " << ie << " " << iz;
                    for (int i = 0; i < w.index.size(); ++i) {
                        std::cout << " (" << w.index[i]
                                  << "," << w.weights[i] << ")";
                    }
                    std::cout << std::endl;
                }
#endif
                tableEntry.krigging.emplace_back(w);
            }
        }
        std::cout << "Krigging weights: " << tableEntry.krigging.size()
                  << " with " << totalEntries << " entries"
                  << " " << 1.0*totalEntries/tableEntry.krigging.size()
                  << " per weight" << std::endl;
    }

    for (auto [table, global] : *tableEntry.globals) {
        std::cout << "TABLE: " << table << " " << global.name << std::endl;
        std::cout << " ENERGY GRID: " << global.oscEnergies.size() << std::endl;
        std::cout << " ZENITH GRID: " << global.oscZenith.size() << std::endl;
        TabulatedNuOscillator::NuOscillatorConfig& config
            = (*tableEntry.configs)[global.nuOscillatorConfig];
        std::cout << " config energies: " << config.energies.size()
                  << std::endl;
    }
}

double RoughZenithPath(double cosz) {
    const double Rd{6371}; //Average Earth Radius in km (average)
    const double Rp{Rd + 30.0}; //Should be thicker than the atmosphere.
    double L = std::sqrt(Rd*Rd*(cosz*cosz-1.0) + Rp*Rp) - Rd*cosz;
    return L;
}

double TableLookup(int enr, int zen,
                   std::vector<double>& table,
                   std::vector<FLOAT_T> energies,
                   std::vector<FLOAT_T> zenith) {
    if ( 0 < zen and zenith.size() <= zen) {
        std::cout << "Zenith angle bin out of bounds: " << zen << std::endl;
        std::exit(EXIT_FAILURE);
    }
    if (enr < 0 or energies.size() <= enr) {
        std::cout << "Energy bin out of bounds: " << enr << std::endl;
        std::exit(EXIT_FAILURE);
    }
    int index = zen*energies.size() + enr;
    if (index < 0 or table.size() <= index) {
        std::cout << "Table index out of bounds: " << index << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return table.at(index);
}

double Krig(std::vector<double>& table, TableEntry::KrigWeight& w) {
    double v = 0.0;
    for (int i = 0; i< w.index.size(); ++i) {
        if (w.weights[i] < 0 or w.weights[i] > 1) {
            std::cout << "Krig weight " << w.weights[i] << std::endl;
        }
        v += w.weights[i]*table.at(w.index[i]);
    }
    return v;
}

void PlotProbabilities(std::string name, std::vector<double> table,
                       std::vector<TableEntry::KrigWeight>& krigging,
                       OscillatorBase* oscillator,
                       std::string flux, int oscInitialFlavor,
                       std::string inter, int oscFinalFlavor,
                       std::vector<FLOAT_T> energies,
                       std::vector<FLOAT_T> zenith,
                       std::vector<FLOAT_T> params) {
    std::cout << "    Plot " << name << std::endl;
    gStyle->SetCanvasDefH(1000);
    gStyle->SetCanvasDefW(1400);

    std::ostringstream title;
    title << "Oscillation probability for " << flux << " to " << inter;

    TGraph2D energyZenithPlot;
    energyZenithPlot.SetNpx(std::min((std::size_t)500,energies.size()));
    energyZenithPlot.SetNpy(std::min((std::size_t)500,zenith.size()));
    TGraph2D energyCosZPlot;
    energyCosZPlot.SetNpx(std::min((std::size_t)500,energies.size()));
    energyCosZPlot.SetNpy(std::min((std::size_t)500,zenith.size()));
    {
        int ezPlot = 0;
        int ecPlot = 0;
        for (int i = 0; i < energies.size(); ++i) {
            for (int j = 0; j < zenith.size(); ++j) {
                double e = energies[i];
                double z = zenith[j];
                double p = TableLookup(i,j,table,energies,zenith);
                energyZenithPlot.SetPoint(ezPlot++, std::log10(e),
                                          RoughZenithPath(z),
                                          p);
                energyCosZPlot.SetPoint(ecPlot++, std::log10(e),
                                        z,
                                        p);
            }
        }
    }
    energyZenithPlot.Draw("colz");
    energyZenithPlot.SetTitle(title.str().c_str());
    energyZenithPlot.GetXaxis()->SetTitle("Energy (log10(GeV))");
    energyZenithPlot.GetYaxis()->SetTitle("Rough Path length (km)");
    gPad->Update();
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
    energyCosZPlot.GetYaxis()->SetTitle("Zenith Cosine");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "OscProbCosZ-" << flux
              << "-" << inter
              << ".png";
        gPad->Print(fName.str().c_str());
    }

    std::cout << "Krig oscillogram" << std::endl;
    TH2D krigEnergyCosZ("krigEnergyCosZ",
                        "Krigged Oscillation Probability",
                        krigging[0].stepE, krigging[0].eMin, krigging[0].eMax,
                        krigging[0].stepZ, krigging[0].zMin, krigging[0].zMax);
    {
        int ezPlot = 0;
        for (TableEntry::KrigWeight& w : krigging) {
            double p = Krig(table,w);
            if (p < 0 or (p-1.0) > 1.0E-7) {
                std::cout << "Probability out of accuracy bounds "
                          << p
                          << " " << p-1.0
                          << std::endl;
            }
            krigEnergyCosZ.SetBinContent(w.ie+1, w.iz+1, p);
        }
    }
    krigEnergyCosZ.Draw("colz");
    krigEnergyCosZ.SetTitle(title.str().c_str());
    krigEnergyCosZ.GetXaxis()->SetTitle("Energy (log10(GeV))");
    krigEnergyCosZ.GetYaxis()->SetTitle("Zenith Cosine");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "OscKrigCosZ-" << flux
              << "-" << inter;
        gPad->Print((fName.str()+".png").c_str());
        TFile output((fName.str()+".root").c_str(),"recreate");
        krigEnergyCosZ.Write();
        output.Close();
    }

    TProfile2D oscProfile("OscProbProfile",title.str().c_str(),
                          20, -1.0, 2.01, 20, -1.0, 1.0);
    for (double e = 0.1; e < 100.0; e *=1.005) {
        for (double z = -1.0; z < 1.0; z += 0.005) {
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
        fName << "OscProfCosZ-" << flux
              << "-" << inter;
        gPad->Print((fName.str()+".png").c_str());
        TFile output((fName.str()+".root").c_str(),"recreate");
        oscProfile.Write();
        output.Close();
    }

    TGraph zenithProbLong;
    {
        int ie = 0;
        std::ostringstream tmp;
        tmp << "Probability at " << int(1000*energies[ie]) << " MeV";
        zenithProbLong.SetTitle(title.str().c_str());
        zenithProbLong.GetYaxis()->SetTitle(tmp.str().c_str());
        zenithProbLong.GetXaxis()->SetTitle("Zenith direction (cosine)");
        for (int i = 0; i < zenith.size(); ++i) {
            double e = energies[ie];
            double z = zenith[i];
            if (z > -0.90) break;
            double l = RoughZenithPath(z);
            double p = TableLookup(0,i,table,energies,zenith);
            zenithProbLong.SetPoint(i, z, p);
        }
    }

    zenithProbLong.Draw("AL*");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "ZenithProbLong-" << flux
              << "-" << inter
              << ".png";
        gPad->Print(fName.str().c_str());
    }

    TGraph zenithProbShort;
    {
        int ie = 0;
        std::ostringstream tmp;
        tmp << "Probability at " << int(1000*energies[ie]) << " MeV";
        zenithProbShort.SetTitle(title.str().c_str());
        zenithProbShort.GetYaxis()->SetTitle(tmp.str().c_str());
        zenithProbShort.GetXaxis()->SetTitle("Zenith direction (cosine)");
        int pnt = 0;
        for (int i = 0; i < zenith.size(); ++i) {
            double e = energies[ie];
            double z = zenith[i];
            double l = RoughZenithPath(z);
            if (1000 < l) continue;
            double p = TableLookup(0,i,table,energies,zenith);
            zenithProbShort.SetPoint(pnt++, z, p);
        }
    }
    zenithProbShort.Draw("AL*");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "ZenithProbShort-" << flux
              << "-" << inter
              << ".png";
        gPad->Print(fName.str().c_str());
    }

    TGraph energyProbLow;
    {
        std::ostringstream tmp;
        tmp << "Probability at cosZ of " << zenith[0];
        energyProbLow.SetTitle(title.str().c_str());
        energyProbLow.GetYaxis()->SetTitle(tmp.str().c_str());
        energyProbLow.GetXaxis()->SetTitle("Energy (GeV)");
    }
    for (int i = 0; i < energies.size(); ++i) {
        double e = energies[i];
        double z = zenith[0];
        double l = RoughZenithPath(z);
        double p = TableLookup(i,0,table,energies,zenith);
        energyProbLow.SetPoint(i, e, p);
        if (e > 0.120) break;
    }
    energyProbLow.Draw("AL*");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "EnergyProbLow-" << flux
              << "-" << inter
              << ".png";
        gPad->Print(fName.str().c_str());
    }

    TGraph energyProbHigh;
    {
        std::ostringstream tmp;
        tmp << "Probability at cosZ of " << zenith[0];
        energyProbHigh.SetTitle(title.str().c_str());
        energyProbHigh.GetYaxis()->SetTitle(tmp.str().c_str());
        energyProbHigh.GetXaxis()->SetTitle("Energy log10(GeV)");
        int pnt = 0;
        for (int i = 0; i < energies.size(); ++i) {
            double e = energies[i];
            if (e < 1.1) continue;
            double z = zenith[0];
            double l = RoughZenithPath(z);
            double p = TableLookup(i,0,table,energies,zenith);
            // energyProbHigh.SetPoint(pnt++, e), p);
            energyProbHigh.SetPoint(pnt++, std::log10(e), p);
            // energyProbHigh.SetPoint(pnt++, 1/e, p);
        }
    }
    energyProbHigh.Draw("AL*");
    gPad->Update();
    {
        std::ostringstream fName;
        fName << "EnergyProbHigh-" << flux
              << "-" << inter
              << ".png";
        gPad->Print(fName.str().c_str());
    }
}

int main(int argc, char** argv) {
    std::string enrSmt{"0.4"};
    std::string enrRes{"0.05"};
    std::string zenSmt{"800"};
    std::string zenRes{"0.0"};
    std::string oscer{"cudaprob3"};

    if (argc > 1) enrRes = argv[1];
    if (argc > 2) zenRes = argv[2];
    if (argc > 3) enrSmt = argv[3];
    if (argc > 4) zenSmt = argv[4];

#ifdef TestOscProb
#warning "Test OscProb"
    if (oscer == "oscprob") {
        std::cout << "Testing OscProb" << std::endl;
        AddTable("./Configs/GUNDAM_OscProb","muon","muon",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_OscProb","muon","electron",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_OscProb","muon","tau",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_OscProb","electron","electron",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_OscProb","electron","muon",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_OscProb","electron","tau",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
    }
#else
#warning "Not testing OscProb"
#endif

#ifdef TestCUDAProb3Fixed
#warning "Test CUDAProb3 with Fixed Production Height"
    if (oscer == "cudaprob3") {
        std::cout << "Testing CUDAProb3 with Fixed Production Height"
                  << std::endl;
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","muon","muon",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
#define ALL_HISTOGRAMS
#ifdef ALL_HISTOGRAMS
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","muon","electron",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","muon","tau",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","electron","electron",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","electron","muon",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_fixed","electron","tau",
                 "./Configs/exampleZenithBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
#endif
    }
#else
#warning "Not testing CUDAProb3 with fixed height"
#endif

#ifdef TestCUDAProb3Height
#warning "Test CUDAProb3 with Production Height"
    if (oscer == "cudaprob3") {
        std::cout << "Testing CUDAProb3 with Production Height"
                  << std::endl;
        AddTable("./Configs/GUNDAM_CUDAProb3_height","muon","muon",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
#define ALL_HISTOGRAMS
#ifdef ALL_HISTOGRAMS
        AddTable("./Configs/GUNDAM_CUDAProb3_height","muon","electron",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_height","muon","tau",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_height","electron","electron",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_height","electron","muon",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
        AddTable("./Configs/GUNDAM_CUDAProb3_height","electron","tau",
                 "./Configs/exampleHeightBinning.root","ProductionHeight_dummy",
                 enrSmt,enrRes,zenSmt,zenRes);
#endif
    }
#else
#warning "Not testing CUDAProb3 with production heights"
#endif

    // The PDG values for oscillation parameters
    std::vector<double> pdgPar = {3.07e-1,  // ss12
                                  5.28e-1,  // ss23
                                  2.18e-2,  // ss13
                                  7.53e-5,  // dms21
                                  2.509e-3, // dms32
                                  -1.601,   // dcp
    };

    std::vector<double> pdgParNO = {3.07e-1,  // ss12
                                    5.28e-1,  // ss23
                                    2.18e-2,  // ss13
                                    7.53e-5,  // dms21
                                    2.509e-3, // dms32
                                    -1.601,   // dcp
                                    1.0       // sign32
    };

    std::vector<double> pdgParIH = {3.07e-1,  // ss12
                                    5.28e-1,  // ss23
                                    2.18e-2,  // ss13
                                    7.53e-5,  // dms21
                                    2.509e-3, // dms32
                                    -1.601,   // dcp
                                    -1.0       // sign32
    };

    std::vector<double> nullOsc = {0.0,  // ss12
                                   0.0,  // ss23
                                   0.0,  // ss13
                                   0.0,  // dms21
                                   0.0,  // dms32
                                   0.0,  // dcp
                                   0.0   // sign32
    };

    std::cout << "TABLES INITIALIZED" << std::endl;

    // Check the update function.
#ifdef DEBUG_NULL_OSCILLATIONS
    std::vector<double> par = nullOsc;
#else
    std::vector<double> par = pdgParIH;
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
        PlotProbabilities(global.name, t.table, t.krigging, config.oscillator,
                          t.flux, global.oscInitialFlavor,
                          t.interaction, global.oscFinalFlavor,
                          config.energies, config.zenith, config.oscParams);
    }

    std::exit(EXIT_SUCCESS);
}
