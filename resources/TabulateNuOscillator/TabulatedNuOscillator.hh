#ifndef TabulateNuOscillator_hh_seen
#define TabulateNuOscillator_hh_seen
// TabulatedNuOscillator is not intended to be included directly into user
// code, and should be accessed vial the defined dlopen/dlsym interface to the
// shared library.
//
// This includes declarations so that TabulatedNuOscillator can be directly
// included in compiled code. This is mainly for debugging since the usual
// interface is via dlopen/dlsym with the shared library (and defined
// interface).
//
// Including TabulatedNuOscillator allows access to "private" internal data,
// but should generally not be done.

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
#ifdef UseCUDAProb3
#warning Including CUDAProb3 in the build
#include <OscProbCalcer/OscProbCalcer_CUDAProb3.h>
#endif

namespace TabulatedNuOscillator {
    struct OscillationParameters {
        int ss12;
        int ss13;
        int ss23;
        int dm21;
        int dm32;
        int dcp;
    };

    /// The payload for a map between the NuOscillator config file used and
    /// values that are used with that config file.  Notice that this will
    /// often contain copies of the values in the global config (which is
    /// indexed by the "Tabular" dial name), but they mean different things.
    /// Several tabular dials can share the same NuOscillatorConfig, but not
    /// all dials need to.
    struct NuOscillatorConfig {
        std::string name;      // The configuration file to use.
        OscillationParameters oscParIndex;
        double oscPath;
        double oscDensity;
        double oscElectronDensity;
        double oscProdHeight;
        // NuOscillator interface values:
        //    -- FLOAT_T is defined in OscillatorConstants.h (no namespace).
#ifdef TABULATED_NUOSCILLATOR_DECONSTRUCTABLE
        // OK when OscillatorBase can be safely deconstructed (it is iffy)
        std::unique_ptr<OscillatorBase> oscillator;
#else
        // Work around OscillatorBase deconstructor bugs.
        OscillatorBase* oscillator;
#endif
        std::vector<FLOAT_T> energies; // The energies for each bin
        std::vector<FLOAT_T> zenith;   // The cosines for each bin
        std::vector<FLOAT_T> oscParams;
    };
    // A map between the NuOscillator config file and the config variables.
    using ConfigLookup = std::map<std::string, NuOscillatorConfig>;
    extern "C" ConfigLookup configLookup;

    // The values associated with a particular tabulated dial.  Notice that
    // this can contain values "shared" with the NuOscillatorConfig, but they
    // mean different things.  Multiple tabulated dials can share the same
    // NuOscillatorConfig, but they do not have to.  Notably, the values here,
    // and in the associated NuOscillatorConfig must match, and there are
    // checks to make sure they do.
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
        int oscEnergySmooth;       // Bins to smooth over for energy
        int oscZenithSmooth;       // Bins to smooth over for zenith
        OscillationParameters oscParIndex;
        // NuOscillator interface values:
        //    -- FLOAT_T is defined in OscillatorConstants.h (no namespace).
        //
        // The flavors are defined by NuOscillators with values of
        // NuOscillator::kElectron==1, NuOscillator::kMuon==2, and
        // NuOscillator::kTau==3.  Anti-neutrinos are specified using a
        // negative value.  The oscInitialFlavor and oscFinalFlavor must have
        // the same sign to be valid.
        int oscInitialFlavor;       // Flavor of the parent (neg. for anti)
        int oscFinalFlavor;         // Flaver of the interacting
        FLOAT_T oscDensity;         // The density along the path (gm/cc)
        FLOAT_T oscElectronDensity; // The electron density (usually 0.5).
        FLOAT_T oscPath;            // The path length for the table (km).
        FLOAT_T oscProdHeight;      // The production height for the table (km).
        std::vector<FLOAT_T> oscEnergies; // Energies for bins
        std::vector<FLOAT_T> oscZenith; // Zenith cosines for bins (optional)
        struct OscWeight {
            std::size_t index;
            const FLOAT_T* address;
            float weight;
        };
        std::vector<OscWeight> weightAddress; // NuOscillator to Tabulate map
    };
    using GlobalLookup = std::map<std::string, TableGlobals>;
    extern "C" GlobalLookup globalLookup;

    void FillInverseEnergyArray(std::vector<FLOAT_T>& energies,
                                double eMin, double eMax);

    void FillLogarithmicEnergyArray(std::vector<FLOAT_T>& energies,
                                    double eMin, double eMax);

    void FillEnergyArray(std::vector<FLOAT_T>& energies,
                         const std::string& type,
                         double eMin, double eMax);

    void FillZenithArray(std::vector<FLOAT_T>& zenith);

    double RoughZenithPath(double cosz);

    bool AlmostEqual(double a, double b);

    void ConfigureNuOscillator(const TableGlobals& globals);
};
#endif
