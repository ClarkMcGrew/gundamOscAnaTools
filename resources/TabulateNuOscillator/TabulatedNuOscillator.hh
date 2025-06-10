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

class TH1;

namespace TabulatedNuOscillator {
    /// The index of the oscillation parameter in the parameter array from
    /// GUNDAM.  The order is controlled by the "PARAMETERS <list>" string
    /// argument in the initializeTable call.  The values are copied into the
    /// correct NuOscillator parameter locations.
    struct OscillationParameters {
        int ss12;
        int ss13;
        int ss23;
        int dm21;
        int dm32;
        int dcp;
    };

    /// The payload for a map between the NuOscillator config file used and
    /// values that are used with that config file.  This describes the
    /// configuration of the NuOscillator library.  Notice that this will
    /// often contain copies of the values in the global config (which is
    /// indexed by the "Tabular" dial name), but they mean different things.
    /// Several dials can share the same NuOscillatorConfig, but not all dials
    /// need to.
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
    // this will contain values "shared" with the NuOscillatorConfig, but they
    // mean different things.  Multiple dials can share the same
    // NuOscillatorConfig, but they do not have to.  Notably, the values here,
    // and in the associated NuOscillatorConfig must match, and there are
    // checks to make sure they do.
    struct TableGlobals {
        std::string name;           // The table name
        std::vector<std::string> arguments; // initialization arguments
        std::string nuOscillatorConfig;     // The configuration file to use.
        std::string oscBinningFile;         // Name of the binning file
        std::string oscBinningHistName;     // Name of the binning histogram
        std::string oscEnergyType; // The binning (edge, average, log, inverse)
        std::string oscZenithType; // The binning (edge, average)

        // Provide a window in 1/E and path length that is used for smoothing.
        // The actual smoothing is provided by the resolution fields.
        double oscEnergySmooth;    // Smoothing 1/E (1/GeV)
        double oscZenithSmooth;    // Smoothing L (km)

        // Provide the resolution that will be appled to the energy
        // (fractional), and the angular resolution.
        double oscEnergyResol;     // Fractional energy resolution
        double oscZenithResol;     // Angular resolution (radians or cosine)
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

    /// Calculate the difference between e2 and e1 based on the energy binning
    /// used.  This is (1/e2 - 1/e1) and logarithmic binning is not used
    /// (deprecated and considered a mistake).
    void energyBinDelta(double e2, double e1);

    /// Fill a vector with the energies that NuOscillator will used to
    /// calculate the oscillation weights.  This fills (the preferred) 1/E
    /// spacing.  The energy resolution is provided, and is used to limit the
    /// overall step size.
    void FillInverseEnergyArray(std::vector<FLOAT_T>& energies,
                                double eMin, double eMax, double eRes);

    /// Fill a vector with energies.  This uses the (deprecated) log energy
    /// step.  It is provided for testing, but shouldn't be used (much).
    void FillLogarithmicEnergyArray(std::vector<FLOAT_T>& energies,
                                    double eMin, double eMax);

    /// Fill a vector with energies for NuOscillator.  This uses either
    /// FillInverseEnergyArray, or FillLogarithmicEnergyArray.
    void FillEnergyArray(std::vector<FLOAT_T>& energies,
                         const std::string& type,
                         double eMin, double eMax, double eRes);

    /// Fill a vector with energies for NuOscillator.
    void FillEnergyArray(std::vector<FLOAT_T>& energies,
                         const std::string& type,
                         TH1* energyBins);

    /// Calculate the approximate "delta" along the zenith angle axis.  The
    /// table spacing is approximately by path length while the bins are
    /// labeled in cos(zenithAngle).  This is the approximate difference in
    /// path length. This returns the absolute value of the change.
    double zenithBinDelta(double c2, double c1);

    /// Fill a vector with energies for NuOscillator.
    void FillZenithArray(std::vector<FLOAT_T>& zenith,
                         const std::string& type,
                         TH1* zenithBins);

    /// Fill a vector with the zenith angle binning.
    void FillZenithArray(std::vector<FLOAT_T>& zenith);

    /// Calculate an approximateion of the path length for a zenith cosine.
    double RoughZenithPath(double cosz);

    /// Check if two doubles are almost equal.  They both need to be within a
    /// few 1E-10 of the average.
    bool AlmostEqual(double a, double b);

    /// Make sure the globals and nuoscillator configs agree.
    void ConfigureNuOscillator(const TableGlobals& globals);
};
#endif
