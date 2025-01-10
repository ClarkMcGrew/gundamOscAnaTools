# Example: Using Prob3PlusPlus to Tabulate Oscillations

This is a simple example of how to use a oscillation weight calculation with GUNDAM.  This example is based on Prob3PlusPlus https://github.com/rogerwendell/Prob3plusplus.git, but the basic principles should apply.

Compiling this code on linux use:

```bash
mkdir build
cd build
cmake ..
make
```

That will grab the Prob3plusplus repository from Roger's github, build it, and then build the libTabulatedProb3PlusPlus.so file that can be included in GUNDAM.  The code is not currently supported on MacOS, but probably could be made to work.

# Controlling TabulatedProb3PlusPlus with the GUNDAM yaml config file

This library is used to implemented a Tabulated dial in the a gundam yaml config file.  This is done under the each parameter set definition.

```yaml
fitterEngineConfig:
  propagatorConfig:
    parameterSetListConfig
    - name: aParameterSetName
      dialSetDefinitions:
        - dialType: Tabulated
          tableConfig:
            name: "A table name"  # must be unique
            libraryPath: "${LIB_PATH}/libTabulatedProb3PlusPlus.so"
            updateFunction: "updateTable"
            binningFunction: "bintable"
            binningVariables:
              - branch-with-nu-type   # Nu MC Code
              - branch-with-nu-energy # Nu energy in GeV
              - branch-with-path      # Pathlength in km
            initFunction: "initializeTable"
            initArguments:
              - "BINS <N>"       # bins per neutrino
              - "MIN_ENERGY <E>" # minimum neutrino energy
              - "PATH <L>"       # nominal path in km
              - "DENSITY <D>"    # the density in ???
              - "PARAMETERS <file>" # file defining fit params
              - "FLUX <file>"    # File with flux table
```

The neutrino MC codes are the usual PDG codes for the neutrinos


| Neutrino Type         | MC Code |
|:----------------------|--------:|
| electron neutrino     |      12 |
| electron antineutrino |     -12 |
| muon neutrino         |      14 |
| muon antineutrino     |     -14 |
| tau neutrino          |      16 |
| tau antineutrino      |     -16 |

# The parameter and flux definition files

The flux and parameter definitions are contained in files that contain single line definitions (for either the flux or the fit parameters).  The basic file format is

```
<NAME> <val> <val> <val>
```

where `<NAME>` is a name defined in the next two subsections, and `<val>`
represents a possible value associated with the name.  A single file can
contain the parameter definitions, and all the flux definitions.  The
parameters need to be defined in one file, but there can be multiple files
defining the fluxes. The definition files use '#' to define comment lines.

## The parameter definitions

The parameter definitions are contained in a text file that specifies how the GUNDAM config file has defined the parameters.  This is the order of the parameters in the input array of doubles.  The oscillation parameters are defined using the strings

* SS12    -- Sin-Squared theta 12.
* SS13    -- Sin-Squared theta 13.
* SS23    -- Sin-Squared theta 23.
* DM21    -- delta M^2 21.
* DM32    -- delta M^2 32.
* DCP     -- delta CP.
* WEIGHT nutype energy(GeV)  -- The flux weight for nutype at energy (repeated).

## The flux tables

The flux table defines the neutrino flux at a particular energy.  It
consists of lines that contain "FLUX" definitions where each definition
contains a neutrino type (mc code), the neutrino energy in GeV, and the
flux at that energy.  Since the flux will always occur in a ratio with
another flux, the units for the flux do not matter.  Comments can be
included by prefixing them with '#'.  For example

```txt
# An electron neutrino flux
FLUX 12 0.1 1.0
FLUX 12 0.2 3.0
FLUX 12 0.4 2.0
FLUX 12 0.5 1.0
FLUX 12 1.5 0.1
```

A flux file can contain multiple neutrino types.
