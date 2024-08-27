#include <iostream>
#include <dlfcn.h>
#include <vector>

int (*initFunc)(const char* name,
                int argc, const char* argv[],
                int bins);
int (*updateFunc)(const char* name,
                  double table[], int bins,
                  const double par[], int npar);
double (*binningFunc)(const char* name,
                      int varc, double varv[],
                      int bins);

int main(int argc, char** argv) {

    void* library = dlopen("./libTabulatedProb3PlusPlus.so", RTLD_LAZY );
    if( library == nullptr ){
        std::cout << "Cannot load library: " << dlerror() << std::endl;
        std::exit(EXIT_FAILURE);
    }

    void* symbol = dlsym(library, "initializeTable");
    if( symbol == nullptr ){
        std::cout << "Initialization function symbol not found: "
                  << "initializeTable"
                  << std::endl;
        std::exit(EXIT_FAILURE); // Exit, not throw!
    }


    initFunc
        = reinterpret_cast<
            int(*)(const char* name,int argc, const char* argv[], int bins)
        >(symbol);

    std::cout << "initializationTable address: " << (void*) initFunc
              << std::endl;
    std::vector<std::string> initFunc_arguments;
    initFunc_arguments.push_back("NEUTRINO_TYPES 6");
    initFunc_arguments.push_back("BINS 100");
    initFunc_arguments.push_back("MIN_ENERGY 0.05");
    initFunc_arguments.push_back("MAX_PATH 295.0");
    initFunc_arguments.push_back("PARAMETERS ${CONFIG_DIR}/200TabulatedAsimov-parameters.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/nue.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/nuebar.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/numu.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/numubar.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/nutau.txt");
    initFunc_arguments.push_back("FLUX ${DATA_DIR}/nutaubar.txt");
    std::vector<const char*> initFunc_argv;
    for (std::string& arg : initFunc_arguments)
        initFunc_argv.push_back(arg.c_str());

    (*initFunc)("name",initFunc_argv.size(),initFunc_argv.data(),-1);

    // Get the update function
    symbol = dlsym(library, "updateTable");
    if( symbol == nullptr ){
        std::cout << "Update function symbol not found: "
                 << "updateTable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    updateFunc
        = reinterpret_cast<
            int(*)(const char* name, double table[],
                   int bins, const double par[], int npar)>(symbol);

    std::cout << "updateTable address: " << (void*) updateFunc
              << std::endl;

    // Get the binning function
    symbol = dlsym(library, "binTable");
    if( symbol == nullptr ){
        std::cout << "Binning function symbol not found: "
                  << "binTable"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }
    binningFunc
        = reinterpret_cast<
            double(*)(const char* name,
                   int varc, double varv[], int bins)>(symbol);

    std::cout << "binTable address: " << (void*) binningFunc
              << std::endl;

    std::exit(EXIT_SUCCESS);
}
