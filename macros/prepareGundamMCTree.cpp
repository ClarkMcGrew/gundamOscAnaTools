#include <iostream>
#include <string>
#include <TClonesArray.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TGraph.h>
#include <vector>
#include <map>



#include "TSystem.h"

using namespace TMath;

const double DELMSQ_31 = 2.515e-3; // In eV^2
const double LOSC = 1300.; // In km

const double THETA_23 = 0.859;
const double THETA_13 = 0.150;

TString FHCnonswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/FHC/nonswap/FD_FHC_nonswap.root";
TString FHCnueswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/FHC/nueswap/FD_FHC_nueswap.root";
TString FHCtauswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/FHC/tauswap/FD_FHC_tauswap.root";

TString RHCnonswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/RHC/nonswap/FD_RHC_nonswap.root";
TString RHCnueswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/RHC/nueswap/FD_RHC_nueswap.root";
TString RHCtauswap = "/storage/shared/DUNE/OA-inputs/CDR2/FDCAF/RHC/tauswap/FD_RHC_tauswap.root";

double POTperYear = 1.1e21;
double FHCnonswapPOT = 1.62824e24;
double FHCnueswapPOT = 1.64546e24;
double FHCtauswapPOT = 5.18551e24;

double RHCnonswapPOT = 3.27608e+24;
double RHCnueswapPOT = 3.24713e+24;
double RHCtauswapPOT = 8.58955e+24;

double scaleCorrection = 40. / 1.13; //1.13 is kt

void fillArrayWithWeights(int size, std::vector<double>& array, int centralIndex, double centralValue) {
    for (int i = 0; i < size; i++) {
        double value = centralValue + (i - centralIndex);
        array.push_back(value);
    }
}

void SplitMCTree_allSplit_func(bool useRHC, int startEntry, int endEntry, const std::string& outputFileSuffix) {

    TString nonswap = "";
    TString nueswap = "";
    TString tauswap = "";
    double nonswapPOT = 0.;
    double nueswapPOT = 0.;
    double tauswapPOT = 0.;

    if (!useRHC) {
        nonswap = FHCnonswap;
        nueswap = FHCnueswap;
        tauswap = FHCtauswap;

        nonswapPOT = FHCnonswapPOT;
        nueswapPOT = FHCnueswapPOT;
        tauswapPOT = FHCtauswapPOT;
    } else {
        nonswap = RHCnonswap;
        nueswap = RHCnueswap;
        tauswap = RHCtauswap;

        nonswapPOT = RHCnonswapPOT;
        nueswapPOT = RHCnueswapPOT;
        tauswapPOT = RHCtauswapPOT;
    }

    TString filename = "";

    TChain* ch1 = new TChain("cafTree");
    ch1->Add(nonswap);
    ch1->Add(nueswap);
    ch1->Add(tauswap);

    // Create the output file with a unique suffix
    std::string lMode = (useRHC ? "RHC" : "FHC");

    gSystem->mkdir(TString::Format("output_directory"), kTRUE);

    TString outputFileName = TString::Format("./output_directory/MCevents_%s_%s_%d_to_%d.root", lMode.c_str(), outputFileSuffix.c_str(), startEntry, endEntry);
    TFile* outputFile = new TFile(outputFileName, "RECREATE");

    int nuPDG = 0;
    int isCC = 0;
    double Ev = 0.0;
    double POTScaledOscweight = 0.0;
    double POTweight = 0.0;
    int isNC = 0;
    int isRHC = 0;

    int Nonswap = 0;
    int Nueswap = 0;
    int Tauswap = 0;

    ch1->SetBranchStatus("*", true);
    ch1->SetBranchAddress("nuPDG", &nuPDG);
    ch1->SetBranchAddress("isCC", &isCC);
    ch1->SetBranchAddress("Ev", &Ev);

    std::vector<std::string> branchPrefixes = {
        "FSILikeEAvailSmearing",
        "SPPLowQ2Suppression",
        "nuenuebar_xsec_ratio",
        "C12ToAr40_2p2hScaling_nubar",
        "C12ToAr40_2p2hScaling_nu",
        "EbFSLepMomShift",
        "BeRPA_E",
        "BeRPA_D",
        "BeRPA_B",
        "BeRPA_A",
        "NR_nubar_p_NC_3Pi",
        "NR_nubar_p_NC_2Pi",
        "NR_nubar_p_NC_1Pi",
        "NR_nubar_n_NC_3Pi",
        "NR_nubar_n_NC_2Pi",
        "NR_nubar_n_NC_1Pi",
        "NR_nubar_p_CC_3Pi",
        "NR_nubar_p_CC_2Pi",
        "NR_nubar_p_CC_1Pi",
        "NR_nubar_n_CC_3Pi",
        "NR_nubar_n_CC_2Pi",
        "NR_nubar_n_CC_1Pi",
        "NR_nu_p_NC_3Pi",
        "NR_nu_p_NC_2Pi",
        "NR_nu_p_NC_1Pi",
        "NR_nu_n_NC_3Pi",
        "NR_nu_n_NC_2Pi",
        "NR_nu_n_NC_1Pi",
        "NR_nu_np_CC_1Pi",
        "NR_nu_p_CC_3Pi",
        "NR_nu_p_CC_2Pi",
        "NR_nu_n_CC_3Pi",
        "NR_nu_n_CC_2Pi",
        "E2p2h_B_nubar",
        "E2p2h_A_nubar",
        "E2p2h_B_nu",
        "E2p2h_A_nu",
        "MKSPP_ReWeight",
        "Mnv2p2hGaussEnhancement",
        "CCQEPauliSupViaKF",
        "FrPiProd_N",
        "FrAbs_N",
        "FrInel_N",
        "FrElas_N",
        "FrCEx_N",
        "MFP_N",
        "FrPiProd_pi",
        "FrAbs_pi",
        "FrInel_pi",
        "FrElas_pi",
        "FrCEx_pi",
        "MFP_pi",
        "FormZone",
        "CV2uBY",
        "CV1uBY",
        "BhtBY",
        "AhtBY",
        "Theta_Delta2Npi",
        "RDecBR1eta",
        "RDecBR1gamma",
        "MvNCRES",
        "MaNCRES",
        "MvCCRES",
        "MaCCRES",
        "EtaNCEL",
        "MaNCEL",
        "VecFFCCQEshape",
        "MaCCQE"
    };

    std::map<std::string, int> nshiftsMap;
    std::map<std::string, double> cvwgtMap;
    std::map<std::string, Double_t*> wgtMap;  // Pointer for array of weights
    std::map<std::string, std::vector<Double_t>> Xarray;

    for (const auto& prefix : branchPrefixes) {


        std::string prefix_centr = prefix + "_cvwgt";
        std::string prefix_wgt = "wgt_"+prefix;

        //cout<<prefix<<endl;

        if (prefix == "EbFSLepMomShift"){

            //cout<<prefix<<endl;
            prefix_centr = prefix + "_cvvar";
            prefix_wgt = "var_"+prefix;
            //cout<<prefix_centr<<" "<<prefix_wgt<<endl;

        }

        ch1->SetBranchStatus((prefix + "_nshifts").c_str(), 1);
        ch1->SetBranchStatus((prefix_centr).c_str(), 1);
        ch1->SetBranchStatus((prefix_wgt).c_str(), 1);

        nshiftsMap[prefix] = 7;
        cvwgtMap[prefix] = 0.0;
        wgtMap[prefix] = nullptr;  // Initialize to nullptr

        ch1->SetBranchAddress((prefix + "_nshifts").c_str(), &nshiftsMap[prefix]);
        ch1->SetBranchAddress((prefix_centr).c_str(), &cvwgtMap[prefix]);

        // Make sure to clean up old memory before allocating new array
        if (wgtMap[prefix] != nullptr) {
            delete[] wgtMap[prefix];  // Free previously allocated memory
            wgtMap[prefix] = nullptr;  // Reset to nullptr to prevent double deletion
        }

        // Only allocate if there are valid shifts
        if (nshiftsMap[prefix] > 0) {
            wgtMap[prefix] = new Double_t[nshiftsMap[prefix]];  // Allocate new memory for the array
            ch1->SetBranchAddress((prefix_wgt).c_str(), wgtMap[prefix]);  // Set the new branch address
        } else {
            std::cerr << "Warning: nshiftsMap[" << prefix << "] is zero or negative." << std::endl;
        }
    }

    TTree* event_tree = ch1->CloneTree(0);  // Cloning structure only (no events)

    event_tree->SetBranchStatus("*", true);

    event_tree->Branch("POTScaledWeight", &POTScaledOscweight);
    event_tree->Branch("POTweight", &POTweight);
    event_tree->Branch("Nonswap", &Nonswap);
    event_tree->Branch("Nueswap", &Nueswap);
    event_tree->Branch("Tauswap", &Tauswap);
    event_tree->Branch("isNC", &isNC);
    event_tree->Branch("isRHC", &isRHC);

    std::vector<TClonesArray *> response_splines;
    std::vector<TBranch *> branch_vector;

    for (const auto& prefix : branchPrefixes) {

        response_splines.push_back(new TClonesArray("TGraph", 1));
        branch_vector.push_back(event_tree->Branch((prefix+"_TGraph").c_str(), "TClonesArray", response_splines.data()[(int)(response_splines.size() -1)], 256000, 0));

    }

    std::cout<<"response_splines size: " << response_splines.size() << std::endl;
    std::cout<<"branch_vector size:  " << branch_vector.size() << std::endl;

    // Get the number of entries in the chain
    int nentries = ch1->GetEntries();

    // Ensure that the specified range does not exceed the number of entries
    if (endEntry > nentries) {
        endEntry = nentries;
    }

    for (int i = startEntry; i < endEntry; ++i) {
        ch1->GetEntry(i);
        // Reset variables for each entry
        Tauswap = 0;
        Nonswap = 0;
        Nueswap = 0;
        isNC = 0;
        isRHC = 0;

        filename = ch1->GetCurrentFile()->GetName();
        if (isCC) {
            if (filename == nonswap) Nonswap = 1;
            if (filename == nueswap) Nueswap = 1;
            if (filename == tauswap) Tauswap = 1;
        } else {
            if (filename == nonswap) isNC = 1;
            if (filename == nonswap) Nonswap = 1;
        }

        if(useRHC) isRHC = 1;

        POTScaledOscweight = scaleCorrection * POTperYear / nonswapPOT;
        POTweight = POTperYear / nonswapPOT;

        int parCount = -1;
        for (const auto& prefix : branchPrefixes) {
            parCount++;

            response_splines[parCount]->Clear();  // Clear the TClonesArray for the new event
            Xarray[prefix].clear();

            // Check for valid shifts before processing
            if (nshiftsMap[prefix] > 0) {
                fillArrayWithWeights(nshiftsMap[prefix], Xarray[prefix], nshiftsMap[prefix] / 2, cvwgtMap[prefix]);

                TClonesArray& arr_tmp = *response_splines[parCount];
                new(arr_tmp[0]) TGraph();  // Construct new TGraph in TClonesArray
                TGraph* graph = (TGraph*)(arr_tmp[0]);

                if (nshiftsMap[prefix] != 7) std::cout<<prefix<<": "<<nshiftsMap[prefix]<<std::endl;

                // Fill the graph with points
                for (size_t j = 0; j < (size_t)nshiftsMap[prefix]; j++) {
                    double x = Xarray[prefix][j];
                    double y = wgtMap[prefix][j];  // Now properly initialized
                    graph->SetPoint(j, x, y);
                    //std::cout << "Filled TGraph for prefix: " << prefix << ", entry: " << i << " x: " << x << " Weight: " << y << std::endl;
                }
                graph->SetName(prefix.c_str());
                graph->Sort();
            } else {
                std::cerr << "Skipping prefix: " << prefix << " due to zero or negative shifts." << std::endl;
            }
        }

        event_tree->Fill();
    }

    std::cout<<"Cleaning!"<<std::endl;

    outputFile->cd();
    outputFile->WriteObject(event_tree, "event_tree");
    outputFile->Close();

    delete ch1;

}

int main(int argc, char* argv[]) {
    // Check if the correct number of arguments is provided
    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <useRHC> <splitCount> <splitID>" << std::endl;
        return 1;
    }

    // Parse the command-line arguments
    bool useRHC = std::stoi(argv[1]);          // Convert string to integer for useRHC (0 for FHC, 1 for RHC)
    int splitCount = std::stoi(argv[2]);       // Total number of splits (jobs)
    int splitID = std::stoi(argv[3]);          // The split (job) ID

    // Create TChain and dynamically determine the total number of entries
    TChain ch1("cafTree");

    // Add files to the chain
    if (useRHC) {
        ch1.Add(RHCnonswap);
        ch1.Add(RHCnonswap);
        ch1.Add(RHCnonswap);
    } else {
        ch1.Add(FHCnonswap);
        ch1.Add(FHCnonswap);
        ch1.Add(FHCnonswap);    }

    // Get the total number of entries
    Long64_t totalEntries = ch1.GetEntries();
    if (totalEntries <= 0) {
        std::cerr << "Error: No entries found in the TChain!" << std::endl;
        return 1;
    }

    // Dynamically calculate the range of entries this job should process
    Long64_t entriesPerSplit = totalEntries / splitCount;
    Long64_t startEntry = splitID * entriesPerSplit;
    Long64_t endEntry = (splitID == splitCount - 1) ? totalEntries : (startEntry + entriesPerSplit);

    // Print information for debugging purposes
    std::cout << "Processing " << (useRHC ? "RHC" : "FHC") << " entries from " << startEntry << " to " << endEntry << std::endl;
    std::string outputFileSuffix = "job_" + std::to_string(splitID);

    // Call the main function to process the ROOT file
    SplitMCTree_allSplit_func(useRHC, startEntry, endEntry, outputFileSuffix);

    return 0;
}
