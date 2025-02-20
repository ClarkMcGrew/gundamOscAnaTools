#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TMath.h"
#include "TVector3.h"

using namespace std;
using namespace TMath;

// Global constants
const double DELMSQ_31 = 2.515e-3; // In eV^2
const double LOSC = 1300.;         // In km
const double THETA_23 = 0.859;
const double THETA_13 = 0.150;

const vector<TString> files = {
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nue_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_numu_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nutau_nueselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nue_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_numu_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nutau_nueselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nuebar_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_numubar_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nutaubar_nueselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nuebar_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_numubar_nueselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nutaubar_nueselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nue_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_numu_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nutau_numuselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nue_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_numu_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nutau_numuselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nuebar_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_numubar_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nutaubar_numuselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nuebar_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_numubar_numuselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nutaubar_numuselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nue_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_numu_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnue_x_nutau_ncselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nue_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_numu_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumu_x_nutau_ncselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nuebar_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_numubar_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnuebar_x_nutaubar_ncselec.root",

    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nuebar_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_numubar_ncselec.root",
    "/gpfs/scratch/uyevarouskay/atm/atm_hdnumubar_x_nutaubar_ncselec.root"
};

TString createOutputFileName(const TString& inputFile) {
    TString outputFileName = "/gpfs/scratch/uyevarouskay/atm/gundam_files_fin_v5/";
    TString baseFileName = gSystem->BaseName(inputFile); // Extract file name from path
    outputFileName += baseFileName;
    outputFileName.ReplaceAll(".root", "_gundam.root");
    gSystem->Exec(Form("mkdir -p %s", outputFileName)); // Ensure the output directory exists
    return outputFileName;
}

void prepareGundamMCTree_atm(const TString& inputFile, int id) {
    // Create a TChain for the current file
    TChain* ch1 = new TChain("cafTree");
    ch1->AddFile(inputFile);

    TString outputFileName = createOutputFileName(inputFile);
    TFile* outputFile = new TFile(outputFileName, "RECREATE");
    
    // Set up branch variables.

    Long_t mcID_[10] = {0};            //mc.nu.id                        Long_t
    Long_t mcGenieID_[10] = {0};       //mc.nu.genieId                    Long_t
    Float_t mcOscW_[10] = {0};         //mc.nu.imp_weight = osc_from_x_w      Float_t 
    Float_t mcGenWeight_[10] = {0};    //mc.nu.genweight =xsec_w*nuE_flux     Float_t
    Int_t mcPDG_[10] = {0};          //mc.nu.pdg           Int_t
    Int_t mcPDGorig_[10] = {0};      //mc.nu.pdgorig       Int_t
    Bool_t mcIsCC_[10] = {0};         //mc.nu.iscc          Bool_t
    Int_t mcTargetPDG_[10] = {0};    //mc.nu.targetPDG     Int_t
    Float_t mcEnu_[10] = {0};          //mc.nu.E             Float_t
    Float_t mcVtxX_[10] = {0};         //mc.nu.vtx.x         Float_t
    Float_t mcVtxY_[10] = {0};         //mc.nu.vtx.y         Float_t
    Float_t mcVtxZ_[10] = {0};         //mc.nu.vtx.z         Float_t
    Float_t mcMomX_[10] = {0};         //mc.nu.momentum.x    Float_t
    Float_t mcMomY_[10] = {0};         //mc.nu.momentum.y    Float_t
    Float_t mcMomZ_[10] = {0};         //mc.nu.momentum.z    Float_t
    Float_t mcQ2_[10] = {0};           //mc.nu.Q2            Float_t
    Float_t mcq0_[10] = {0};           //mc.nu.q0            Float_t
   
    Float_t recoVtxX_[10] = {0};    //common.ixn.pandora.vtx.x            Float_t
    Float_t recoVtxY_[10] = {0};    //common.ixn.pandora.vtx.y            Float_t
    Float_t recoVtxZ_[10] = {0};    //common.ixn.pandora.vtx.z            Float_t
   
    //Float_t recoEnuCalo_[10] = {0};    //common.ixn.pandora.Enu.calo         Float_t
    Float_t recoEnuLepCalo_[10] = {0};    //common.ixn.pandora.Enu.lep_calo     Float_t reco mu energy
    //Float_t recoEnuMuRange_[10] = {0};    //common.ixn.pandora.Enu.mu_range     Float_t
    //Float_t recoEnuMuMcs_[10] = {0};    //common.ixn.pandora.Enu.mu_mcs       Float_t
    Float_t recoEnuEcalo_[10] = {0};    //common.ixn.pandora.Enu.e_calo       Float_t reco e energy
   
    Float_t recoEnuCVNnue_[10] = {0};    //common.ixn.pandora.nuhyp.cvn.nue     Float_t
    Float_t recoEnuCVNnumu_[10] = {0};    //common.ixn.pandora.nuhyp.cvn.numu    Float_t
    Float_t recoEnuCVNnc_[10] = {0};    //common.ixn.pandora.nuhyp.cvn.nc      Float_t


    //electron 
    Float_t recoEle_dirX_[10] = {0}; //common.ixn.pandora[0].dir.heshw.x
    Float_t recoEle_dirY_[10] = {0}; //common.ixn.pandora[0].dir.heshw.Y()
    Float_t recoEle_dirZ_[10] = {0}; //common.ixn.pandora[0].dir.heshw.Z() 

    //muon
    Float_t recoMu_dirX_[10] = {0}; //common.ixn.pandora[0].dir.lngtrk.x 
    Float_t recoMu_dirY_[10] = {0}; //common.ixn.pandora[0].dir.lngtrk.Y()
    Float_t recoMu_dirZ_[10] = {0}; //common.ixn.pandora[0].dir.lngtrk.Z()

    // Set up branches for the variables
    ch1->SetBranchAddress("mc.nu.id", &mcID_);
    ch1->SetBranchAddress("mc.nu.genieIdx", &mcGenieID_);
    ch1->SetBranchAddress("mc.nu.imp_weight", &mcOscW_);
    ch1->SetBranchAddress("mc.nu.genweight", &mcGenWeight_);
    ch1->SetBranchAddress("mc.nu.pdg", &mcPDG_);
    ch1->SetBranchAddress("mc.nu.iscc", &mcIsCC_);
    ch1->SetBranchAddress("mc.nu.targetPDG", &mcTargetPDG_);
    ch1->SetBranchAddress("mc.nu.E", &mcEnu_);
    ch1->SetBranchAddress("mc.nu.vtx.x", &mcVtxX_);
    ch1->SetBranchAddress("mc.nu.vtx.y", &mcVtxY_);
    ch1->SetBranchAddress("mc.nu.vtx.z", &mcVtxZ_);
    ch1->SetBranchAddress("mc.nu.momentum.x", &mcMomX_);
    ch1->SetBranchAddress("mc.nu.momentum.y", &mcMomY_);
    ch1->SetBranchAddress("mc.nu.momentum.z", &mcMomZ_);
    ch1->SetBranchAddress("mc.nu.Q2", &mcQ2_);
    ch1->SetBranchAddress("mc.nu.q0", &mcq0_);

    // Additional reconstructed variables
    ch1->SetBranchAddress("common.ixn.pandora.vtx.x", &recoVtxX_);
    ch1->SetBranchAddress("common.ixn.pandora.vtx.y", &recoVtxY_);
    ch1->SetBranchAddress("common.ixn.pandora.vtx.z", &recoVtxZ_);

    ch1->SetBranchAddress("common.ixn.pandora.Enu.lep_calo", &recoEnuLepCalo_);
    ch1->SetBranchAddress("common.ixn.pandora.Enu.e_calo", &recoEnuEcalo_);

    ch1->SetBranchAddress("common.ixn.pandora.nuhyp.cvn.nue", &recoEnuCVNnue_);
    ch1->SetBranchAddress("common.ixn.pandora.nuhyp.cvn.numu", &recoEnuCVNnumu_);
    ch1->SetBranchAddress("common.ixn.pandora.nuhyp.cvn.nc", &recoEnuCVNnc_);

    //electron direction
    ch1->SetBranchAddress("common.ixn.pandora.dir.heshw.x", &recoEle_dirX_);
    ch1->SetBranchAddress("common.ixn.pandora.dir.heshw.y", &recoEle_dirY_);
    ch1->SetBranchAddress("common.ixn.pandora.dir.heshw.z", &recoEle_dirZ_);

    //muon direction ()
    ch1->SetBranchAddress("common.ixn.pandora.dir.lngtrk.x", &recoMu_dirX_); //lngtrk: refers to 
    ch1->SetBranchAddress("common.ixn.pandora.dir.lngtrk.y", &recoMu_dirY_); //a vector associated with the long
    ch1->SetBranchAddress("common.ixn.pandora.dir.lngtrk.z", &recoMu_dirZ_); //track of the neutrino in the detector.


    //ch1->SetBranchAddress("common.ixn.pandora.Enu.calo", &recoEnuCalo_);
    //ch1->SetBranchAddress("common.ixn.pandora.Enu.mu_range", &recoEnuMuRange_);
    //ch1->SetBranchAddress("common.ixn.pandora.Enu.mu_mcs", &recoEnuMuMcs_);


    // Declare new variables for the branches
    Long_t mcID = 0, mcGenieID = 0;
    Float_t mcOscW = 0, mcGenWeight = 0, mcEnu = 0;
    Int_t mcPDG = 0, mcPDGorig = 0, mcIsCC = 0, mcTargetPDG = 0;
    Float_t mcVtxX = 0, mcVtxY = 0, mcVtxZ = 0;
    Float_t mcMomX = 0, mcMomY = 0, mcMomZ = 0;
    Float_t mcQ2 = 0, mcq0 = 0;
    Float_t recoVtxX = 0, recoVtxY = 0, recoVtxZ = 0;
    //Float_t recoEnuCalo = 0, 
    Float_t recoEnuLepCalo = 0, recoEnuEcalo = 0;
    //Float_t recoEnuMuMcs = 0, recoEnuMuRange = 0;
    Float_t recoEnuCVNnue = 0, recoEnuCVNnumu = 0, recoEnuCVNnc = 0;
    Float_t recoEle_dirX = 0, recoEle_dirY = 0, recoEle_dirZ = 0;
    Float_t recoMu_dirX = 0, recoMu_dirY = 0, recoMu_dirZ = 0;

    Int_t sample_idx = -1;
    Float_t recoEle_thetaZ = 0;
    Float_t recoMu_thetaZ = 0;

    Int_t iniNue = 0, iniNuebar = 0, iniNumu = 0, iniNumubar = 0;
    Int_t intrNue = 0, intrNuebar = 0, intrNumu = 0, intrNumubar = 0, intrNutau = 0, intrNutaubar = 0;
    Int_t detNue = 0, detNumu = 0, detNC = 0;

    TTree* event_tree = new TTree();

    //initial neutrino flavor
    event_tree->Branch("iniNue", &iniNue, "iniNue/I");
    event_tree->Branch("iniNuebar", &iniNuebar, "iniNuebar/I");
    event_tree->Branch("iniNumu", &iniNumu, "iniNumu/I");
    event_tree->Branch("iniNumubar", &iniNumubar, "iniNumubar/I");

    //interacted neutrino flavor
    event_tree->Branch("intrNue", &intrNue, "intrNue/I");
    event_tree->Branch("intrNuebar", &intrNuebar, "intrNuebar/I");
    event_tree->Branch("intrNumu", &intrNumu, "intrNumu/I");
    event_tree->Branch("intrNumubar", &intrNumubar, "intrNumubar/I");
    event_tree->Branch("intrNutau", &intrNutau, "intrNutau/I");
    event_tree->Branch("intrNutaubar", &intrNutaubar, "intrNutaubar/I");
    
    //detected neutrino flavor
    event_tree->Branch("detNue", &detNue, "detNue/I");
    event_tree->Branch("detNumu", &detNumu, "detNumu/I");
    event_tree->Branch("detNC", &detNC, "detNC/I");

    event_tree->Branch("sample_idx", &sample_idx, "sample_idx/I");
    event_tree->Branch("recoEle_thetaZ", &recoEle_thetaZ, "recoEle_thetaZ/F");
    event_tree->Branch("recoMu_thetaZ", &recoMu_thetaZ, "recoMu_thetaZ/F");

    event_tree->Branch("mcID", &mcID, "mcID/L");
    event_tree->Branch("mcGenieID", &mcGenieID, "mcGenieID/L");
    event_tree->Branch("mcOscW", &mcOscW, "mcOscW/F");
    event_tree->Branch("mcGenWeight", &mcGenWeight, "mcGenWeight/F");
    event_tree->Branch("mcPDG", &mcPDG, "mcPDG/I");
    event_tree->Branch("mcPDGorig", &mcPDGorig, "mcPDGorig/I");
    event_tree->Branch("mcIsCC", &mcIsCC, "mcIsCC/I");
    event_tree->Branch("mcTargetPDG", &mcTargetPDG, "mcTargetPDG/I");
    event_tree->Branch("mcEnu", &mcEnu, "mcEnu/F");
    event_tree->Branch("mcVtxX", &mcVtxX, "mcVtxX/F");
    event_tree->Branch("mcVtxY", &mcVtxY, "mcVtxY/F");
    event_tree->Branch("mcVtxZ", &mcVtxZ, "mcVtxZ/F");
    event_tree->Branch("mcMomX", &mcMomX, "mcMomX/F");
    event_tree->Branch("mcMomY", &mcMomY, "mcMomY/F");
    event_tree->Branch("mcMomZ", &mcMomZ, "mcMomZ/F");
    event_tree->Branch("mcQ2", &mcQ2, "mcQ2/F");
    event_tree->Branch("mcq0", &mcq0, "mcq0/F");

    // Add branches for the reconstructed variables
    event_tree->Branch("recoVtxX", &recoVtxX, "recoVtxX/F");
    event_tree->Branch("recoVtxY", &recoVtxY, "recoVtxY/F");
    event_tree->Branch("recoVtxZ", &recoVtxZ, "recoVtxZ/F");

    event_tree->Branch("recoEnuLepCalo", &recoEnuLepCalo, "recoEnuLepCalo/F");
    event_tree->Branch("recoEnuEcalo", &recoEnuEcalo, "recoEnuEcalo/F");

    event_tree->Branch("recoEnuCVNnue", &recoEnuCVNnue, "recoEnuCVNnue/F");
    event_tree->Branch("recoEnuCVNnumu", &recoEnuCVNnumu, "recoEnuCVNnumu/F");
    event_tree->Branch("recoEnuCVNnc", &recoEnuCVNnc, "recoEnuCVNnc/F");

    //electron direction
    event_tree->Branch("recoEle_dirX", &recoEle_dirX, "recoEle_dirX/F");
    event_tree->Branch("recoEle_dirY", &recoEle_dirY, "recoEle_dirY/F");
    event_tree->Branch("recoEle_dirZ", &recoEle_dirZ, "recoEle_dirZ/F");

    //muon direction
    event_tree->Branch("recoMu_dirX", &recoMu_dirX, "recoMu_dirX/F");
    event_tree->Branch("recoMu_dirY", &recoMu_dirY, "recoMu_dirY/F");
    event_tree->Branch("recoMu_dirZ", &recoMu_dirZ, "recoMu_dirZ/F");

    //event_tree->Branch("recoEnuCalo", &recoEnuCalo, "recoEnuCalo/F");
    //event_tree->Branch("recoEnuMuRange", &recoEnuMuRange, "recoEnuMuRange/F");
    //event_tree->Branch("recoEnuMuMcs", &recoEnuMuMcs, "recoEnuMuMcs/F");

    Long64_t totalEntries = ch1->GetEntries();
    cout << "Total entries in chain for file " << inputFile << ": " << totalEntries << endl;

    // Loop over the events
    for (Long64_t i = 0; i < totalEntries; ++i) {
        ch1->GetEntry(i);

        sample_idx = id;  // Reset sample index for this event.
    
        mcID = mcID_[0];            //mc.nu.id                        Long_t
        mcGenieID = mcGenieID_[0];       //mc.nu.genieId                    Long_t
        mcOscW = mcOscW_[0];         //mc.nu.imp_weight = osc_from_x_w      
        mcGenWeight = mcGenWeight_[0];    //mc.nu.genweight =xsec_w*nuE_flux     Float_t

        mcPDG = mcPDG_[0];          //mc.nu.pdg           Int_t
        mcPDGorig = mcPDGorig_[0];      //mc.nu.pdgorig       Int_t
        mcIsCC = mcIsCC_[0];         //mc.nu.iscc          Bool_t
        mcTargetPDG = mcTargetPDG_[0];    //mc.nu.targetPDG     Int_t
        mcEnu = mcEnu_[0];          //mc.nu.E             Float_t
        mcVtxX = mcVtxX_[0];         //mc.nu.vtx.x         Float_t
        mcVtxY = mcVtxY_[0];         //mc.nu.vtx.y         Float_t
        mcVtxZ = mcVtxZ_[0];         //mc.nu.vtx.z         Float_t
        mcMomX = mcMomX_[0];         //mc.nu.momentum.x    Float_t
        mcMomY = mcMomY_[0];         //mc.nu.momentum.y    Float_t
        mcMomZ = mcMomZ_[0];         //mc.nu.momentum.z    Float_t
        mcQ2 = mcQ2_[0];           //mc.nu.Q2            Float_t
        mcq0 = mcq0_[0];           //mc.nu.q0            Float_t
       
        recoVtxX = recoVtxX_[0];    //common.ixn.pandora.vtx.x            Float_t
        recoVtxY = recoVtxY_[0];    //common.ixn.pandora.vtx.y            Float_t
        recoVtxZ = recoVtxZ_[0];    //common.ixn.pandora.vtx.z            Float_t
       
        //recoEnuCalo = recoEnuCalo_[0];    //common.ixn.pandora.Enu.calo         Float_t
        recoEnuLepCalo = recoEnuLepCalo_[0];    //common.ixn.pandora.Enu.lep_calo     Float_t
        //recoEnuMuRange = recoEnuMuRange_[0];    //common.ixn.pandora.Enu.mu_range     Float_t
        //recoEnuMuMcs = recoEnuMuMcs_[0];    //common.ixn.pandora.Enu.mu_mcs       Float_t
        recoEnuEcalo = recoEnuEcalo_[0];    //common.ixn.pandora.Enu.e_calo       Float_t
       
        recoEnuCVNnue = recoEnuCVNnue_[0];    //common.ixn.pandora.nuhyp.cvn.nue     Float_t
        recoEnuCVNnumu = recoEnuCVNnumu_[0];    //common.ixn.pandora.nuhyp.cvn.numu    Float_t
        recoEnuCVNnc = recoEnuCVNnc_[0];

        recoMu_dirX = recoMu_dirX_[0];
        recoMu_dirY = recoMu_dirY_[0];
        recoMu_dirZ = recoMu_dirZ_[0];

        recoEle_dirX = recoEle_dirX_[0];
        recoEle_dirY = recoEle_dirY_[0];
        recoEle_dirZ = recoEle_dirZ_[0];

        // from MaCH3 code https://github.com/DUNE/MaCh3_DUNE/blob/develop/samplePDFDUNE/samplePDFDUNEAtm.cpp
        // .Unit() normalizes the vector to make it a unit vector (i.e., it has a magnitude of 1).
        //Zenith angle is calculated as the angle between the momentum vector and the vertical axis (Y-axis)
        //The zenith cosine is just the Y-component of the unit vector p_y

        TVector3 RecoNuEleMomentumVector;
        TVector3 RecoNuMuMomentumVector;

        RecoNuEleMomentumVector = (TVector3(recoEle_dirX,recoEle_dirY,recoEle_dirZ)).Unit();      
        RecoNuMuMomentumVector = (TVector3(recoMu_dirX,recoMu_dirY,recoMu_dirZ)).Unit();      

        recoEle_thetaZ = RecoNuEleMomentumVector.Y();
        recoMu_thetaZ =  RecoNuMuMomentumVector.Y();

        // For debugging, print a few entries.
        //if (i < 5) {
        //    cout << "Entry " << i << ": nuPDG = " << mcPDG
        //         << ", isCC = " << mcIsCC
        //         << ", Ev = " << recoEnuLepCalo << endl;
        //}

        if (sample_idx == 1){
            iniNue = 1;
            intrNue = 1;
            detNue = 1;
        }

        if (sample_idx == 2){
            iniNue = 1;
            intrNumu = 1;
            detNue = 1;
        }

        if (sample_idx == 3){
            iniNue = 1;
            intrNutau = 1;
            detNue = 1;
        }

        if (sample_idx == 4){
            iniNumu = 1;
            intrNue = 1;
            detNue = 1;
        }

        if (sample_idx == 5){
            iniNumu = 1;
            intrNumu = 1;
            detNue = 1;
        }

        if (sample_idx == 6){
            iniNumu = 1;
            intrNutau = 1;
            detNue = 1;
        }

        if (sample_idx == 7){
            iniNuebar = 1;
            intrNuebar = 1;
            detNue = 1;
        }

        if (sample_idx == 8){
            iniNuebar = 1;
            intrNumubar = 1;
            detNue = 1;
        }

        if (sample_idx == 9){
            iniNuebar = 1;
            intrNutaubar = 1;
            detNue = 1;
        }

        if (sample_idx == 10){
            iniNumubar = 1;
            intrNuebar = 1;
            detNue = 1;
        }

        if (sample_idx == 11){
            iniNumubar = 1;
            intrNumubar = 1;
            detNue = 1;
        }

        if (sample_idx == 12){
            iniNumubar = 1;
            intrNutaubar = 1;
            detNue = 1;
        }
   
        //Numu selection

        if (sample_idx == 13){
            iniNue = 1;
            intrNue = 1;
            detNumu = 1;
        }

        if (sample_idx == 14){
            iniNue = 1;
            intrNumu = 1;
            detNumu = 1;
        }

        if (sample_idx == 15){
            iniNue = 1;
            intrNutau = 1;
            detNumu = 1;
        }

        if (sample_idx == 16){
            iniNumu = 1;
            intrNue = 1;
            detNumu = 1;
        }

        if (sample_idx == 17){
            iniNumu = 1;
            intrNumu = 1;
            detNumu = 1;
        }

        if (sample_idx == 18){
            iniNumu = 1;
            intrNutau = 1;
            detNumu = 1;
        }

        if (sample_idx == 19){
            iniNuebar = 1;
            intrNuebar = 1;
            detNumu = 1;
        }

        if (sample_idx == 20){
            iniNuebar = 1;
            intrNumubar = 1;
            detNumu = 1;
        }

        if (sample_idx == 21){
            iniNuebar = 1;
            intrNutaubar = 1;
            detNumu = 1;
        }

        if (sample_idx == 22){
            iniNumubar = 1;
            intrNuebar = 1;
            detNumu = 1;
        }

        if (sample_idx == 23){
            iniNumubar = 1;
            intrNumubar = 1;
            detNumu = 1;
        }

        if (sample_idx == 24){
            iniNumubar = 1;
            intrNutaubar = 1;
            detNumu = 1;
        }

        // NC selection

        if (sample_idx == 25){
            iniNue = 1;
            intrNue = 1;
            detNC = 1;
        }

        if (sample_idx == 26){
            iniNue = 1;
            intrNumu = 1;
            detNC = 1;
        }

        if (sample_idx == 27){
            iniNue = 1;
            intrNutau = 1;
            detNC = 1;
        }

        if (sample_idx == 28){
            iniNumu = 1;
            intrNue = 1;
            detNC = 1;
        }

        if (sample_idx == 29){
            iniNumu = 1;
            intrNumu = 1;
            detNC = 1;
        }

        if (sample_idx == 30){
            iniNumu = 1;
            intrNutau = 1;
            detNC = 1;
        }

        if (sample_idx == 31){
            iniNuebar = 1;
            intrNuebar = 1;
            detNC = 1;
        }

        if (sample_idx == 32){
            iniNuebar = 1;
            intrNumubar = 1;
            detNC = 1;
        }

        if (sample_idx == 33){
            iniNuebar = 1;
            intrNutaubar = 1;
            detNC = 1;
        }

        if (sample_idx == 34){
            iniNumubar = 1;
            intrNuebar = 1;
            detNC = 1;
        }

        if (sample_idx == 35){
            iniNumubar = 1;
            intrNumubar = 1;
            detNC = 1;
        }

        if (sample_idx == 36){
            iniNumubar = 1;
            intrNutaubar = 1;
            detNC = 1;
        }

        event_tree->Fill();
    }

    // Write the output file.
    outputFile->cd();
    event_tree->Write("event_tree");
    outputFile->Close();

    delete ch1;
    delete outputFile;
}

int main(int argc, char* argv[]) {

    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <SLURM_ARRAY_TASK_ID>" << endl;
        return 1;
    }

    int task_id = stoi(argv[1]);  // Get SLURM_ARRAY_TASK_ID
    if (task_id < 0 || task_id >= files.size()) {
        cerr << "Error: SLURM_ARRAY_TASK_ID out of range!" << endl;
        return 1;
    }

    prepareGundamMCTree_atm(files[task_id], task_id + 1);  // Process corresponding file based on task_id
    return 0;

}
