void run(std::string filename1, std::string DialNames, TCanvas *canv) 
{
    // Open the file and get the tree
    TFile *file = new TFile(filename1.c_str());
    TTree *tree = (TTree*)file->Get("event_tree");

    // Use the TClonesArray for each branch
    TClonesArray *SPLINE_1 = new TClonesArray(DialNames.c_str());
    tree->SetBranchAddress(DialNames.c_str(), &SPLINE_1);
        
    int nentries = tree->GetEntries();
    std::cout << "Processing " << nentries << " entries for branch: " << DialNames << std::endl;

    // Reuse the graph and spline objects
    //TGraph *graphSpline;
    //TSpline3 *Tspline = nullptr;

    for (int i = 9; i < 10; ++i)  // Process 3 events only
    {

            // Reuse the graph and spline objects
        TGraph *graphSpline;
        TSpline3 *Tspline = nullptr;

        tree->GetEntry(i);
        //cout<<SPLINE_1->GetEntries()<<endl;
        graphSpline = (TGraph*)((SPLINE_1)->At(0));

        // Skip if the graph is empty or has only one knot
        if (!graphSpline || graphSpline->GetN() <= 1) {
            std::cerr << "Skipping event " << i << ": graph has only one knot or is null." << std::endl;
            continue;
        }

        // Set canvas pad and grid
        canv->cd();
        gPad->SetGrid();
        gPad->SetLeftMargin(0.2);  // Set margin to fit the Y-axis label

        // Configure graph appearance
        graphSpline->SetLineStyle(kSolid);
        graphSpline->SetLineColor(kBlue);
        graphSpline->SetMarkerStyle(5);
        graphSpline->SetMarkerColor(kRed);
        graphSpline->SetTitle(Form("%s - Event %i", DialNames.c_str(), i));
        graphSpline->GetYaxis()->SetTitle("weight");
        graphSpline->GetXaxis()->SetTitle("parameter variation");

        // Create and configure the spline
        if (Tspline) delete Tspline;  // Avoid memory leaks by deleting the old spline
        Tspline = new TSpline3(Form("%s - Event %i", DialNames.c_str(), i), graphSpline);
        Tspline->SetLineColor(kRed);

        // Draw the graph and spline
        graphSpline->Draw("ALP");
        Tspline->Draw("SAME");

        // No need for frequent canvas updates; only update if necessary
        canv->Update();
    canv->Print("DialsDUNE.pdf");

    }


    // Cleanup
    delete SPLINE_1;
    delete file;
}

void drawSplines(std::string filename1) 
{
    // List of branch prefixes
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

    // Create a canvas and open the PDF file once
    TCanvas *canv = new TCanvas("canv", "canv", 341, 424);
    canv->Print("DialsDUNE.pdf[");  // Open the PDF file

    // Loop over the branch prefixes and run the function
    for (const auto& prefix : branchPrefixes) {
        run(filename1, prefix + "_TGraph", canv);
    }

    // Close the PDF file
    canv->Print("DialsDUNE.pdf]");
    delete canv;  // Clean up canvas
}
