/*
Implementation of the DigiChecker class
*/

// Includers from project files
//
#include "DigiChecker.hh"

// Includers from ROOT
//
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TTree.h"

DigiChecker::DigiChecker(const std::string& LegacyDigiFileName, const std::string& NewDigiFileName)
  : LLabel(std::getenv("LegacyLabel") ? std::string(std::getenv("LegacyLabel"))
                                      : "Run3Digitization"),
    NLabel(std::getenv("Run4Label") ? std::string(std::getenv("Run4Label")) : "Run4Digitization"),
    LDFileName(LegacyDigiFileName),
    NDFileName(NewDigiFileName)
{
  PrintLabels();
  LDFile = TFile::Open(LegacyDigiFileName.c_str(), "READ");
  NDFile = TFile::Open(NewDigiFileName.c_str(), "READ");
  OutFile = TFile::Open("ATLSRPCDigiValidation.root", "RECREATE");

  if (!LDFile || LDFile->IsZombie()) {
    std::cerr << "Error while opening file " << LegacyDigiFileName << std::endl;
    std::abort();
  }
  else if (!NDFile || NDFile->IsZombie()) {
    std::cerr << "Error while opening file " << NewDigiFileName << std::endl;
    std::abort();
  }
  else {
    std::cout << "DigiChecker found valid files" << std::endl;
  }
};

DigiChecker::~DigiChecker()
{
  LDFile->Close();
  NDFile->Close();
  OutFile->Close();
  delete LDFile;
  delete NDFile;
  delete OutFile;
};

void DigiChecker::CreateXYMap(const std::string& aFileName,
                              const std::string& HitType /*sim or digi*/)
{
  TTree* tree = nullptr;
  if (aFileName == LDFileName)
    tree = (TTree*)LDFile->Get("MuonHitTest");
  else if (aFileName == NDFileName)
    tree = (TTree*)NDFile->Get("MuonHitTest");
  else
    std::cerr << "DigiChecker::CreateXYSimMap invalid file name: " << aFileName << std::endl;

  std::vector<float>*x = nullptr, *y = nullptr;
  std::vector<char>* stationindex = nullptr;

  if (HitType == "sim") {
    tree->SetBranchAddress("RpcSimHitsGlobPosX", &x);
    tree->SetBranchAddress("RpcSimHitsGlobPosY", &y);
    tree->SetBranchAddress("RpcSimHits_stationIndex", &stationindex);
  }
  else if (HitType == "digi") {
    tree->SetBranchAddress("Digits_RPC_globalPosX", &x);
    tree->SetBranchAddress("Digits_RPC_globalPosY", &y);
    tree->SetBranchAddress("Digits_RPC_stationIndex", &stationindex);
  }
  else {
    std::cerr << "DigiChecker::CreateXYSimMap invalid hittype: " << HitType << std::endl;
  }

  int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kBlack};
  TMultiGraph* mg = new TMultiGraph();
  std::vector<TGraph*> graphs(7, nullptr);
  for (int i = 0; i < 7; i++) {
    graphs[i] = new TGraph();
    graphs[i]->SetMarkerStyle(20);
    graphs[i]->SetMarkerSize(0.5);
    graphs[i]->SetMarkerColor(colors[i]);
  }

  // Loop over tree entries
  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);

    // Loop over hits inside the vectors
    size_t nHits = x->size();
    for (size_t j = 0; j < nHits; j++) {
      int st = static_cast<int>((*stationindex)[j]);  // Convert char to int
      if (st >= 0 && st < 7) {  // Ensure station index is valid
        graphs[st]->SetPoint(graphs[st]->GetN(), (*x)[j], (*y)[j]);
      }
      else {
        std::cerr << "DigiChecker::CreateXYSimMap error with RPC station index: " << st
                  << std::endl;
      }
    }
  }

  // Add graphs to the multigraph
  for (int i = 0; i < 7; i++) {
    if (graphs[i]->GetN() > 0) {  // skip graphs with no points (avoid error)
      mg->Add(graphs[i], "P");  // "P" means points
    }
  }

  // Draw the scatter plot
  std::string canvas_name = "XYMap_" + HitType + "Hit";
  std::string mg_name = "XYMap_" + HitType + "Hit;X [mm];Y [mm]";
  TCanvas* c1 = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
  c1->SetLeftMargin(0.15);
  mg->Draw("A");
  mg->SetTitle(mg_name.c_str());

  // Add legend
  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  for (int i = 0; i < 7; i++) {
    leg->AddEntry(graphs[i], Form("Station index %d", i), "p");
  }
  leg->Draw();

  c1->Update();

  // Write canvas to output file
  OutFile->cd();
  c1->Write();
};

void DigiChecker::CreateYZMap(const std::string& aFileName,
                              const std::string& HitType /*sim or digi*/)
{
  TTree* tree = nullptr;
  if (aFileName == LDFileName)
    tree = (TTree*)LDFile->Get("MuonHitTest");
  else if (aFileName == NDFileName)
    tree = (TTree*)NDFile->Get("MuonHitTest");
  else
    std::cerr << "DigiChecker::CreateYZMap invalid file name: " << aFileName << std::endl;

  std::vector<float>*z = nullptr, *y = nullptr;
  std::vector<char>* etaindex = nullptr;

  if (HitType == "sim") {
    tree->SetBranchAddress("RpcSimHitsGlobPosZ", &z);
    tree->SetBranchAddress("RpcSimHitsGlobPosY", &y);
    tree->SetBranchAddress("RpcSimHits_stationEta", &etaindex);
  }
  else if (HitType == "digi") {
    tree->SetBranchAddress("Digits_RPC_globalPosZ", &z);
    tree->SetBranchAddress("Digits_RPC_globalPosY", &y);
    tree->SetBranchAddress("Digits_RPC_stationEta", &etaindex);
  }
  else {
    std::cerr << "DigiChecker::CreateYZMap invalid hittype: " << HitType << std::endl;
  }

  // Define a color palette for eta indices (-20 to 20)
  const int min_eta = -20, max_eta = 20;
  const int num_eta_bins = max_eta - min_eta + 1;
  int colors[num_eta_bins];

  // Fill the color array with ROOT predefined colors
  int colorList[] = {kRed,   kBlue,     kGreen,    kMagenta,  kCyan,     kOrange,
                     kPink,  kViolet,   kTeal,     kSpring,   kAzure,    kYellow,
                     kBlack, kGray + 1, kGray + 2, kGray + 3, kGray + 4, kGray + 5};
  for (int i = 0; i < num_eta_bins; i++) {
    colors[i] = colorList[i % (sizeof(colorList) / sizeof(colorList[0]))];
  }

  TMultiGraph* mg = new TMultiGraph();

  // Create a vector of TGraph pointers for each eta index (-20 to 20)
  std::vector<TGraph*> graphs(num_eta_bins, nullptr);
  for (int i = 0; i < num_eta_bins; i++) {
    graphs[i] = new TGraph();
    graphs[i]->SetMarkerStyle(20);  // Circular markers
    graphs[i]->SetMarkerColor(colors[i]);
  }

  Long64_t nEntries = tree->GetEntries();
  for (Long64_t i = 0; i < nEntries; i++) {
    tree->GetEntry(i);

    // Loop over hits inside the vectors
    size_t nHits = z->size();
    for (size_t j = 0; j < nHits; j++) {
      int eta = (*etaindex)[j];  // Directly using int
      if (eta >= min_eta && eta <= max_eta) {  // Ensure eta index is valid
        int eta_bin = eta - min_eta;  // Shift to zero-based index
        graphs[eta_bin]->SetPoint(graphs[eta_bin]->GetN(), (*z)[j], (*y)[j]);
      }
    }
  }

  // Add graphs to the multigraph
  for (int i = 0; i < num_eta_bins; i++) {
    if (graphs[i]->GetN() > 0) {  // skip graphs with no points (avoid error)
      mg->Add(graphs[i], "P");  // "P" means points
    }
  }

  // Draw the scatter plot
  std::string canvas_name = "ZYMap_" + HitType + "Hit";
  std::string mg_name = "ZYMap_" + HitType + "Hit;Z;Y";
  TCanvas* c1 = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 800, 600);
  mg->Draw("A");
  mg->SetTitle(mg_name.c_str());
  // Add legend
  TLegend* leg = new TLegend(0.7, 0.3, 0.9, 0.7);
  for (int i = 0; i < num_eta_bins; i++) {
    if (graphs[i]->GetN() > 0) leg->AddEntry(graphs[i], Form("Eta %d", i + min_eta), "p");
  }
  leg->Draw();

  c1->Update();

  // Write canvas to output file
  OutFile->cd();
  c1->Write();
};
