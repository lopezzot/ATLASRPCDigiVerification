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
  const int stationidx = 11;

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

  // int colors[] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kBlack};
  int colors[] = {kRed,   kBlue,   kGreen,  kMagenta, kCyan,  kOrange,
                  kBlack, kViolet, kYellow, kGray,    kSpring};
  TMultiGraph* mg = new TMultiGraph();
  std::vector<TGraph*> graphs(stationidx, nullptr);
  for (int i = 0; i < stationidx; i++) {
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
      if (st >= 0 && st < stationidx) {  // Ensure station index is valid
        graphs[st]->SetPoint(graphs[st]->GetN(), (*x)[j], (*y)[j]);
      }
      else {
        std::cerr << "DigiChecker::CreateXYSimMap error with RPC station index: " << st
                  << std::endl;
      }
    }
  }

  // Add graphs to the multigraph
  for (int i = 0; i < stationidx; i++) {
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
  for (int i = 0; i < stationidx; i++) {
    leg->AddEntry(graphs[i], Form("Station index %d", i), "p");
  }
  leg->Draw();

  c1->Update();

  // Write canvas to output file
  OutFile->cd();
  c1->Write();
};
