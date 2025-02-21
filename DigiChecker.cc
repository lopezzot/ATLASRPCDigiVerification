/*
Implementation of the DigiChecker class
*/

// Includers from project files
//
#include "DigiChecker.hh"

DigiChecker::DigiChecker(const std::string& LegacyDigiFileName, const std::string& NewDigiFileName)
  : LLabel(std::getenv("LegacyLabel") ? std::string(std::getenv("LegacyLabel"))
                                      : "Run3Digitization"),
    NLabel(std::getenv("Run4Label") ? std::string(std::getenv("Run4Label")) : "Run4Digitization")
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
