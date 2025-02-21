/*
Definition of the DigiChecker class
*/

// Includers from ROOT
//
#include "TFile.h"

// Includers from STL
//
#include <iostream>

class DigiChecker
{
  public:
    // RAII implementation, i.e. input files
    // lifecycle is bound to the lifetime of this class
    //
    DigiChecker() = delete;
    DigiChecker(const std::string& LegacyDigiFileName, const std::string& NewDigiFileName);
    ~DigiChecker();

  private:
    TFile* LDFile = nullptr;  // LegacyDigiFile
    TFile* NDFile = nullptr;  // NewDigiFile
    TFile* OutFile = nullptr;  // Output file
};
