/*
Definition of the DigiChecker class
*/

// Includers from ROOT
//
#include "TFile.h"

// Includers from STL
//
#include <cstdlib>
#include <iostream>
#include <string>

class DigiChecker
{
  public:
    // RAII implementation, i.e. input files
    // lifecycle is bound to the lifetime of this class
    //
    DigiChecker() = delete;
    DigiChecker(const std::string& LegacyDigiFileName, const std::string& NewDigiFileName);
    ~DigiChecker();
    void PrintLabels() const;
    void CreateXYMap(const std::string& aFileName, const std::string& HitType /*sim or digi*/);

  private:
    TFile* LDFile = nullptr;  // LegacyDigiFile
    TFile* NDFile = nullptr;  // NewDigiFile
    TFile* OutFile = nullptr;  // Output file
    std::string LLabel;  // Label for legacy digitization plots
    std::string NLabel;  // Label for run4 digitization plots
    std::string LDFileName;
    std::string NDFileName;
};

inline void DigiChecker::PrintLabels() const
{
  std::cout << "DigiChecker using plot labels: " << LLabel << " " << NLabel << std::endl;
};
