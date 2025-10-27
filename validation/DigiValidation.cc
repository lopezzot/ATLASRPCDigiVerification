/*
Description: Code to validate ATLAS Muon RPC Digitization
tools and compare the Run3 legacy validation tools with new
Run4 one.
Author: Lorenzo Pezzotti
Usage: Assumes two root input files, one produced with the
legacy digitization code and one with the new digitization code.
*/

// Includers from ROOT
//
#include "TFile.h"

// Includers from STL
//
#include <iostream>

// Includers from project files
//
#include "DigiChecker.hh"

int main(int argc, char* argv[])
{
  if (argc < 3) {  // Ensure at least two arguments are provided
    std::cerr << "Usage: " << argv[0] << " <LegacyDigiFileName> <NewDigiFileName>" << std::endl;
    return 1;
  }
  std::string LDFileName = argv[1];  // Legacy Digi File Name
  std::string NDFileName = argv[2];  // New Digi File Name
  DigiChecker theDigiChecker(LDFileName, NDFileName);
  theDigiChecker.CreateXYMap(NDFileName, "sim");
  theDigiChecker.CreateXYMap(NDFileName, "digi");

  std::cout << "bye bye" << std::endl;
  return 0;
}
