#pragma once

#include <string>
#include <unordered_set>
#include <vector>

// From @daoyu15
// RXXR is the minimal requirement of FCS
// RXR or RR or RXXXR will not be cleaved by furin
// FCS found in proteins are normally RXRR or RRXR
// Some viral FCS are RRXRR
// R will not be cleaved by furin
// Occasionally for the RRXRR sites the second and third R can be replaced by
// other residues, basic is the norm

enum class FCSPatterns {
  RXXR,
  RXRROrRRXR,
  RRXRR,
};

std::unordered_set<std::string>
getFcsPatterns(FCSPatterns fcsPattern, bool almost, bool debug = false) {

  std::vector<std::string> listR = {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"};
  std::vector<std::string> patterns;
  std::vector<std::string> almostPatterns;
  std::vector<std::string> fcsPatterns;
  std::unordered_set<std::string> fcsPatternsSet;

  patterns = listR;

  if (fcsPattern == FCSPatterns::RXXR) {
    for (size_t i = 0; i < patterns.size(); i++)
      for (size_t j = 0; j < patterns.size(); j++) {
        if (almost) {
          for (size_t l = 0; l < 3; l++) {
            auto almostPattern =
                patterns[i] + std::string("XXXXXX") + patterns[j];
            almostPattern[l] = 'X';
            fcsPatterns.push_back(almostPattern);
          }
          for (size_t l = 0; l < 3; l++) {
            auto almostPattern =
                patterns[i] + std::string("XXXXXX") + patterns[j];
            almostPattern[9 + l] = 'X';
            fcsPatterns.push_back(almostPattern);
          }
        } else {
          fcsPatterns.push_back(patterns[i] + std::string("XXXXXX") +
                                patterns[j]);
        }
      }
  } else if (fcsPattern == FCSPatterns::RXRROrRRXR) {
    for (size_t i = 0; i < patterns.size(); i++)
      for (size_t j = 0; j < patterns.size(); j++) {
        for (size_t k = 0; k < patterns.size(); k++) {
          auto patternLeft =
              patterns[i] + std::string("XXX") + patterns[j] + patterns[k];
          auto patternRight =
              patterns[i] + patterns[j] + std::string("XXX") + patterns[k];

          if (almost) {
            for (size_t l = 0; l < 3; l++) {
              auto almostPatternLeft = patternLeft;
              auto almostPatternRight = patternRight;
              almostPatternLeft[l] = 'X';
              almostPatternRight[9 + l] = 'X';
              fcsPatterns.push_back(almostPatternLeft);
              fcsPatterns.push_back(almostPatternRight);
            }
            for (size_t l = 0; l < 6; l++) {
              auto almostPatternLeft = patternLeft;
              auto almostPatternRight = patternRight;
              almostPatternLeft[6 + l] = 'X';
              almostPatternRight[l] = 'X';
              fcsPatterns.push_back(almostPatternLeft);
              fcsPatterns.push_back(almostPatternRight);
            }
          } else {
            fcsPatterns.push_back(patternLeft);
            fcsPatterns.push_back(patternRight);
          }
        }
      }
  } else if (fcsPattern == FCSPatterns::RRXRR) {
  }

  if (debug) {
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  if (fcsPattern == FCSPatterns::RXXR ||
      fcsPattern == FCSPatterns::RXRROrRRXR) {
    for (size_t i = 0; i < 12; i++) {
      std::vector<std::string> tempFcsPatterns;
      for (auto pattern : fcsPatterns) {
        if (pattern[i] == 'X') {
          pattern[i] = 'U';
          tempFcsPatterns.push_back(pattern);
          pattern[i] = 'C';
          tempFcsPatterns.push_back(pattern);
          pattern[i] = 'A';
          tempFcsPatterns.push_back(pattern);
          pattern[i] = 'G';
          tempFcsPatterns.push_back(pattern);
        } else {
          tempFcsPatterns.push_back(pattern);
        }
      }
      fcsPatterns = tempFcsPatterns;
      if (debug) {
        std::cout << "updated pattern " << std::endl;
        for (auto pattern : fcsPatterns) {
          std::cout << pattern << " ";
        }
        std::cout << std::endl;
      }
    }
  } else if (fcsPattern == FCSPatterns::RRXRR) {
  }

  if (debug) {
    std::cout << "final pattern " << std::endl;
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  for (auto pattern : fcsPatterns) {
    fcsPatternsSet.insert(pattern);
  }

  return fcsPatternsSet;
}
