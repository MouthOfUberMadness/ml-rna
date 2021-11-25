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

std::vector<std::string> replaceX(const std::vector<std::string> &input,
                                  size_t sizePattern, bool debug = false) {
  std::vector<std::string> output = input;
  for (size_t i = 0; i < sizePattern; i++) {
    std::vector<std::string> temp;
    for (auto pattern : output) {
      if (pattern[i] == 'X') {
        pattern[i] = 'U';
        temp.push_back(pattern);
        pattern[i] = 'C';
        temp.push_back(pattern);
        pattern[i] = 'A';
        temp.push_back(pattern);
        pattern[i] = 'G';
        temp.push_back(pattern);
      } else {
        temp.push_back(pattern);
      }
    }
    output = temp;
    if (debug) {
      std::cout << "updated pattern " << std::endl;
      for (auto pattern : output) {
        std::cout << pattern << " ";
      }
      std::cout << std::endl;
    }
  }
  return output;
}

std::unordered_set<std::string> getRXXRPatterns(bool almost,
                                                bool debug = false) {
  std::vector<std::string> listR = {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"};
  std::vector<std::string> patterns = listR;
  std::vector<std::string> fcsPatterns;

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

  if (debug) {
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  fcsPatterns = replaceX(fcsPatterns, 12, debug);

  if (debug) {
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  std::unordered_set<std::string> fcsPatternsSet;
  for (auto pattern : fcsPatterns) {
    fcsPatternsSet.insert(pattern);
  }
  return fcsPatternsSet;
}

std::unordered_set<std::string> getRRXROrRXRRFcsPatterns(bool almost,
                                                         bool debug = false) {
  std::vector<std::string> listR = {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"};
  std::vector<std::string> patterns = listR;
  std::vector<std::string> fcsPatterns;

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

  if (debug) {
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  fcsPatterns = replaceX(fcsPatterns, 12, debug);

  if (debug) {
    for (auto pattern : fcsPatterns) {
      std::cout << pattern << " ";
    }
    std::cout << std::endl;
  }

  std::unordered_set<std::string> fcsPatternsSet;
  for (auto pattern : fcsPatterns) {
    fcsPatternsSet.insert(pattern);
  }
  return fcsPatternsSet;
}

std::unordered_set<std::string> getFcsPatterns(FCSPatterns fcsPattern,
                                               bool almostExclusive,
                                               bool debug = false) {
  std::unordered_set<std::string> fcsPatternsSet;
  std::unordered_set<std::string> almostFcsPatternsSet;

  if (fcsPattern == FCSPatterns::RXXR) {
    fcsPatternsSet = getRXXRPatterns(false, debug);
    if (almostExclusive)
      almostFcsPatternsSet = getRXXRPatterns(true, debug);
  } else if (fcsPattern == FCSPatterns::RXRROrRRXR) {
    fcsPatternsSet = getRRXROrRXRRFcsPatterns(false, debug);
    if (almostExclusive)
      almostFcsPatternsSet = getRRXROrRXRRFcsPatterns(true, debug);
  }

  if (almostExclusive) {
    for (auto pattern : fcsPatternsSet) {
      almostFcsPatternsSet.erase(pattern);
    }
    return almostFcsPatternsSet;
  } else {
    return fcsPatternsSet;
  }
}
