#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SequencePositionIterators.h>
#include <Bpp/Seq/SequenceWalker.h>

#include <Bpp/Phyl/Tree/PhyloTreeTools.h>

#include <Bpp/Phyl/Io/Newick.h>

#include "analysis.h"

struct StructuredHistogram {
  StructuredHistogram(size_t s, size_t sS1, size_t eS1, size_t sS2, size_t eS2,
                      size_t e) {
    start = s;
    startS1 = sS1;
    endS1 = eS1;
    startS2 = sS2;
    endS2 = eS2;
    end = e;
    data.resize(5);
    std::fill(data.begin(), data.end(), 0);
  }
  void add(size_t idx) {
    if ((start <= idx) && (idx < startS1)) {
      data[0]++;
    } else if ((startS1 <= idx) && (idx <= endS1)) {
      data[1]++;
    } else if ((endS1 < idx) && (idx <= startS2)) {
      data[2]++;
    } else if ((startS2 < idx) && (idx <= endS2)) {
      data[3]++;
    } else if ((endS2 < idx) && (idx <= end)) {
      data[4]++;
    }
  }
  void add(const StructuredHistogram &histogram) {
    for (size_t i = 0; i < histogram.data.size(); i++) {
      data[i] += histogram.data[i];
    }
  }
  void reset() { std::fill(data.begin(), data.end(), 0); }

  void printNormalized(const StructuredHistogram &histogram) {
    std::vector<double> normalized(data.size());
    double sum = 0;
    for (size_t i = 0; i < normalized.size(); i++) {
      // TODO: check division by 0
      normalized[i] = ((double)data[i]) / ((double)histogram.data[i]);
      sum += normalized[i];
    }

    std::cout << "linear probabilities | " << normalized[0] << " |S1 "
              << normalized[1] << " |FCS " << normalized[2] << " |S2 "
              << normalized[3] << " | " << normalized[4] << " | ";

    for (size_t i = 0; i < normalized.size(); i++) {
      normalized[i] /= sum;
    }

    std::cout << "  --  " << sum << " * | " << normalized[0] << " |S1 "
              << normalized[1] << " |FCS " << normalized[2] << " |S2 "
              << normalized[3] << " | " << normalized[4] << " | ";
    std::cout << std::endl;
  }

  void print() {
    std::cout << "histogram | " << data[0] << " |S1 " << data[1] << " |FCS "
              << data[2] << " | "
              << " |S2 " << data[3] << " | " << data[4] << " | " << std::endl;
  }

  size_t start;
  size_t startS1;
  size_t endS1;
  size_t startS2;
  size_t endS2;
  size_t end;
  std::vector<size_t> data;
};

struct WindowHistogram {
  WindowHistogram(size_t size, size_t windowRadius_) {
    windowRadius = windowRadius_;
    data.resize(size);
    std::fill(data.begin(), data.end(), 0);
  }
  void add(size_t idx) {
    for (size_t i = idx - windowRadius; i <= (idx + windowRadius); i++) {
      data[i]++;
    }
  }
  void add(const WindowHistogram &histogram) {
    for (size_t i = 0; i < histogram.data.size(); i++) {
      data[i] += histogram.data[i];
    }
  }

  void reset() { std::fill(data.begin(), data.end(), 0); }

  void print(size_t anchor, size_t numSamples) {
    size_t size = data.size();
    size_t delta = size / numSamples;
    size_t steps = anchor / delta;
    size_t start = anchor - steps * delta;

    size_t cursor = start;
    while (cursor < size) {
      if (cursor == anchor)
        std::cout << "[" << data[cursor] << "] ";
      else
        std::cout << data[cursor] << " ";
      cursor += delta;
    }
  }

  size_t windowRadius;
  std::vector<size_t> data;
};

void toAlignedMapping(const bpp::AlignedSequenceContainer &alignedSequences,
                      size_t seqId, std::vector<size_t> &mapping) {
  auto alpha = std::make_shared<bpp::RNA>();
  size_t numSites = alignedSequences.getNumberOfSites();
  size_t cursor = 0;

  mapping.clear();

  while (cursor < numSites) {
    bpp::Site site = alignedSequences.getSite(cursor);
    auto nucleotide = alpha.get()->intToChar(site[seqId]);

    if (nucleotide != "N") {
      mapping.push_back(cursor);
    }
    cursor++;
  }
}

size_t findStartingCodon(const bpp::AlignedSequenceContainer &alignedSequences,
                         size_t seqId) {
  auto alpha = std::make_shared<bpp::RNA>();
  size_t numSites = alignedSequences.getNumberOfSites();
  std::string aCodon = "AUGGAGAG";

  for (size_t i = 0; i < (numSites - aCodon.size()); i++) {
    std::string codon;
    for (size_t j = 0; j < aCodon.size(); j++) {
      bpp::Site site = alignedSequences.getSite(i + j);
      auto nucleotide = alpha.get()->intToChar(site[seqId]);
      codon.append(nucleotide);
    }
    if (codon == aCodon) {
      return i;
    }
  }
  return 0;
}

bool matchesWildcard(const std::string &s1, const std::string &s2) {
  int size = std::min(s1.size(), s2.size());
  for (int i = 0; i < size; i++) {
    if ((s1[i] != s2[i]) && (s1[i] != 'X') && (s2[i] != 'X')) {
      return false;
    }
  }
  return true;
}

template <typename H, typename P>
bool processSequence(const bpp::AlignedSequenceContainer &alignedSequences,
                     const P &patternList, H &histogram, H &supportHistogram,
                     size_t seqId, size_t readingFrame) {
  bool found = false;
  auto alpha = std::make_shared<bpp::RNA>();
  size_t cursor = readingFrame;
  size_t numSites = alignedSequences.getNumberOfSites();
  size_t fullStrideSize = 3; // Codon sized stride
  size_t fullBufferSize = patternList.front().size();

  // numSites = 400;
  while (cursor < numSites) {

    std::string buffer;
    bool isBufferFull = false;
    size_t internalCursor = cursor;

    // Fill the buffer starting from "cursor".
    while (!isBufferFull && (internalCursor < numSites)) {
      bpp::Site site = alignedSequences.getSite(internalCursor);
      auto nucleotide = alpha.get()->intToChar(site[seqId]);

      if (nucleotide != "N")
        buffer.append(nucleotide);

      isBufferFull = (buffer.size() == fullBufferSize);

      internalCursor++;
    }

    if (isBufferFull) {
      // std::cout << "buffer = " << buffer << " at " << cursor << std::endl;

      // Check if the buffer match the pattern list.
      auto it = std::find_if(patternList.begin(), patternList.end(),
                             [&buffer](const std::string &s) {
                               return matchesWildcard(s, buffer);
                             });
      found = it != patternList.end();
      if (found) {
        // std::cout << "buffer = " << buffer << " match " << *it << " at "
        //           << cursor << std::endl;
        histogram.add(cursor);
      }
    }

    // Take a stride of size 3 ignoring gaps ("N")
    size_t stride = 0;
    bool isStrideFull = false;
    while (!isStrideFull && (cursor < numSites)) {
      bpp::Site site = alignedSequences.getSite(cursor);
      auto nucleotide = alpha.get()->intToChar(site[seqId]);

      if (nucleotide != "N") {
        supportHistogram.add(cursor);
        stride++;
      }

      isStrideFull = (stride == fullStrideSize);

      cursor++;
    }
  }

  return found;
}

void structureAnalysis(const bpp::AlignedSequenceContainer &sites,
                       size_t startCodon, size_t startS1, size_t endS1,
                       size_t startS2, size_t endS2) {
  size_t numSequences = sites.getNumberOfSequences();
  size_t numSites = sites.getNumberOfSites();
  size_t end = numSites - 1;
  std::cout << "start = " << 0 << " start S1 = " << startS1
            << " end S1 = " << endS1 << " start S2 = " << startS2
            << " end S2 = " << endS2 << " end = " << end << std::endl;

  // size_t window = 20;
  // WindowHistogram gHistogram(numSites, window);
  // WindowHistogram histogram(numSites, window);
  StructuredHistogram gHistogram(0, startS1, endS1, startS2, endS2, end);
  StructuredHistogram gSupportHistogram(0, startS1, endS1, startS2, endS2, end);
  StructuredHistogram histogram(0, startS1, endS1, startS2, endS2, end);
  StructuredHistogram supportHistogram(0, startS1, endS1, startS2, endS2, end);

  std::vector<std::string> listR = {"CGU", "CGC", "CGA", "CGG", "AGA", "AGG"};
  std::vector<std::string> listA = {"UUA", "UUG", "CUU", "CUC", "CUA", "CUG"};
  std::vector<std::string> patternList = listR;

  std::vector<std::string> almostPatternList;
  for (auto pattern : patternList) {
    for (size_t i = 0; i < 3; i++) {
      auto almostPattern = pattern;
      almostPattern[i] = 'X';
      almostPatternList.push_back(almostPattern);
    }
  }

  std::vector<std::string> duetPatternList;
  for (size_t i = 0; i < listR.size(); i++)
    for (size_t j = 0; j < listR.size(); j++) {
      duetPatternList.push_back(patternList[i] + std::string("XXX") +
                                patternList[j]);
    }

  std::vector<std::string> almostDuetPatternList;
  for (auto almostPattern : almostPatternList)
    for (auto pattern : patternList) {
      almostDuetPatternList.push_back(almostPattern + std::string("XXX") +
                                      pattern);
      almostDuetPatternList.push_back(pattern + std::string("XXX") +
                                      almostPattern);
    }

  // std::vector<std::string> finalPatternList = almostPatternList;
  // std::cout << "almost patterns" << std::endl;
  // for (auto element : almostPatternList) {
  //   std::cout << element << " ";
  // }
  // std::cout << std::endl;

  // std::vector<std::string> finalPatternList = duetPatternList;
  // std::cout << "duet patterns" << std::endl;
  // for (auto element : duetPatternList) {
  //   std::cout << element << " ";
  // }
  // std::cout << std::endl;

  std::vector<std::string> finalPatternList = almostDuetPatternList;
  std::cout << "almost duet patterns" << std::endl;
  // for (auto element : almostDuetPatternList) {
  //   std::cout << element << " ";
  // }
  // std::cout << std::endl;

  for (size_t i = 0; i < numSequences; i++) {
    bool debug = false;
    auto name = sites.getSequencesNames()[i];
    if (name.rfind("MT835139.1", 0) ==
        0) { // pos=0 limits the search to the prefix
      // s starts with prefix
      debug = true;
    }

    if (debug)
      std::cout << "name = " << name << std::endl;

    histogram.reset();
    supportHistogram.reset();

    size_t readingFrame = startCodon;

    if (debug)
      std::cout << "reading frame = " << readingFrame << std::endl;

    processSequence(sites, finalPatternList, histogram, supportHistogram, i,
                    readingFrame);

    if (i == 0 && false) {
      histogram.print();
      histogram.printNormalized(supportHistogram);
      supportHistogram.print();
    }

    gHistogram.add(histogram);
    gSupportHistogram.add(supportHistogram);
  }

  // gHistogram.print();
  gHistogram.printNormalized(gSupportHistogram);
  // gSupportHistogram.print();
  std::cout << std::endl;
}

void analysePhylogeneticTree(const std::string &filename, size_t startCodon,
                             size_t startS1, size_t endS1, size_t startS2,
                             size_t endS2) {
  size_t lastindex = filename.find_last_of(".");
  std::string rawname = filename.substr(0, lastindex);
  std::string sequencesName = rawname + std::string(".aln");

  bpp::Newick treeIO;
  auto phyloTree = treeIO.readPhyloTree(filename);

  bpp::Fasta Fst;
  auto alpha = std::make_shared<bpp::RNA>();
  auto sites = std::shared_ptr<bpp::AlignedSequenceContainer>(
      Fst.readAlignment(sequencesName, alpha.get()));

  structureAnalysis(*sites, startCodon, startS1, endS1, startS2, endS2);
}