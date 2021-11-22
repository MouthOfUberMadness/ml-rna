#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SequencePositionIterators.h>
#include <Bpp/Seq/SequenceWalker.h>

#include <Bpp/Phyl/Tree/PhyloTreeTools.h>

#include <Bpp/Phyl/Io/Newick.h>

#include "analysis.h"

void processSequence(const bpp::AlignedSequenceContainer &alignedSequences,
                     size_t seqId) {
  //   auto walker = bpp::SequenceWalker(sequence);

  auto alpha = std::make_shared<bpp::RNA>();
  std::string codon;
  size_t cursor = 0;
  size_t cursorCodon = 0;
  bool isTheEnd = false;

  size_t numSites = alignedSequences.getNumberOfSites();
  size_t numBins = 20;

  std::vector<size_t> histogram(numBins);
  std::string aCodon = "CGG";
//   std::string aCodon = "UCC";
  while (cursor < numSites) {

    bpp::Site site = alignedSequences.getSite(cursor);
    if (seqId < site.size()) {
      auto nucleotide = alpha.get()->intToChar(site[seqId]);
      //   if (cursor < 100)
      //     std::cout << "cursor " << cursor << " " << nucleotide << std::endl;
      if (nucleotide != "N")
        codon.append(nucleotide);
    }

    if (codon.size() == 3) {
      if (codon == aCodon) {
        size_t bin = std::floor(((double)cursor / (double)numSites) * numBins);
        histogram[bin]++;
      }

      codon.clear();
    }

    cursor++;
  }

  for (size_t bin : histogram) {
    std::cout << bin << std::endl;
  }
}

void analysePhylogeneticTree(const std::string &filename) {
  size_t lastindex = filename.find_last_of(".");
  std::string rawname = filename.substr(0, lastindex);
  std::string sequencesName = rawname + std::string(".aln");

  bpp::Newick treeIO;
  auto phyloTree = treeIO.readPhyloTree(filename);

  bpp::Fasta Fst;
  auto alpha = std::make_shared<bpp::RNA>();
  auto sites = std::shared_ptr<bpp::AlignedSequenceContainer>(
      Fst.readAlignment(sequencesName, alpha.get()));

  processSequence(*sites, 0);
}