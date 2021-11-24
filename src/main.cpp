
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/Fasta.h>

#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RE08.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>

#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>

#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>
#include <Bpp/Phyl/Io/Newick.h>

#include "analysis.h"

void buildPhylogeneticTreeFromAlignement(const std::string &filename) {
  std::string nameSeq = filename;
  bpp::Fasta Fst;

  auto alpha = std::make_shared<bpp::RNA>();

  auto sites = std::shared_ptr<bpp::AlignedSequenceContainer>(
      Fst.readAlignment(nameSeq, alpha.get()));

  for (size_t i = 0; i < sites->getNumberOfSequences(); i++) {
    bool debug = false;
    auto name = sites->getSequencesNames()[i];
    std::cout << "name = " << name << std::endl;
  }

  auto model = std::make_shared<bpp::T92>(alpha.get(), 3., .1);

  auto rdist = std::make_shared<bpp::ConstantRateDistribution>();

  auto rootFreqs = std::make_shared<bpp::GCFrequencySet>(alpha.get());

  // Compute ML distance between sequences
  bpp::DistanceEstimation distEstimation(model, rdist, sites.get(), 1, false);

  // Method of clustering
  auto distMethod = std::make_shared<bpp::BioNJ>();
  distMethod->outputPositiveLengths(true);

  // Build tree
  bpp::ParameterList parametersToIgnore;
  bool optimizeBrLen = false;
  std::string param = bpp::OptimizationTools::DISTANCEMETHOD_INIT;
  double tolerance = 0.000001;
  unsigned int tlEvalMax = 1000000;
  unsigned int verbose = 2;
  bpp::TreeTemplate<bpp::Node> *tree =
      bpp::OptimizationTools::buildDistanceTree(
          distEstimation, *distMethod, parametersToIgnore, optimizeBrLen, param,
          tolerance, tlEvalMax, 0, 0, verbose);

  // Optimize likelihood on this tree

  /* Convert Tree in PhyloTree */
  auto phyloT = bpp::PhyloTreeTools::buildFromTreeTemplate(*tree);
  bpp::ParametrizablePhyloTree parTree(*phyloT);

  /* Create Substitution Process */
  auto subProc = bpp::NonHomogeneousSubstitutionProcess::
      createHomogeneousSubstitutionProcess(model, rdist, &parTree, rootFreqs);

  /* Put objects in DataFlow */
  bpp::Context context;
  auto lik = std::make_shared<bpp::LikelihoodCalculationSingleProcess>(
      context, *sites, *subProc);

  bpp::SingleProcessPhyloLikelihood ntl(context, lik);

  std::cout << std::endl << std::endl;
  std::cout << "Initial value " << ntl.getValue() << std::endl;
  bpp::OutputStream *profiler =
      new bpp::StlOutputStream(new std::ofstream("profile.txt", ios::out));
  bpp::OutputStream *messenger =
      new bpp::StlOutputStream(new std::ofstream("messages.txt", ios::out));

  bpp::OptimizationTools::optimizeNumericalParameters2(
      ntl, ntl.getParameters(), 0, 0.0001, 10000, messenger, profiler, false,
      false, 1, bpp::OptimizationTools::OPTIMIZATION_NEWTON);

  std::cout << "Final value " << ntl.getValue() << std::endl;

  ntl.getParameters().printParameters(cout);

  /* Update Substitution Process parameters */

  subProc->matchParametersValues(ntl.getParameters());

  subProc->getParameters().printParameters(cout);

  bpp::Newick treeWriter;
  size_t lastindex = filename.find_last_of(".");
  std::string rawname = filename.substr(0, lastindex);
  treeWriter.writePhyloTree(parTree, rawname + std::string(".dnd"));
  treeWriter.writeTree(*tree, rawname + std::string("_unoptim.dnd"));
}

int main(int argc, char *argv[]) {
  std::string filename;
  std::string option;
  std::string readingFrame;
  std::string startS1;
  std::string endS1;
  std::string startS2;
  std::string endS2;
  if (argc == 8) {
    readingFrame = argv[1];
    startS1 = argv[2];
    endS1 = argv[3];
    startS2 = argv[4];
    endS2 = argv[5];
    option = argv[6];
    filename = argv[7];
  } else {
    return 0;
  }

  if (option == "--build-phyl") {
    buildPhylogeneticTreeFromAlignement(filename);
  } else if (option == "--analyse-phyl") {
    analysePhylogeneticTree(filename, std::stoi(readingFrame),
                            std::stoi(startS1), std::stoi(endS1),
                            std::stoi(startS2), std::stoi(endS2));
  }
}