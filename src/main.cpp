
#include <cstdlib>
#include <iomanip>
#include <iostream>

#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Alphabet/RNA.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/RE08.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>

#include <Bpp/Phyl/Likelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h>

#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/BppOSubstitutionModelFormat.h>

int main(int argc, char *argv[])
{
    std::string nameSeq = "../data/sarbecoviruses.aln";
    // std::string nameSeq = "../data/lysozymeLarge.fasta";
    // std::string nameSeq = "../data/lysozymeLargeGaps.fasta";
    // std::string nameSeq = "../data/lysozymeLarge1Gap.fasta";
    bpp::Fasta Fst;

    auto alpha = std::make_shared<bpp::RNA>();

    auto sites = std::shared_ptr<bpp::AlignedSequenceContainer>(Fst.readAlignment(nameSeq, alpha.get()));

    auto model = std::make_shared<bpp::T92>(alpha.get(), 3., .1);

    auto cast1 = dynamic_cast<bpp::AbstractReversibleNucleotideSubstitutionModel *>(model.get());
    if (cast1 != nullptr)
    {
        std::cout << "Cast1 succeeded" << std::endl;
        auto cast2 = dynamic_cast<bpp::NucleotideReversibleSubstitutionModel *>(model.get());
        if (cast2 != nullptr)
        {
            std::cout << "Cast2 succeeded" << std::endl;
        }
        else
        {
            std::cout << "Cast2 failed" << std::endl;
        }

        auto cast3 = dynamic_cast<bpp::ReversibleSubstitutionModel *>(cast1);
        if (cast3 != nullptr)
        {
            std::cout << "Cast3 succeeded" << std::endl;
            auto cast4 = dynamic_cast<bpp::NucleotideReversibleSubstitutionModel *>(cast3);
            if (cast4 != nullptr)
            {
                std::cout << "Cast4 succeeded" << std::endl;
            }
            else
            {
                std::cout << "Cast4 failed" << std::endl;
            }
        }
        else
        {
            std::cout << "Cast3 failed" << std::endl;
        }
    }
    else
    {
        std::cout << "Cast1 failed" << std::endl;
    }

    // return 0;

    // auto modelBase = std::make_shared<bpp::T92>(alpha.get(), 3., .1);
    // auto model = std::make_shared<bpp::RE08Nucleotide>(reinterpret_cast<bpp::NucleotideReversibleSubstitutionModel *>(modelBase.get()));
    std::cout << "Construct done." << std::endl;
    auto rdist = std::make_shared<bpp::ConstantRateDistribution>();

    auto rootFreqs = std::make_shared<bpp::GCFrequencySet>(alpha.get());

    // Compute ML distance between sequences
    bpp::DistanceEstimation distEstimation(model, rdist, sites.get(), 1, false);

    // Method of clustering
    auto distMethod = std::make_shared<bpp::BioNJ>();
    distMethod->outputPositiveLengths(true);

    // Build tree
    bpp::ParameterList parametersToIgnore;
    bpp::TreeTemplate<bpp::Node> *tree = bpp::OptimizationTools::buildDistanceTree(distEstimation, *distMethod, parametersToIgnore);

    // bool optimizeBrLen = false;
    // std::string param = bpp::OptimizationTools::DISTANCEMETHOD_INIT;
    // double tolerance = 0.000001;
    // unsigned int tlEvalMax = 1000000;
    // unsigned int verbose = 4;
    // bpp::TreeTemplate<bpp::Node> *tree = bpp::OptimizationTools::buildDistanceTree(distEstimation, *distMethod, parametersToIgnore,
    //                                                                                optimizeBrLen,
    //                                                                                param,
    //                                                                                tolerance,
    //                                                                                tlEvalMax,
    //                                                                                0,
    //                                                                                0,
    //                                                                                verbose);

    // Optimize likelihood on this tree

    /* Convert Tree in PhyloTree */
    auto phyloT = bpp::PhyloTreeTools::buildFromTreeTemplate(*tree);
    bpp::ParametrizablePhyloTree parTree(*phyloT);

    /* Create Substitution Process */
    auto subProc = bpp::NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, &parTree, rootFreqs);

    /* Put objects in DataFlow */
    bpp::Context context;
    auto lik = std::make_shared<bpp::LikelihoodCalculationSingleProcess>(context, *sites, *subProc);

    bpp::SingleProcessPhyloLikelihood ntl(context, lik);

    std::cout << std::endl
              << std::endl;
    std::cout << "Initial value " << ntl.getValue() << std::endl;
    bpp::OutputStream *profiler = new bpp::StlOutputStream(new std::ofstream("profile.txt", ios::out));
    bpp::OutputStream *messenger = new bpp::StlOutputStream(new std::ofstream("messages.txt", ios::out));

    bpp::OptimizationTools::optimizeNumericalParameters2(ntl, ntl.getParameters(), 0,
                                                         0.0001, 10000, messenger, profiler, false, false, 1,
                                                         bpp::OptimizationTools::OPTIMIZATION_NEWTON);

    std::cout << "Final value " << ntl.getValue() << std::endl;

    ntl.getParameters().printParameters(cout);

    /* Update Substitution Process parameters */

    subProc->matchParametersValues(ntl.getParameters());

    subProc->getParameters().printParameters(cout);

    bpp::Newick treeWriter;
    treeWriter.writePhyloTree(parTree, "sarbecoviruses_optim.dnd");
    treeWriter.writeTree(*tree, "sarbecoviruses.dnd");
    return 0;
}