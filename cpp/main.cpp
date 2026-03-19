//
// Created by martinalavanya on 08/10/2025.
//

#include <cassert>
#include <iostream>

#include "BedGraphRow.h"
#include "Parser.h"
#include "Averager.h"
#include "GTFRow.h"
#include "Integrator.h"
#include <chrono>
#include <thread>
int main(int argc, char* argv[]) {
	auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::string> chromosomes;
    // default values (if not provided by user)
    int position_tolerance = 5;
    int min_length = 10;
    double coverage_tolerance = 0.8;
    double min_coverage = 0.05;
    std::string directory;
    int cores = 10;
    bool stranded = false;

    std::cout
    << "\n "

    "_____ _    ____ _____ ____  _____ ____  \n"
    "|  ___/ \\  / ___|_   _|  _ \\| ____|  _ \\ \n"
    "| |_ / _ \\ \\___ \\ | | | | | |  _| | |_) |\n"
    "|  _/ ___ \\ ___) || | | |_| | |___|  _ < \n"
    "|_|/_/   \\_\\____/ |_| |____/|_____|_| \\_\\"
    << "\n\nCopyright (c) 2025 Martina Lavanya Lehmann\n"
    << std::endl;
    std::cout << "Usage: fastder [options]\n\n"
                << "Options:\n"
                << "  --dir <path> ...             [REQUIRED] Path to directory with input files (relative path to build directory or absolute path). \n"
                << "                               Example: --dir ../../data/test_exon_skipping \n\n"
                << "  --chr <chr1> <chr2> ...      List of chromosomes to process. Default = ALL.\n"
                << "                               Example: --chr chr1 chr2 or --chr 1 2 \n\n"
                << "  --min-coverage <float>       Coverage threshold to qualify as an expressed region (ER), in [CPM]. \n"
                   "                               Normalization is done in-place by library size. \n"
                   "                               Default = 0.05 CPM.\n"
                << "                               Example: --min-coverage 0.25\n\n"
                << "  --position-tolerance <int>   Maximum permitted position deviation of splice junction and ER coordinates, in [nt]. Default = 5 nt.\n"
                << "                               Example: --position-tolerance 5\n\n"
                << "  --coverage-tolerance <float> Permitted coverage deviation within a stitched ER, as a proportion (e.g. 0.1 = 10 %).\n"
                << "                               The value is not strictly bound in (0,1). Default = 0.7.\n"
                << "                               Example: --coverage-tolerance 0.8\n\n"
                << "  --cores <int>                Number of cores that fastder may use. Default = 10 cores.\n"
                << "                               Example: --cores 23\n\n"
                << "  --stranded                   Enable stranded mode. Expects paired BedGraph files per sample\n"
                << "                               with 'plus'/'strand+' or 'minus'/'strand-' in filename.\n"
                << "                               Processes each strand independently and outputs strand-specific GTF.\n\n"
                << "Example:\n"
                << "  ./fastder --dir ../data --chr chr1 chr2 --position-tolerance 5 "
                 "--min-coverage 0.05 --coverage-tolerance 0.7 --cores 23\n"
                << std::endl;

    // parse command-line arguments
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "--chr") // --chr chr1 chr2
        {
            ++i;
            // the next argument starts with --
            while (i < argc && std::string(argv[i]).rfind("--", 0) == std::string::npos)
            {
                std::string chrom = argv[i];
                if (chrom.find("chr") == std::string::npos)
                {
                    chrom = "chr" + chrom; // add "chr" prefix if user passes --chr 1 2 3
                }
                chromosomes.push_back(chrom);
                ++i;
            }
            // decrement i again
            --i;

        }
        else if (arg == "--position-tolerance")
        {
            position_tolerance = atoi(argv[++i]);
        }
        else if (arg == "--cores")
        {
            cores = atoi(argv[++i]);
        }
        else if (arg == "--min-length")
        {
            min_length = atoi(argv[++i]);
        }
        else if (arg == "--min-coverage")
        {
            min_coverage = std::stod(argv[++i]);
        }

        else if (arg == "--coverage-tolerance")
        {
            coverage_tolerance = std::stod(argv[++i]);
        }

        else if (arg == "--dir")
        {
            directory = argv[++i];
        }
        else if (arg == "--stranded")
        {
            stranded = true;
        }
        else
        {
            std::cout << "[ERROR] Unknown argument '" << argv[i] << "'" << std::endl;
            return 1;
        }
    }

    // exit if no directory is provided
    if (directory.empty())
    {
        std::cerr << "[ERROR] No directory with input files specified! Quitting..." << std::endl;
        return 1;
    }

    // parse files
    std::cout << "[INFO] Expecting to parse MM, RR, BedGraph and Metadata CSV files from " << directory << std::endl;
    if (stranded) std::cout << "[INFO] Stranded mode enabled." << std::endl;
    Parser parser(directory, chromosomes, cores, stranded);
    parser.search_directory();

    // print parsing duration
    auto end_parsing = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parsing = end_parsing - start;
    std::cout << "[INFO] Parsing took " << elapsed_parsing.count() << " seconds.\n \n";

    // build output path
    std::string prefix = (directory.back() == '/') ? "FASTDER_RESULT_POS_TOL_" : "/FASTDER_RESULT_POS_TOL_";
    std::string output_path = directory + prefix + std::to_string(position_tolerance) + "_MIN_COV_" + std::to_string(min_coverage) + "_COV_TOL_" + std::to_string(coverage_tolerance) + "_MIN_LENGTH_" + std::to_string(min_length) + ".gtf";

    Integrator integrator(coverage_tolerance, position_tolerance);

    if (stranded)
    {
        // --- Plus strand ---
        std::cout << "\n[INFO] === Processing plus (+) strand ===" << std::endl;
        Averager avg_plus(cores);
        avg_plus.compute_mean_coverage(parser.all_per_base_coverages_plus);
        avg_plus.find_ERs(min_coverage, min_length);

        auto plus_sjs = Integrator::filter_sjs_by_strand(parser.mm_chrom_sj, parser.rr_all_sj, true);
        Integrator int_plus(coverage_tolerance, position_tolerance);
        int_plus.stitch_up(avg_plus.expressed_regions, plus_sjs, parser.rr_all_sj);
        for (auto& ser : int_plus.stitched_ERs) ser.strand = "+";

        // --- Minus strand ---
        std::cout << "\n[INFO] === Processing minus (-) strand ===" << std::endl;
        Averager avg_minus(cores);
        avg_minus.compute_mean_coverage(parser.all_per_base_coverages_minus);
        avg_minus.find_ERs(min_coverage, min_length);

        auto minus_sjs = Integrator::filter_sjs_by_strand(parser.mm_chrom_sj, parser.rr_all_sj, false);
        Integrator int_minus(coverage_tolerance, position_tolerance);
        int_minus.stitch_up(avg_minus.expressed_regions, minus_sjs, parser.rr_all_sj);
        for (auto& ser : int_minus.stitched_ERs) ser.strand = "-";

        // --- Merge results ---
        integrator.stitched_ERs.insert(integrator.stitched_ERs.end(),
            int_plus.stitched_ERs.begin(), int_plus.stitched_ERs.end());
        integrator.stitched_ERs.insert(integrator.stitched_ERs.end(),
            int_minus.stitched_ERs.begin(), int_minus.stitched_ERs.end());

        // Summary
        std::cout << "\n[INFO] Plus-strand Stitched ERs: " << int_plus.stitched_ERs.size() << std::endl;
        std::cout << "[INFO] Minus-strand Stitched ERs: " << int_minus.stitched_ERs.size() << std::endl;
        std::cout << "[INFO] Total Stitched ERs (both strands): " << integrator.stitched_ERs.size() << std::endl;

        for (auto& c : plus_sjs)
            std::cout << "[INFO] Plus-strand splice junctions in " << c.first << " : " << c.second.size() << std::endl;
        for (auto& c : minus_sjs)
            std::cout << "[INFO] Minus-strand splice junctions in " << c.first << " : " << c.second.size() << std::endl;
    }
    else
    {
        // --- Unstranded mode (original pipeline) ---
        Averager averager(cores);
        averager.compute_mean_coverage(parser.all_per_base_coverages);

        averager.find_ERs(min_coverage, min_length);

        integrator.stitch_up(averager.expressed_regions, parser.mm_chrom_sj, parser.rr_all_sj);

        // SUMMARY OF SPLICE JUNCTIONS
        for (auto& c : parser.mm_chrom_sj)
        {
            std::cout << "[INFO] Splice junctions in " << c.first << " : " << c.second.size() << std::endl;
        }

        // SUMMARY OF EXPRESSED REGIONS
        for (auto& c : averager.expressed_regions)
        {
            std::cout << "[INFO] Expressed Regions in " << c.first << ": " << c.second.size() << std::endl;
        }

        // SUMMARY OF STITCHED EXPRESSED REGIONS
        std::unordered_map<std::string, int> ser_counts;
        for (auto& stitched_er : integrator.stitched_ERs)
        {
            ser_counts[stitched_er.chrom]++;
        }

        std::cout << "[INFO] Total Stitched ERs across all chromosomes: " << integrator.stitched_ERs.size() << std::endl;
        for (auto& c : ser_counts)
        {
            std::cout << "[INFO] Stitched ERs in " << c.first << ": " << c.second << std::endl;
        }
    }

    // convert to GTF format
    integrator.write_to_gtf(output_path);

    // print duration
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "[INFO] Total duration (wallclock time): " << elapsed.count() << " seconds.\n";

    return 0;
}
