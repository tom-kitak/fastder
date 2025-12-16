//
// Created by martinalavanya on 24.09.25.
//

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <unordered_map>
#include <cassert>
#include "BedGraphRow.h"
#include "SJRow.h"
#include "Averager.h"
#include <future>

Averager::Averager(int threads_)
{
    nof_threads = threads_;
}

//compute mean coverage vector across samples
void Averager::compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages)
{

    if (all_per_base_coverages.empty())
    {
        std::cerr << "[ERROR] No per-base coverages were computed...";
        return;
    }

    // reserve storage
    mean_coverage.reserve(all_per_base_coverages[0].size());

    std::vector<std::string> chroms;
    for (auto& item : all_per_base_coverages.at(0))
    {
        chroms.push_back(item.first);
    }

    // store threads
    std::vector<std::thread> threads;
    threads.reserve(nof_threads);
    std::atomic_int next_index{0};

    std::cout << "[INFO] Provided #samples = " << all_per_base_coverages.size() << ", #chromosomes = "<< all_per_base_coverages[0].size() << std::endl;

    // parallel iteration over chromosomes. pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (unsigned int t = 0; t < nof_threads; ++t)
    {
        // capture everything
        threads.emplace_back([&]()
        {
            while (true)
            {
                unsigned int idx = next_index++;
                if (idx >= chroms.size()) break; // leave loop after reaching the last chromosome
                const std::string chrom = chroms[idx];
                std::vector<double> chrom_mean_vector;

                chrom_mean_vector.reserve(all_per_base_coverages[0].at(chrom).size());

                // compute mean coverage across samples
                for (unsigned int i = 0; i < all_per_base_coverages[0][chrom].size(); i++)
                {
                    double sum = 0.0;
                    // iterate over all positions i across samples j
                    for (unsigned int j = 0; j < all_per_base_coverages.size(); j++)
                    {
                        // all_per_base_coverages[sample_nr][chromosome][base_pair_position]
                        sum += all_per_base_coverages[j][chrom][i]; //sample j, chromosome chrom, position i
                    }
                    chrom_mean_vector.push_back(sum / all_per_base_coverages.size()); // (sum over coverage at bp i) / (#nof samples)
                }
                mean_coverage[chrom] = std::move(chrom_mean_vector);
            }
        });
    }
    // join threads
        for (auto& thr: threads) {
            thr.join();
        }
}

// identify ERs with coverage > threshold and length > min_length bp
void Averager::find_ERs(double threshold, int min_length)
{
    // define workers --> use one thread per chromosome
    std::unordered_map<std::string, std::future<std::vector<BedGraphRow>>> workers;
    workers.reserve(mean_coverage.size()); //pre-allocate space for workers

    //if length > 5 and coverage at bp < 5 -> ER has ended, append
    // if coverage at bp > 5 --> add to current ER
    // if coverage < 5 && length < 5: don't count as ER, reset start and average
    //pair.first = chromosome, pair.second = vector of per base coverages of that chromosome
    for (auto& coverage : mean_coverage)
    {
        std::string chrom = coverage.first;
        workers[chrom] = std::async(std::launch::async, [&, chrom]
        {
            int start = 0;
            double current_sum = 0;
            int count = 0;
            std::vector<BedGraphRow> chrom_expressed_regions;
            for (unsigned int i = 0; i <= mean_coverage[chrom].size(); i++)
            {
                // append the last expressed region if it's long enough
                if (i == mean_coverage[chrom].size()) {
                    if ((i - start) > min_length)
                    {
                        //TODO last ER can be inclusive or exclusive depending on mean_coverage[chrom][i] <= threshold, but this is currently not encoded
                        double current_avg = current_sum / (i - start);
                        BedGraphRow expressed_region = BedGraphRow(chrom, start, i, current_avg);
                        chrom_expressed_regions.push_back(expressed_region);
                    }
                    break;
                }

                double coverage =  mean_coverage[chrom][i];
                // coverage is less than threshold
                if (mean_coverage[chrom][i] <= threshold)
                {
                    // region at least 5 bp long, append to results
                    if ((i - start) > min_length)
                    {
                        double current_avg = current_sum / (i - start);

                        // assert(count == i - start);
                        BedGraphRow expressed_region = BedGraphRow(chrom, start, i, current_avg);
                        chrom_expressed_regions.push_back(expressed_region);

                    }
                    //region too short to be appended, reset start and current avg
                    start = i + 1;
                    count = 0;
                    current_sum = 0;

                }
                // add to current ER
                else if (coverage > threshold)
                {
                    current_sum += coverage;
                    ++count;
                }

            } // end inner for loop
                return chrom_expressed_regions;
            });
        }



        // single-threaded merge
        // TODO later use >1 thread per chromosome
        for (auto& pair : mean_coverage)
        {
            std::string chrom = pair.first;
            expressed_regions[chrom] = workers[chrom].get(); //get result
        }

}