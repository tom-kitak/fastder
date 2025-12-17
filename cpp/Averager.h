//
// Created by marti on 08/10/2025.
//
#pragma once

#include "BedGraphRow.h"
#include <vector>
#include <iostream>
#include <unordered_map>
#include <mutex>
#ifndef FASTDER_AVERAGE_H
#define FASTDER_AVERAGE_H

#endif //FASTDER_AVERAGE_H

class Averager {

    public:
        Averager(int threads_);
        void compute_mean_coverage(std::vector<std::unordered_map<std::string, std::vector<double>>>& all_per_base_coverages);
        void find_ERs(double threshold, int min_length);

        // matrix consisting of multiple vectors, where each row in a vector is of type BedGraphRow,
        // a BedGraphRow is a bin of nucleotides with the same read count
        // a vector of BedGraphRows corresponds to one sample
        //matrix where each vector contains the normalized count per base from ONE sample
        int nof_threads;
        std::vector<std::string> chroms;
        std::mutex map_mutex;
        std::unordered_map<std::string, std::vector<double>> mean_coverage; //key = chromosome, value = vector of per-base mean coverage of this chromosome
        std::unordered_map<std::string, std::vector<BedGraphRow>> expressed_regions; //key = chromosome, value = vector of BedGraphRows with coverage > threshold
        // store all the individual sample maps in a vector (since sample identity doesn't matter anymore later on)



};