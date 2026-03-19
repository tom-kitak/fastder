//
// Created by marti on 08/10/2025.
//

#ifndef FASTDER_PARSE_H
#define FASTDER_PARSE_H
#include "SJRow.h"
#include <unordered_map>
#include <map>
#include <filesystem>
#include <fstream>
#include <unordered_set>

namespace fs = std::filesystem;

class Parser {
public:
    Parser(std::string path_, std::vector<std::string> chromosomes_, int cores_, bool stranded_ = false);
    void search_directory();

    // read in individual file types
    void read_all_bedgraphs(std::vector<std::string> bedgraph_files, unsigned int nof_threads,
                            std::vector<std::unordered_map<std::string, std::vector<double>>>& target_coverages);
    std::vector<BedGraphRow> read_bedgraph(const std::string& filename, uint64_t& library_size) const;
    void read_mm(std::string filename);
    void read_rr(std::string filename);
    void read_url_csv(std::string filename);
    void fill_up(std::vector<std::string> bedgraph_files);

    static void compute_per_base_coverage(const BedGraphRow& row, std::unordered_map<std::string, std::vector<double>>& per_base_coverage);

    // TODO add function get_rail_id_from_filename(filename)?
    bool stranded;
    unsigned int user_cores;
    std::string path;
    std::vector<std::string> chromosomes_vec; // for fast iteration
    std::unordered_set<std::string> chromosomes_set; // for fast check if chromosome is included
    std::vector<std::vector<BedGraphRow>> all_bedgraphs; //TODO maybe change to unordered map with key = sample id, value = bedgraph of the sample?
    std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages; //NOT ordered by chromosomes
    // stranded mode: separate coverage collections per strand
    std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages_plus;
    std::vector<std::unordered_map<std::string, std::vector<double>>> all_per_base_coverages_minus;
    // store RR info for each splice junctions
    std::vector<SJRow> rr_all_sj;

    // store relevant Market Matrix (MM)
    std::map<std::string, std::vector<uint64_t>> mm_chrom_sj; // <chrom, sj_id> ordered by sj_id, map of sj occurring in samples part of the user input

    std::vector<std::pair<unsigned int, std::string>> rail_id_to_ext_id; // <rail_id, external_id> for all samples in the dataset
    // later sorted by rail_id to receive rank (= mm_id)

    std::unordered_set<unsigned int> mm_ids; // unordered map for fast mm_id lookup


    const std::unordered_set<std::string> permitted_chromosomes =  {
        "chr1",
         "chr2",
         "chr3",
         "chr4",
         "chr5",
         "chr6",
         "chr7",
         "chr8",
        "chr9",
         "chr10",
         "chr11",
         "chr12",
         "chr13",
         "chr14",
         "chr15",
         "chr16",
         "chr17",
         "chr18",
         "chr19",
        "chr20",
         "chr21",
         "chr22",
         "chrX",
    };



};


#endif //FASTDER_PARSE_H