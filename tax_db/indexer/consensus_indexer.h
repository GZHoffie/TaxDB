#ifndef BUCKET_MAP_CONSENSUS_INDEXER_H
#define BUCKET_MAP_CONSENSUS_INDEXER_H

#include <fstream>
#include <bitset>
#include <string>
#include <iostream>

#include <seqan3/io/sequence_file/input.hpp>
#include "../utils.h"
#include "./indexer.h"
#include "./bucket_indexer.h"

/**
 * Just for experimenting purposes: far from optimized and generalizable.
 */

// Assume we use k=9 for k-mer length.
#define BIT_VECTOR_LENGTH 262144


class consensus_indexer {
private:

    std::vector<std::bitset<BIT_VECTOR_LENGTH>> kmers;
    std::vector<std::vector<std::tuple<unsigned int, unsigned int, bool>>> kmers_info;

    unsigned int ref_bucket_length = 65536, seq_bucket_length = 4396;
    unsigned int ref_overlap_length = 4096, seq_overlap_length = 300;

    void _read_reference(const seqan3::dna4_vector& sequence, unsigned int species_id) {
        unsigned int bucket_index = 0;

        auto operation = [&](const seqan3::bitpacked_sequence<seqan3::dna4>& seq, unsigned int id) {
            std::bitset<BIT_VECTOR_LENGTH> bucket_kmers;
            for (auto && value : seq | seqan3::views::kmer_hash(seqan3::ungapped{9})) {
                bucket_kmers.set(value);
            }

            // push them into the kmer storage
            kmers.push_back(bucket_kmers);
            std::vector<std::tuple<unsigned int, unsigned int, bool>> info_init;
            auto info = std::make_tuple(species_id, bucket_index, true);
            info_init.push_back(info);
            kmers_info.push_back(info_init);

            bucket_index++;
        };

        iterate_through_buckets(sequence, species_id, ref_bucket_length, ref_overlap_length, operation, true);
    }

    void _read_sequence(const seqan3::dna4_vector& sequence, unsigned int species_id) {
        unsigned int bucket_index = 0;

        auto operation = [&](const seqan3::bitpacked_sequence<seqan3::dna4>& seq, unsigned int id) {
            std::bitset<BIT_VECTOR_LENGTH> bucket_kmers_orig;
            for (auto && value : seq | seqan3::views::kmer_hash(seqan3::ungapped{9})) {
                bucket_kmers_orig.set(value);
            }

            std::bitset<BIT_VECTOR_LENGTH> bucket_kmers_revcomp;
            for (auto && value : seq | seqan3::views::complement | std::views::reverse | seqan3::views::kmer_hash(seqan3::ungapped{9})) {
                bucket_kmers_revcomp.set(value);
            }

            // for each kmer bucket in the database, find if there is a bucket that covers most of the k-mers in current bucket.
            unsigned int max_coverage = 0, best_bucket = 0;
            bool orig = true;
            for (unsigned int i = 0; i < kmers.size(); i++) {
                unsigned int coverage_orig = (bucket_kmers_orig & kmers[i]).count();
                unsigned int coverage_revcomp = (bucket_kmers_revcomp & kmers[i]).count();
                if (coverage_orig > max_coverage) {
                    max_coverage = coverage_orig; best_bucket = i; orig = true;
                }
                if (coverage_revcomp > max_coverage) {
                    max_coverage = coverage_revcomp; best_bucket = i; orig = false;
                }
            }

            seqan3::debug_stream << "Sub-bucket " << bucket_index << " is best covered by bucket " << best_bucket << ", coverage: " 
                                 << max_coverage << "; " << ((float) max_coverage) / bucket_kmers_orig.count() << ", original: " << orig << ".\n";
            // Check if the max coverage is large enough.
            if (((float) max_coverage) / bucket_kmers_orig.count() >= 0.8) {
                // if enough k-mers are covered, merge the small bucket into the large bucket.
                seqan3::debug_stream << "Merge into bucket " << best_bucket << "!\n";
            }

            bucket_index++;
        };

        iterate_through_buckets(sequence, species_id, seq_bucket_length, seq_overlap_length, operation, true);
    }


public:
    consensus_indexer() {}

    void read_reference_from_file(std::filesystem::path const & fasta_file_name, unsigned int species_id) {
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        for (auto && record : reference_genome) {
            seqan3::dna4_vector seq(record.sequence().begin(), record.sequence().end());
            _read_reference(seq, species_id);
        }
    }

    void read_sequence_from_file(std::filesystem::path const & fasta_file_name, unsigned int species_id) {
        seqan3::sequence_file_input reference_genome{fasta_file_name};
        for (auto && record : reference_genome) {
            seqan3::dna4_vector seq(record.sequence().begin(), record.sequence().end());
            _read_sequence(seq, species_id);
        }
    }

};

#endif