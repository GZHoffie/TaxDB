#ifndef BUCKET_MAP_Q_GRAM_MAP_H
#define BUCKET_MAP_Q_GRAM_MAP_H


#include "../indexer/bucket_fm_indexer.h"
#include "../utils.h"
#include "./mapper.h"
#include "./quality_filter.h"

#include <string>
#include <vector>
#include <bitset>
#include <cmath>
#include <chrono>
#include <fstream>
#include <ranges>
#include <iostream>
#include <algorithm>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/views/minimiser.hpp>

using seqan3::operator""_shape;


template<unsigned int NUM_BUCKETS>
class bucket_filter {
    /**
     * @brief Novel data structure for fast bit-parallel pop-count of multiple BitSet objects.
     *        Eliminating those BitSets with more than `num_fault_tolerance` zeros.
     *        This can help filter out the correct buckets fast.
     */
private:
    unsigned int num_fault_tolerance;
    unsigned int allowed_max_candidate_buckets;
    std::vector<std::bitset<NUM_BUCKETS>> filters;
    
    // adaptive filter
    std::bitset<NUM_BUCKETS> eliminated;
    unsigned int num_matches;

    std::vector<unsigned int> _set_bits(unsigned int index) {
        /**
         * @brief Find all set bits in filters[index]
         * @param index indicates which filter we want to check.
         */
        std::vector<unsigned int> res;
        if (index >= num_fault_tolerance || index < 0) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The input index "
                                 << index << " exceeds the fault tolerance level." << '\n';
            return res;
        } else if (filters[index].none()) {
            return res;
        }
        for (int i = filters[index]._Find_first(); i < NUM_BUCKETS; i = filters[index]._Find_next(i)) {
            res.push_back(i);
        }
        return res;
    }

    unsigned int _index_in_circular_list(unsigned int index) {
        /**
         * @brief we record the filters in a circular list to avoid frequent copying during elimination of buckets.
         * @param index an unsigned int between 0 and `num_fault_tolerance - 1`.
         */
        return (index + num_matches) % num_fault_tolerance;
    }

    void _eliminate_worst_buckets() {
        /**
         * @brief If sufficient number of k-mers are missed, eliminate that bucket.
         */
        // record the eliminated buckets
        eliminated &= filters[_index_in_circular_list(0)];

        // empty the `num_matches`-th filter
        filters[_index_in_circular_list(0)].reset();
        num_matches += 1;

        //seqan3::debug_stream << "Eliminated worst bucket! Current num matches: " << num_matches << "\n";
    }

public:
    bucket_filter(unsigned int fault, unsigned int num_buckets) {
        num_fault_tolerance = fault;
        for (int i = 0; i < num_fault_tolerance; i++) {
            std::bitset<NUM_BUCKETS> filter;
            filters.push_back(filter);
        }
        eliminated.set();
        num_matches = 0;
    }

    void reset() {
        for (int i = 0; i < num_fault_tolerance; i++) {
            filters[i].reset();
        }
        eliminated.set();
        num_matches = 0;
    }

    bool read(std::bitset<NUM_BUCKETS>& input) {
        /**
         * @brief Read a set of bits and filter out the indices at which they are 0.
         * @param input a bitset with size NUM_BUCKETS indicating whether a certain k-mer appear
         *              in each bucket. (0 at position i indicates that the k-mer doesn't exist 
         *              in bucket i)
         * @returns a boolean value, true if a unique best bucket is found.
         */
        // ignore the eliminated buckets
        input &= eliminated;

        // Record how many time we see a 1 for each bucket
        for (unsigned int i = num_fault_tolerance - 1; i > 0; i--) {
            filters[_index_in_circular_list(i)] |= (filters[_index_in_circular_list(i-1)] & input);
        }
        // Modify filters[0]
        filters[_index_in_circular_list(0)] |= input;

        // check the number of best buckets observed
        unsigned int num_best_buckets = filters[_index_in_circular_list(num_fault_tolerance - 1)].count();
        if (num_best_buckets > 0) {
            // eliminate the worst buckets
            _eliminate_worst_buckets();
            if (num_best_buckets == 1) return true;
        }
        return false;
    }

    std::pair<std::vector<unsigned int>, unsigned int> best_results() {
        /**
         * @brief Return the buckets that contains the most number of k-mers,
         *        as well as the number of k-mers that are missed
         */
        std::vector<unsigned int> res;
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            res = _set_bits(_index_in_circular_list(i));
            if (!res.empty()) {
                return std::make_pair(res, num_matches + i + 1);
            }
        }
        return std::make_pair(res, 0);
    }


    int _check_bucket(unsigned int bucket_id) {
        /**
         * @brief Check how many errors actually present in the true bucket.
         */
        for (int i = num_fault_tolerance-1; i >= 0; i--) {
            if (filters[i][bucket_id]) {
                return i;
            }
        }
        return 0;
    }
};

template<unsigned int NUM_BUCKETS>
class distinguishability_filter {
    /**
     * @brief A filter that filter out Q-grams that appear in most of the buckets.
     */
private:
    unsigned int threshold;

public:
    std::vector<unsigned int> zeros;

    distinguishability_filter(float distinguishability) {
        /**
         * @brief Initializer of the filter.
         * @param distinguishability the percentage of zeros in the bitset for each Q-gram.
         */
        threshold = (unsigned int) (distinguishability * NUM_BUCKETS);
    }

    void read(const std::vector<std::bitset<NUM_BUCKETS>>& q_grams_index) {
        int valid_q_grams = 0;
        for (auto &i: q_grams_index) {
            unsigned int ones = i.count();
            if (NUM_BUCKETS - ones > threshold) {
                valid_q_grams++;
            }
            if (ones == 0) {
                // the q-gram doesn't appear at all
                zeros.push_back(NUM_BUCKETS);
            } else {
                zeros.push_back(NUM_BUCKETS - ones);
            }
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of Q-grams with distinguishability >= " << ((float) threshold) / NUM_BUCKETS << ": " 
                             << valid_q_grams << " (" << ((float) valid_q_grams) / q_grams_index.size() * 100 << "%).\n";
    }

    bool is_highly_distinguishable(unsigned int kmer_hash) {
        return zeros[kmer_hash] >= threshold;
    }
};






template<unsigned int NUM_BUCKETS>
class q_gram_mapper : public mapper {
private:
    // q-gram index related information
    seqan3::shape q_gram_shape;
    std::vector<std::bitset<NUM_BUCKETS>> q_grams_index;
    unsigned int q;
    int k; // length of k-mer used for query.
    unsigned int size;

    // bucket index related information
    std::vector<std::string> bucket_to_species;
    unsigned int bucket_length;
    unsigned int read_length;
    unsigned int num_segment_samples;

    // mapper related information
    unsigned int num_samples;
    unsigned int num_fault_tolerance;
    unsigned int allowed_max_candidate_buckets;

    // filter that filter out the most possible bucket
    bucket_filter<NUM_BUCKETS>* filter;

    // Q-gram filters for map efficiency
    distinguishability_filter<NUM_BUCKETS>* dist_filter;
    // quality filter for avoiding sequencing errors
    unsigned int min_base_quality;
    std::vector<bool> high_quality_kmers;

    // sampling k-mers
    Sampler* kmer_sampler;
    Sampler* segment_sampler;

    // Name/species/genus information
    std::unordered_map<std::string, int> name_to_species;
    std::unordered_map<std::string, int> name_to_genus;

    void _read_dictionary(std::filesystem::path dict_path) {
        std::ifstream is(dict_path);
        std::string name;
        unsigned int species_id, genus_id;
        while (is) {
            is >> name >> species_id >> genus_id;
            name_to_species[name] = species_id;
            name_to_genus[name] = genus_id;
        }
        seqan3::debug_stream << "[INFO]\t\tSuccessfully loaded dictionary from path " << dict_path << ".\n";
        is.close();
    }

    std::bitset<NUM_BUCKETS> _bitset_from_bytes(const std::vector<char>& buf) {
        /**
         * @brief Convert 8-byte chars to bitset.
         * Adopted from https://stackoverflow.com/a/7463972.
         */
        assert(buf.size() == ((NUM_BUCKETS + 7) >> 3));
        std::bitset<NUM_BUCKETS> result;
        for (unsigned int j = 0; j < NUM_BUCKETS; j++)
            result[j] = ((buf.at(j>>3) >> (j & 7)) & 1);
        return result;
    }

    /**
     * @brief Find all high quality k-mers.
     * 
     * @param quality the vector of quality alphabets.
     * fills in `high_quality_kmers`: A vector with the i-th element representing whether
     *                                the i-th k-mer is high-quality.
     * TODO: support spaced k-mer.
     */
    void _high_quality_kmers(const std::vector<seqan3::phred94>& quality) {
        unsigned int consecutive_high_quality_base = 0, index = 0;
        high_quality_kmers.clear();
        for (auto & qual: quality) {
            if (qual.to_rank() >= min_base_quality) {
                consecutive_high_quality_base++;
                if (consecutive_high_quality_base >= q) {
                    high_quality_kmers.push_back(true);
                    continue;
                }
            } else {
                consecutive_high_quality_base = 0;
            }
            index++;
            if (index >= q) {
                high_quality_kmers.push_back(false);
            }
        }
        assert(high_quality_kmers.size() == quality.size() - q + 1);
    }


    std::vector<int> _determine_candidate_species(const std::unordered_map<int, int>& counter_orig,
                                                  const std::unordered_map<int, int>& counter_rev_comp) {
        /**
         * @brief find the candidate species, given the outputs of the bucket mapping procedure.
         * 
         * @param counter_orig the candidate species given original string,
         * @param counter_rev_comp the candidate species given the reverse complement of the string
         * 
         * @return a vector of all species that have the highest number of votes
         * */                                                    

        // Find the highest number of votes
        int most_votes = 0;
        std::vector<int> res;
        for (auto &entry : counter_orig) {
            if (entry.second >= most_votes) {
                res.clear();
                res.push_back(entry.first);
                most_votes = entry.second;
            }
        }
        for (auto &entry : counter_rev_comp) {
            if (entry.second >= most_votes) {
                res.clear();
                res.push_back(entry.first);
                most_votes = entry.second;
            }
        }
        return res;
    
    }


public:
    q_gram_mapper(unsigned int bucket_len, unsigned int segment_len, unsigned int segment_samples, seqan3::shape shape, 
                  unsigned int num_bundled_kmers,
                  unsigned int samples, unsigned int fault, float distinguishability,
                  unsigned int quality_threshold = 35, 
                  unsigned int num_candidate_buckets = 30) : mapper() {
        // initialize private variables
        bucket_length = bucket_len;
        read_length = segment_len;
        num_segment_samples = segment_samples;
        
        q_gram_shape = shape;
        size = std::ranges::size(shape);
        q = shape.count();
        k = num_bundled_kmers;
        
        num_samples = samples;
        num_fault_tolerance = fault;
        allowed_max_candidate_buckets = num_candidate_buckets;

        // initialize filter
        filter = new bucket_filter<NUM_BUCKETS>(num_fault_tolerance, num_candidate_buckets);
        dist_filter = new distinguishability_filter<NUM_BUCKETS>(distinguishability);
        min_base_quality = quality_threshold * q;

        // Initialize sampler
        kmer_sampler = new Sampler(num_samples);
        segment_sampler = new Sampler(num_segment_samples);
    }

    ~q_gram_mapper() {
        delete filter;
        delete dist_filter;
        delete kmer_sampler;
        delete segment_sampler;
    }


    void load(std::filesystem::path const & dict_directory, std::filesystem::path const & index_directory, const std::string& indicator) {
        /**
         * @brief Look for index and pattern file inside the index_directory,
         *        read the files and store the values in class attribute.
         * @param dict_dicrectory the directory that stores the name to species/genus id file.
         * @param index_directory the directory that stores the q_gram_count_file.
         */
        // q_grams_index should be empty
        if (!q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is not empty. Terminating load.\n";
            return;
        }
        // initialize q_gram index
        int total_q_grams = (int) pow(4, q);
        int num_chars_per_q_gram = (NUM_BUCKETS + 7) >> 3;

        // Read the index file
        std::ifstream file(index_directory / (indicator + ".qgram"));
        if (file) {
            Timer clock;
            clock.tick();

            // read several bytes from file at a time
            std::vector<char> buffer(num_chars_per_q_gram);
            for (unsigned int i = 0; i < total_q_grams; i++) {
                file.read(&buffer[0], sizeof(unsigned char) * num_chars_per_q_gram);
                q_grams_index.push_back(_bitset_from_bytes(buffer));
            }

            dist_filter->read(q_grams_index);

            // Complete the read
            clock.tock();
            seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for loading index files: " 
                                 << clock.elapsed_seconds() << " s." << '\n';
            seqan3::debug_stream << "[INFO]\t\t" << "Successfully loaded " 
                                 << index_directory / (indicator + ".qgram") << "." << '\n';
        }

        // read the bucket_id file
        // NOTE: this is specific to the taxonomy dataset.
        std::ifstream bucket_id(index_directory / (indicator + ".bucket_id"));
        std::string name;
        while (std::getline(bucket_id, name)) {
            std::string species_name(name.substr(name.rfind("|") + 1));
            bucket_to_species.push_back(species_name);
        }

        // read the dictionary
        _read_dictionary(dict_directory);
    }


    std::pair<std::vector<unsigned int>, unsigned int> query(const std::ranges::input_range auto& q_gram_hash) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param q_gram_hash the vector containing all hash values of q-grams in the
         *                    query sequence.
         * @returns a vector of integers indicating the possible regions that the sequence
         *          may belong to, along with an integer showing how many k-mers are missed.
         */
        // q_grams_index should not be empty
        if (q_grams_index.empty()) {
            seqan3::debug_stream << "[ERROR]\t\t" << "The q-gram index is empty. Cannot accept query.\n";
            std::vector<unsigned int> res;
            return std::make_pair(res, 0);
        }
        // Reset the filter
        filter->reset();

        // randomly shuffle the order of k-mers to be queried
        std::vector<unsigned int> k_mer_index(q_gram_hash.size() - 2 * k);
        std::iota(k_mer_index.begin(), k_mer_index.end(), k);
        std::shuffle(k_mer_index.begin(), k_mer_index.end(), std::mt19937{std::random_device{}()});


        // insert those samples into the filter
        unsigned int num_kmers_queried = 0;
        for (auto index : k_mer_index) {
            std::bitset<NUM_BUCKETS> hashes;
            hashes.set();
            //seqan3::debug_stream << "Queried k-mer: ";
            for (int i = -k; i <= k; i++) {
                //seqan3::debug_stream << q_gram_hash[index + i] << " ";
                // use bitwise AND to find bucket that contains all of the k-mers
                hashes &= q_grams_index[q_gram_hash[index + i]];
            }
            //seqan3::debug_stream << hashes.count();
            //seqan3::debug_stream << "\n";
            num_kmers_queried++;
            //seqan3::debug_stream << "Queried kmer num: " << num_kmers_queried << "\n";
            if (filter->read(hashes) && num_kmers_queried >= 20) break;
        }
        //seqan3::debug_stream << filter->best_results() << " " << num_kmers_queried << "\n";
        if (((float) std::get<1>(filter->best_results())) / num_kmers_queried < 0.5) {
            std::vector<unsigned int> res;
            return std::make_pair(res, 0);
        }
        return filter->best_results();
        //return filter->all_results();
        //return filter->ok_results();
    }

    std::tuple<std::vector<unsigned int>, std::vector<unsigned int>, unsigned int> 
    query_sequence(const std::vector<seqan3::dna4>& sequence,
                   const std::vector<seqan3::phred94>& quality) {
        /**
         * @brief From `q_grams_index`, determine where the sequence may be coming from.
         * @param sequence A dna4 vector containing the query sequence.
         * @param quality the read quality of the sequence.
         * @returns two vectors of integers indicating the possible regions that the sequence
         *          may belong to. The first represent the buckets that the read can be mapped 
         *          directly, and the second is the ones where reads can be mapped after
         *          reverse complementing.
         */
        // initialize results
        std::vector<unsigned int> candidates_orig;
        std::vector<unsigned int> candidates_rev_comp;
        
        // get satisfactory k-mers that passes through quality filter and distinguishability filter
        auto kmers = sequence | seqan3::views::kmer_hash(q_gram_shape);
        int num_kmers = kmers.size();

        auto [buckets_orig, vote_orig] = query(kmers);

        // query the reverse complements of the sampled k-mers
        auto samples_rev_comp = sequence | std::views::reverse | seqan3::views::complement | seqan3::views::kmer_hash(q_gram_shape);
        //std::vector<std::vector<unsigned int>> samples_rev_comp_vec(samples_rev_comp.begin(), samples_rev_comp.end());
        auto [buckets_rev, vote_rev] = query(samples_rev_comp);

        // TODO: take the number of missed k-mers into consideration
        if (vote_rev >= vote_orig) candidates_rev_comp = buckets_rev;
        if (vote_orig >= vote_rev) candidates_orig = buckets_orig;
        
        return std::make_tuple(candidates_orig, candidates_rev_comp, std::max(vote_rev, vote_orig));
    }



    std::vector<std::vector<int>>
    map(std::filesystem::path const & sequence_file) {
        /**
         * @brief Read a query fastq file and output the ids of the sequence that are mapped 
         *        to each file.
         */
        // benchmarking percentage of mapped reads.
        unsigned int mapped_reads = 0;
        unsigned int num_buckets_orig = 0, num_buckets_rev_comp = 0;
        
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};

        // initialize returning result
        std::vector<std::vector<int>> res;

        Timer clock;
        clock.tick();
        for (auto & rec : fin) {
            // count the species that appears the most in the prediction.
            std::unordered_map<int, int> counter_orig;
            std::unordered_map<int, int> counter_rev_comp;
            std::unordered_set<int> votes;

            // initialize return value
            std::vector<int> predicted_species;

            // sample segments from the whole sequence
            // skip if the read is too short
            if (rec.sequence().size() <= read_length) {
                res.push_back(predicted_species);
                continue;
            }
            segment_sampler->sample_deterministically(rec.sequence().size() - read_length - 1);
            for (auto i : segment_sampler->samples) {
                seqan3::dna4_vector segment_sequence(rec.sequence().begin() + i, 
                                                     rec.sequence().begin() + i + read_length);
                std::vector<seqan3::phred94> segment_quality(rec.base_qualities().begin() + i, 
                                                             rec.base_qualities().begin() + i + read_length);
                auto [buckets_orig, buckets_rev_comp, num_votes] = query_sequence(segment_sequence, segment_quality);

                // count the votes, for original sequence and the reverse complement.
                votes.clear();
                for (auto bucket : buckets_orig) {
                    votes.emplace(name_to_species[bucket_to_species[bucket]]);
                }
                for (auto vote : votes) {
                    counter_orig[vote] += num_votes;
                }

                votes.clear();
                for (auto bucket : buckets_rev_comp) {
                    votes.emplace(name_to_species[bucket_to_species[bucket]]);
                }
                for (auto vote : votes) {
                    counter_rev_comp[vote] += num_votes;
                }

                if (!buckets_orig.empty() || !buckets_rev_comp.empty()) break;
            }

            predicted_species = _determine_candidate_species(counter_orig, counter_rev_comp);
            num_buckets_orig += predicted_species.size();
            res.push_back(predicted_species);


            if (!predicted_species.empty()) {
                ++mapped_reads;
            }
            ++num_records;
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket mapping: " 
                             << time << " s (" << time * 1000 * 1000 / num_records << " μs/seq).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of reads that have at least one candidate bucket: " 
                             << mapped_reads << "  (" << ((float)mapped_reads) / num_records * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets a read is mapped to: " 
                             << ((float) num_buckets_orig) / mapped_reads << ".\n";
        return res;
    }


    std::vector<std::vector<int>> _query_file(std::filesystem::path sequence_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Read a query fastq file and output the bucket ids each query belongs to.
         * TODO: include the quality information for fastq.
         */
        std::vector<std::vector<int>> res;
        seqan3::sequence_file_input<_phred94_traits> fin{sequence_file};
        Timer clock;
        clock.tick();
 
        for (auto & rec : fin) {
            auto [buckets_orig, buckets_rev_comp] = query_sequence(rec.sequence(), rec.base_qualities());
            res.push_back(buckets_orig);
        }
        clock.tock();
        float time = clock.elapsed_seconds();
        seqan3::debug_stream << "[BENCHMARK]\t" << "Elapsed time for bucket query: " 
                             << time << " s (" << time * 1000 * 1000 / res.size() << " μs/seq).\n";
        return res;
    }

    void _check_ground_truth(const std::vector<std::vector<int>>& query_results, std::filesystem::path ground_truth_file) {
        /**
         * * This function is just for benchmarking.
         * @brief Check the performance of query results against the ground truth.
         *        Output necessary information.
         * @note The ground truth file must be the one generated together with sequence
         *       fastq file.
         * TODO: also check the exact location.
         */
        std::ifstream is(ground_truth_file);
        int bucket, exact_location;
        std::string cigar;

        // Test statistics
        int correct_map = 0;
        int total_bucket_numbers = 0;
        std::map<int, int> bucket_number_map;

        for (int i = 0; i < query_results.size(); i++) {
            is >> bucket >> exact_location >> cigar;
            std::vector<int> buckets = query_results[i];

            if (std::find(buckets.begin(), buckets.end(), bucket) != buckets.end()) {
                correct_map++;
            }
            total_bucket_numbers += buckets.size();
            ++bucket_number_map[buckets.size()];
        }

        // output information
        seqan3::debug_stream << "[BENCHMARK]\t" << "Total number of sequences: " 
                             << query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Correct bucket predictions: " 
                             << correct_map << " (" << ((float) correct_map) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Average number of buckets returned: " 
                             << ((float) total_bucket_numbers) / query_results.size() << ".\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences that have no candidate bucket: " 
                             << bucket_number_map[0] << " (" << ((float) bucket_number_map[0]) / query_results.size() * 100 << "%).\n";
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of uniquely mapped sequences: " 
                             << bucket_number_map[1] << " (" << ((float) bucket_number_map[1]) / query_results.size() * 100 << "%).\n";
        int small_bucket_numbers = 0;
        for (int i = 1; i <= 5; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 5 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
        for (int i = 6; i <= 10; i++) {
            small_bucket_numbers += bucket_number_map[i];
        }
        seqan3::debug_stream << "[BENCHMARK]\t" << "Number of sequences mapped to <= 10 buckets: " 
                             << small_bucket_numbers << " (" << ((float) small_bucket_numbers) / query_results.size() * 100 << "%).\n";
    }

    void reset() {
        /**
         * @brief Release the memory, mainly used by the q_grams_index.
         * TODO: release other variables.
         */
        q_grams_index.clear();
        q_grams_index.shrink_to_fit();
    }
};


#endif