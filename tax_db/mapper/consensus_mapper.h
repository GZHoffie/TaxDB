#include "../utils.h"
#include "./q_gram_mapper.h"

template<unsigned int NUM_KMERS>
class consensus_mapper : public mapper {
private:
    unsigned int genome_length_lower_bound;
    unsigned int k;
    seqan3::shape kmer_shape;

    // structures that stores the k-mers
    std::vector<std::vector<std::bitset<NUM_KMERS>>> kmers_storage;
    std::unordered_map<unsigned int, unsigned int> genus_id_to_bucket;
    
    // metadata
    std::unordered_map<std::string, unsigned int> name_to_species_id;
    std::unordered_map<std::string, unsigned int> name_to_genus_id;

    std::bitset<NUM_KMERS> _sequence_to_kmer_bits(const seqan3::dna4_vector & sequence) {
        std::bitset<NUM_KMERS> res;
        
    }


};