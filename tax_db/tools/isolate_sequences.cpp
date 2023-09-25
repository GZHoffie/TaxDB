
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>

class sequence_isolator {
private:
    std::unordered_map<std::string, unsigned int> name_to_species;
    std::unordered_map<std::string, unsigned int> name_to_genus;

    std::pair<unsigned int, unsigned int> _get_id(std::string sequence_name) {
        /**
         * @brief Hard coding for DB.fa. Extracts the sequence name from identifiers in the fasta file.
         * @returns a pair of integers, corresponding to the genus ID and species ID respectively.
         */
        std::string species_name(sequence_name.substr(sequence_name.rfind("|") + 1));
        return std::make_pair(name_to_genus[species_name], name_to_species[species_name]);
    }

public:
    void read_dictionary(std::filesystem::path dict_path) {
        /**
         * @brief Read the dictionary that maps sequence name to NCBI id.
         */
        std::ifstream is(dict_path);
        std::string name;
        unsigned int species_id, genus_id;
        while (is) {
            is >> name >> species_id >> genus_id;
            name_to_species[name] = species_id;
            name_to_genus[name] = genus_id;
        }
        is.close();
    }

    void isloate_sequences(const std::filesystem::path& fasta_file, const std::filesystem::path& output_dir, 
                           const std::vector<unsigned int>& species_of_interest) {
        
        using namespace seqan3::literals;
 
        using types = seqan3::type_list<std::vector<seqan3::dna5>, std::string>;
        using fields = seqan3::fields<seqan3::field::seq, seqan3::field::id>;
        using sequence_record_type = seqan3::sequence_record<types, fields>;


        seqan3::sequence_file_input reference_genome{fasta_file};
        unsigned int index = 0;
        for (auto && record : reference_genome) {
            auto ids = _get_id(record.id());
            auto species_id = std::get<1>(ids);
            if (!species_of_interest.empty() && std::find(species_of_interest.begin(), species_of_interest.end(), species_id) == species_of_interest.end()) {
                // not among the species of interest
                continue;
            }

            auto fasta_file = output_dir / (std::to_string(species_id) + ".fasta");
 
            // FASTA format detected, std::ofstream opened for file
            seqan3::sequence_file_output fout{fasta_file};
            sequence_record_type rec{std::move(record.sequence()), std::move(record.id())};
            fout.push_back(rec);
        }
    }
};

    



int main() {
    sequence_isolator s;
    s.read_dictionary("/home/zhenhao/data/taxonomy/genome_id_genus_lookup.txt");
    std::vector<unsigned int> species_of_interest{485870, 2745495};
    s.isloate_sequences("/home/zhenhao/data/taxonomy/DB.fa", "/home/zhenhao/data/taxonomy/isolates", species_of_interest);
}