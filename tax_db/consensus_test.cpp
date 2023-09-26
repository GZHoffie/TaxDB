#include "indexer/consensus_indexer.h"
#include "./utils.h"

#include <map>
#include <chrono>




int main() {

    std::filesystem::path genome_path("/home/zhenhao/data/taxonomy/isolates");
    std::filesystem::path index_path("/home/zhenhao/TaxDB/tax_db/benchmark/index/");
    std::filesystem::path query_path("/home/zhenhao/data/taxonomy/mock/sampled.fastq");
    std::filesystem::path dict_path("/home/zhenhao/data/taxonomy/genome_id_genus_lookup.txt");
    std::string indicator("DB");

    consensus_indexer i;
    i.read_reference_from_file(genome_path / "2745495.fasta", 2745495);
    i.read_sequence_from_file(genome_path / "485870.fasta", 485870);


    
    

}