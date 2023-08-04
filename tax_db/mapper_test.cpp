#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

void _output(std::filesystem::path output_file_path, std::vector<std::vector<std::string>> map_res) {
    // output results to the output file
    std::ofstream output(output_file_path);
    for (auto res : map_res) {
        for (auto name : res) {
            output << name << " ";
        }
        output << "\n";
    }
}


int main() {
    // parameters
    unsigned int bucket_len = 65536;
    unsigned int read_len = 300;
    seqan3::shape shape{0b111111111_shape};
    unsigned int samples = 30;
    unsigned int segment_samples = 20;
    unsigned int fault = 10; 
    float distinguishability = 0.6;

    std::filesystem::path genome_path("/home/zhenhao/data/taxonomy/DB.fa");
    std::filesystem::path index_path("/home/zhenhao/TaxDB/tax_db/benchmark/index/");
    std::filesystem::path query_path("/home/zhenhao/data/taxonomy/mock/com31.merged.fastq");
    std::string indicator("DB");

    // do the indexing 33147 19560
    auto indexer = new bucket_hash_indexer<33371>(bucket_len, read_len, shape, shape);
    indexer->index(genome_path, index_path, indicator);

    // perform the mapping
    auto mapper = new q_gram_mapper<33371>(bucket_len, read_len, segment_samples, shape, samples, fault, distinguishability);
    mapper->load(index_path, indicator);
    auto res = mapper->map(query_path);
    _output("/home/zhenhao/TaxDB/com31.output", res);
    
    
    delete indexer;
    delete mapper;



    
    

}