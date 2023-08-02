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
    unsigned int samples = 20;
    unsigned int segment_samples = 50;
    unsigned int fault = 20;
    float distinguishability = 0.5;

    std::filesystem::path genome_path("/home/zhenhao/data/taxonomy/DB_long_sampled.fa");
    std::filesystem::path index_path("/home/zhenhao/TaxDB/tax_db/benchmark/index/");
    std::filesystem::path query_path("/home/zhenhao/data/taxonomy/test/DB_long_12K_0.1_10K.fastq");
    std::string indicator("DB_long_sampled");

    // do the indexing 33147 19560
    auto indexer = new bucket_hash_indexer<19562>(bucket_len, read_len, shape, shape);
    indexer->index(genome_path, index_path, indicator);

    // perform the mapping
    auto mapper = new q_gram_mapper<19562>(bucket_len, read_len, segment_samples, shape, samples, fault, distinguishability);
    mapper->load(index_path, indicator);
    auto res = mapper->map(query_path);
    _output("/home/zhenhao/TaxDB/DB_12K_0.1_10K.output", res);
    
    
    delete indexer;
    delete mapper;



    
    

}