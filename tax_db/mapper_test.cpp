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
    unsigned int segment_samples = 30;
    unsigned int fault = 10;
    float distinguishability = 0.8;

    // perform the mapping
    auto mapper = new q_gram_mapper<19762>(bucket_len, read_len, segment_samples, shape, samples, fault, distinguishability);
    mapper->load("/home/zhenhao/bucket-map/bucket_map/benchmark/index/", "DB_part_bucketmap");
    auto res = mapper->map("/home/zhenhao/data/taxonomy/test/DB_12K_0.1_10K.fastq");
    _output("/home/zhenhao/TaxDB/DB_12K_0.1_10K.output", res);
    
    
    
    delete mapper;



    
    

}