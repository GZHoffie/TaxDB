#include "indexer/bucket_hash_indexer.h"
#include "mapper/q_gram_mapper.h"
#include "tools/short_read_simulator.h"
#include "./utils.h"

#include <map>
#include <chrono>

int main() {
    // parameters
    unsigned int bucket_len = 65536;
    unsigned int read_len = 310;
    seqan3::shape shape{0b111111111_shape};
    unsigned int samples = 30;
    unsigned int fault = 20;
    float distinguishability = 0.5;

    // perform the mapping
    auto mapper = new q_gram_mapper<33376>(bucket_len, read_len, shape, samples, fault, distinguishability);
    mapper->load("/home/zhenhao/bucket-map/bucket_map/benchmark/index/", "DB_bucketmap");
    auto res = mapper->map("/home/zhenhao/data/taxonomy/test/DB_300_0.1_10K.fastq");
    
    
    
    
    delete mapper;



    
    

}