#include "tools/short_read_simulator.h"


int main() {
    std::filesystem::path data_path = "/home/zhenhao/data/taxonomy";
    std::filesystem::path genome_file = data_path / "DB_long.fa";

    int bucket_length = 65536;
    int read_length = 12000;

    short_read_simulator sim(bucket_length, read_length, 0.01, 0.001, 0.001);
    sim.read(genome_file);
    sim.generate_fastq_file(data_path / "test", "DB_long_12K_0.1_10K", 10000);
}