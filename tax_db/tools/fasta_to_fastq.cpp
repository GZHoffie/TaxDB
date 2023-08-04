
#include <seqan3/io/sequence_file/all.hpp>


int main() {
    std::filesystem::path fasta_file = "/home/zhenhao/data/taxonomy/mock/_com1.merged.fasta";
    std::filesystem::path output_file = "/home/zhenhao/data/taxonomy/mock/_com1.merged.fastq";
    
    seqan3::sequence_file_input reference_genome{fasta_file};
    std::ofstream output{output_file};
    unsigned int index = 0;
    for (auto && record : reference_genome) {
        output << "@" << index << " " << record.id() << "\n";
        for (auto nt : record.sequence()) {
            output << nt.to_char();
        }
        // insert a quality string.
        output << "\n+\n" << std::string(record.sequence().size(), 'E') << "\n";
        index++;
    }
}