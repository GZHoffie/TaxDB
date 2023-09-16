
#include <seqan3/io/sequence_file/all.hpp>

std::string _filename_abbreviate(std::string name) {
    std::string::size_type pos = name.find('_');
    if (pos != std::string::npos) {
        return name.substr(0, pos);
    } else {
        return name;
    }
}


int fasta_to_fastq(const std::filesystem::path& fasta_file, std::ofstream& output, unsigned int limit = 0) {
    seqan3::sequence_file_input reference_genome{fasta_file};
    unsigned int index = 0;
    for (auto && record : reference_genome) {
        output << "@" << index << " " << record.id() << "\n";
        for (auto nt : record.sequence()) {
            output << nt.to_char();
        }
        // insert a quality string.
        output << "\n+\n" << std::string(record.sequence().size(), 'E') << "\n";
        index++;
        if (index >= limit && limit != 0) break;
    }
    return index;
}


void sample_from_directory(const std::filesystem::path& fasta_dir, const std::filesystem::path& metadata_file, 
                           const std::filesystem::path& output_fastq_file, const std::filesystem::path& output_ground_truth_file, 
                           unsigned int limit = 0) {
    // read the metadata file
    std::unordered_map<std::string, std::string> filename_to_id;
    std::ifstream metadata(metadata_file);
    std::string line, filename, id;
    unsigned int line_index = 0, col_index = 0;
    while(getline(metadata, line)) {
		std::stringstream str(line);
        line_index++;
        if (line_index == 1) continue;

		
        getline(str, filename, ',');
        getline(str, id, ',');
        filename_to_id.emplace(_filename_abbreviate(filename), id);
    }

    // Sample from each fasta file
    std::ofstream output{output_fastq_file};
    std::ofstream gt{output_ground_truth_file};
    for (const auto & entry : std::filesystem::directory_iterator(fasta_dir)) {
        std::cout << entry.path().filename() << std::endl;
        int num_data = fasta_to_fastq(entry.path(), output, limit);

        // write the ground truth to the ground truth file
        for (int i = 0; i < num_data; i++) {
            gt << filename_to_id[_filename_abbreviate(entry.path().filename())] << "\n";
        }
    }
}


int main() {
    std::filesystem::path fasta_dir = "/home/zhenhao/data/taxonomy/mock/refs_in/refs_in/";
    std::filesystem::path gt_dir = "/home/zhenhao/data/taxonomy/mock/test_dataset_metadata.csv";
    std::filesystem::path output_file = "/home/zhenhao/data/taxonomy/mock/sampled.fastq";
    std::filesystem::path output_ground_truth = "/home/zhenhao/data/taxonomy/mock/sampled_labels.fastq";

    
    unsigned int limit = 1000;
    sample_from_directory(fasta_dir, gt_dir, output_file, output_ground_truth, limit);
}