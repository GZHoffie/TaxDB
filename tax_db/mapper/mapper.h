#ifndef BUCKET_MAP_MAPPER_H
#define BUCKET_MAP_MAPPER_H

class mapper {
    /**
     * @brief Virtual class for all mapper (specific to bucketMap), which map short reads
     *        to their corresponding bucket.
     */
public:
    mapper() {}

    unsigned int num_records = 0;
    
    /**
     * @brief Read a query fastq file and output the ids of the sequence that are mapped 
     *        to each bucket.
     * 
     * @returns A vector, where each element correspond to a read. Each element is a vector of integers,
     *          representing the predicted genus it belongs to.
     */
    virtual std::vector<std::vector<int>>
    map(std::filesystem::path const & sequence_file) = 0;

    /**
     * @brief Release the memory storing the sequences and index.
     */
    virtual void reset() = 0;
};

#endif