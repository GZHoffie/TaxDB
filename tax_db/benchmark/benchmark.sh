#!/bin/bash

GENOME_FILE="/home/zhenhao/data/taxonomy/DB.fa"
BENCHMARK_PATH="/home/zhenhao/bucket-map/bucket_map/benchmark"
QUERY_FILE="/home/zhenhao/data/taxonomy/test/DB_300_0.1_10K.fastq"
INDICATOR="DB"

# initialize benchmark directory if it doesnt exist
mkdir -p ${BENCHMARK_PATH}
mkdir -p "${BENCHMARK_PATH}/index"
mkdir -p "${BENCHMARK_PATH}/log"
mkdir -p "${BENCHMARK_PATH}/output"

# Run index benchmarking
#./benchmark_index.sh ${GENOME_FILE} ${BENCHMARK_PATH} ${INDICATOR}

# Run map benchmarking
./benchmark_map.sh ${QUERY_FILE} ${BENCHMARK_PATH} ${INDICATOR}
