FASTA_FILE=/home/zhenhao/data/taxonomy/DB.fa
BUCKET_LEN=65536
cd build
cmake ../tax_db/ -DBM_FASTA_FILE=${FASTA_FILE} -DBM_BUCKET_LEN=${BUCKET_LEN}