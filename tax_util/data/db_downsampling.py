from Bio import SeqIO
from tax_util.data.tax_dictionary import TaxonomyDict

class TaxDBDownSampler:
    def __init__(self, numSamples=1) -> None:
        self.numSamples = numSamples
        self.sampledSpecies = {}
    
    def sampleDatabase(self, fastaFile, taxDict, outputDir):
        self.numGenomes = 0
        sequences = []
        total_sequences = 0
        sampled_sequences = 0
        for record in SeqIO.parse(fastaFile, "fasta"):
            # Get the NCBI id of the species
            info = record.id.split('|')
            name = info[-1].strip()
            ncbiID = taxDict.getID(name)
            genus = taxDict.getGenus(name)
            total_sequences += 1

            # take the sequence if it belongs to a genus that has fewer than self.numSamples genomes
            if genus not in self.sampledSpecies:
                self.sampledSpecies[genus] = set([ncbiID])
            if ncbiID in self.sampledSpecies[genus] or len(self.sampledSpecies[genus]) < self.numSamples:
                self.sampledSpecies[genus].add(ncbiID)
                sequences.append(record)
                sampled_sequences += 1
            
        
        # Write the sampled genome to a new fasta file
        with open(outputDir, "w") as output_handle:
            SeqIO.write(sequences, output_handle, "fasta")

        print(sampled_sequences, "sequences sampled out of", total_sequences, "references. The output fasta is written to", outputDir)
            

if __name__ == "__main__":
    sampler = TaxDBDownSampler()
    td = TaxonomyDict()
    td.readLookupTable("/home/zhenhao/data/taxonomy/genome_id_lookup.txt")
    sampler.sampleDatabase("/home/zhenhao/data/taxonomy/DB.fa", td, "/home/zhenhao/data/taxonomy/DB_sampled.fa")