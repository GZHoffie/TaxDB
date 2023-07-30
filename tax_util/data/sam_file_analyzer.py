from tax_util.data.tax_dictionary import TaxonomyDict

class SamFileAnalyzer:
    def __init__(self, taxDict) -> None:
        self.groundTruth = []
        self.taxDict = taxDict

    def readGroundTruthFile(self, groundTruthFile):
        with open(groundTruthFile, 'r') as f:
            for line in f.readlines():
                if line.startswith('@'):
                    # Get the NCBI id of the species where each read is coming from
                    _, ref = line.split(' ')
                    info = ref.split('|')
                    ncbiID = self.taxDict.getID(info[-1].strip())
                    self.groundTruth.append(ncbiID)
    

    def analyzeSamFile(self, samFile):
        correct_genus = 0
        correct_species = 0
        with open(samFile, 'r') as f:
            for line in f.readlines():
                if not line.startswith('@'):
                    index = int(line.split(' ')[0])
                    ref = line.split('\t')[2]
                    info = ref.split('|')
                    ncbiID = self.taxDict.getID(info[-1].strip())
                    if ncbiID == self.groundTruth[index]:
                        correct_species += 1

                    predictedGenus = self.taxDict.getGenus(ncbiID)
                    correctGenus = self.taxDict.getGenus(self.groundTruth[index])
                    #print(predictedGenus, correctGenus)
                    if predictedGenus == correctGenus:
                        correct_genus += 1
                        
        print(correct_species, correct_genus, "out of", index + 1, "are correct")

if __name__ == "__main__":
    td = TaxonomyDict()
    td.readLookupTable("/home/zhenhao/data/taxonomy/genome_id_lookup.txt")
    sa = SamFileAnalyzer(td)
    sa.readGroundTruthFile("/home/zhenhao/data/taxonomy/test/DB_300_0.1_10K.fastq")
    sa.analyzeSamFile("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bucketmap_map_full.sam")
    sa.analyzeSamFile("/home/zhenhao/bucket-map/bucket_map/benchmark/output/bucketmap_map.sam")



