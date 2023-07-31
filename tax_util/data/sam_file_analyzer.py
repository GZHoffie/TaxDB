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
    

    def analyzeOutputFile(self, samFile):
        # Sensitivity
        correct_species_prediction = 0
        correct_genus_prediction = 0

        # Precision
        num_correct_species_predictions = 0
        num_correct_genus_predictions = 0
        num_predictions = 0

        index = 0
        with open(samFile, 'r') as f:
            for line in f.readlines():
                correct_species = False
                correct_genus = False
                groundTruth = self.groundTruth[index]
                #print(self.taxDict.getID(self.groundTruth[index]))
                groundTruthGenus = self.taxDict.getGenus(groundTruth)
                for prediction in line.split(' '):
                    if prediction != '\n':
                        num_predictions += 1
                    ncbiID = self.taxDict.getID(prediction)
                    if ncbiID == groundTruth:
                        correct_species = True
                        num_correct_species_predictions += 1
                    
                    predictedGenus = self.taxDict.getGenus(ncbiID)
                    if predictedGenus == groundTruthGenus:
                        correct_genus = True
                        num_correct_genus_predictions += 1
                
                correct_species_prediction += 1 if correct_species else 0
                correct_genus_prediction += 1 if correct_genus else 0


                
                index += 1

        print("Total number of predictions:", num_predictions)
        print("Genus-level sensitivity:", correct_genus_prediction / index)
        print("Species-level sensitivity:", correct_species_prediction / index)
        print("Genus-level precision:", num_correct_genus_predictions / num_predictions)
        print("Species-level precision:", num_correct_species_predictions / num_predictions)
        

if __name__ == "__main__":
    td = TaxonomyDict()
    td.readLookupTable("/home/zhenhao/data/taxonomy/genome_id_lookup.txt")
    sa = SamFileAnalyzer(td)
    sa.readGroundTruthFile("/home/zhenhao/data/taxonomy/test/DB_12K_0.1_10K.fastq")
    sa.analyzeOutputFile("/home/zhenhao/TaxDB/DB_12K_0.1_10K.output")



