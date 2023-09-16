from tax_util.data.tax_dictionary import TaxonomyDict
import collections
import pandas as pd

class SamFileAnalyzer:
    def __init__(self, taxDict) -> None:
        self.groundTruth = []
        self.taxDict = taxDict
        self.baseline = None

    def readGroundTruthFastq(self, groundTruthFile):
        with open(groundTruthFile, 'r') as f:
            for line in f.readlines():
                if line.startswith('@'):
                    # Get the NCBI id of the species where each read is coming from
                    _, ref = line.split(' ')
                    info = ref.split('|')
                    ncbiID = self.taxDict.getID(info[-1].strip())
                    self.groundTruth.append(ncbiID)
    
    def readGroundTruthFile(self, groundTruthFile):
        with open(groundTruthFile, 'r') as f:
            for line in f.readlines():
                ncbiID = int(line)
                self.groundTruth.append(ncbiID)
    
    def readBaseline(self, baselineResults):
        self.baseline = pd.read_csv(baselineResults)
    

    def analyzeOutputFile(self, samFile):
        # Sensitivity
        correct_species_prediction = collections.Counter()
        correct_genus_prediction = collections.Counter()

        num_species = collections.Counter()
        num_genus = collections.Counter()

        # Precision
        num_correct_species_predictions = collections.Counter()
        num_correct_genus_predictions = collections.Counter()
        
        num_species_predictions = collections.Counter()
        num_genus_predictions = collections.Counter()

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
                        num_species_predictions[groundTruth] += 1
                        num_genus_predictions[groundTruthGenus] += 1
                    else:
                        continue
                    ncbiID = int(prediction)
                    if ncbiID == groundTruth:
                        correct_species = True
                        num_correct_species_predictions[groundTruth] += 1
                    
                    predictedGenus = self.taxDict.getGenus(ncbiID)
                    if predictedGenus == groundTruthGenus:
                        correct_genus = True
                        num_correct_genus_predictions[groundTruthGenus] += 1
                
                correct_species_prediction[groundTruth] += 1 if correct_species else 0
                correct_genus_prediction[groundTruthGenus] += 1 if correct_genus else 0

                num_species[groundTruth] += 1
                num_genus[groundTruthGenus] += 1


                
                index += 1

        print("Total number of reads:", index)

        for species in correct_species_prediction:
            correct_species_prediction[species] /= num_species[species]
            num_correct_species_predictions[species] /= num_species_predictions[species]
        
        for genus in correct_genus_prediction:
            correct_genus_prediction[genus] /= num_genus[genus]
            num_correct_genus_predictions[genus] /= num_genus_predictions[genus]
        

        print("Genus-level sensitivity:", correct_genus_prediction)
        print("Species-level sensitivity:", correct_species_prediction)
        print("Genus-level precision:", num_correct_genus_predictions)
        print("Species-level precision:", num_correct_species_predictions)

        sensitivities = []
        precisions = []
        f1s = []

        for species in self.baseline['species_id']:
            s = correct_species_prediction[species]
            p = num_correct_species_predictions[species]
            sensitivities.append(s)
            precisions.append(p)
            f1s.append(2 * s * p / (s + p))
        
        self.baseline['BM_Sensitivity'] = sensitivities
        self.baseline['BM_Precision'] = precisions
        self.baseline['BM_F1'] = f1s
        print(self.baseline)#[:, 'species_id', 'sensitivity', 'precision', 'f1', 'BM_Sensitivity', 'BM_Precision', 'BM_F1'])


        print()
        for genus in correct_genus_prediction:
            print(str(genus) + "\t\t" + str(correct_genus_prediction[genus]) + "\t\t\t" + str(num_correct_genus_predictions[genus]))

        print(num_species)
        

if __name__ == "__main__":
    td = TaxonomyDict()
    td.readLookupTable("/home/zhenhao/data/taxonomy/genome_id_lookup.txt")
    sa = SamFileAnalyzer(td)
    sa.readBaseline("/home/zhenhao/data/taxonomy/mock/kraken2_main_database_in_db_results_species.csv")
    sa.readGroundTruthFile("/home/zhenhao/data/taxonomy/mock/sampled_labels.fastq")
    sa.analyzeOutputFile("/home/zhenhao/TaxDB/sampled.output")



