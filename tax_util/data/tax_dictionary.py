from ete3 import NCBITaxa

class TaxonomyDict:
    """
    A dictionary that is able to
        - Translate names in the dataset to NCBI ids,
        - Get its information in the taxonomy tree.
    """
    def __init__(self) -> None:
        self.name2id = {}
        self.ncbi = NCBITaxa()
    
    def readLookupTable(self, fileName):
        """
        Read a txt file, where each line contains the name of the spiecies
        and its corresponding NCBI id.
        Args:
            - fileName (str): the name of the text file.
        """
        with open(fileName, 'r') as file:
            for line in file.readlines():
                name, ncbiID = line.split('\t')
                self.name2id[name] = int(ncbiID.strip())

    def getID(self, name):
        """
        Given the name of the species, return its NCBI id. If the name is unknown,
        return None.
        Args:
            - name (str): name of the spiecies in the dataset, which is the same as
              the text file that is previously read.
        """
        if name in self.name2id:
            return self.name2id[name]
        else:
            return None
    
    def getGenus(self, ID):
        """
        Find the genus of a given species ID.
        """
        if ID is None:
            return None
        
        # Get the lineage of the species
        lineage = self.ncbi.get_lineage(ID)
        ranks = self.ncbi.get_rank(lineage)
        #print(ranks)
        for i in ranks:
            if ranks[i] == 'genus':
                return i
        
        return None

        
    


if __name__ == "__main__":
    td = TaxonomyDict()
    td.readLookupTable("/home/zhenhao/data/taxonomy/genome_id_lookup.txt")
    print(td.getGenus('NC_014722.1'))
    #print(td.name2id)