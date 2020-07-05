"""
Class to obtain amino acid properties

"""


class AminoAcid:
    # a class to obtain properties of an amino acid
    allAA = ['A', 'V', 'L', 'I', 'M', 'F', 'Y', 'W',
             'S', 'T', 'N', 'Q', 'G', 'P', 'C', 'H',
             'K', 'R', 'D', 'E']
    polarUnchargedAA = ['S', 'T', 'N', 'Q']
    specialAA = ['G', 'P', 'C']
    polarChargedAA = ['H', 'K', 'R', 'D', 'E']
    posChargedAA = ['H', 'K', 'R']
    negChargedAA = ['D', 'E']
    # using -ve/+ve hydropathy index from Kyte and Doolittle 1982
    hydrophobicAA = ['A', 'C', 'I', 'L', 'M', 'F', 'V']
    hydrophilicAA = ['D', 'E', 'G', 'H', 'K', 'N', 'P',
                     'Q', 'R', 'S', 'T', 'W', 'Y']

    def __init__(self, code):
        if code.upper() not in self.allAA:
            raise Exception('Invalid one letter amino acid code!')
        self.letter = code.upper()

    def is_hydrophobic(self):
        if self.letter in self.hydrophobicAA:
            return True
        else:
            return False

    def is_hydrophilic(self):
        if self.letter in self.hydrophilicAA:
            return True
        else:
            return False

    def is_polarUncharged(self):
        if self.letter in self.polarUnchargedAA:
            return True
        else:
            return False

    def is_special(self):
        if self.letter in self.specialAA:
            return True
        else:
            return False

    def is_polarCharged(self):
        if self.letter in self.polarChargedAA:
            return True
        else:
            return False

    def is_posCharged(self):
        if self.letter in self.posChargedAA:
            return True
        else:
            return False

    def is_negCharged(self):
        if self.letter in self.negChargedAA:
            return True
        else:
            return False
