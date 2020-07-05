"""
Models of yeast proteins that contain structural and
evolutionary residue-level annotations.

"""
import numpy as np
from Bio import Seq
from amino_acid import AminoAcid


def is_non_synonymous(cdn1, cdn2):
    """Do two codons code for different amino acids."""
    assert len(cdn1) == 3 and len(cdn2) == 3, 'ERROR with codons'
    if Seq.translate(cdn1) == Seq.translate(cdn2):
        return False
    else:
        return True


class Residue:
    def __init__(self, cdns, sasa, dsasa, aaLetter):
        self.id = ''                # Identifier notation for the residue in the PDB file (to match with dssp output: format examples MET1:A (AaPosition:chain))
        self.indexInSeq = np.nan    # Index of the residue with respect to the sequence of the entire protein as saved in the model (0:prot length)
        self.tags = ['AllResidues']
        self.rsa = sasa             # Relative solvent accessibility in uncomplexed state
        self.dRsa = dsasa           # Delta RSA
        self.aa = AminoAcid(aaLetter)
        if len(cdns) != 9:
            raise ValueError('Wrong number of codons. Expected 9')
        if not all([len(c) == 3 for c in cdns.values()]):
            raise ValueError('codons need to be 3 letters long')
        self.codons = cdns                              # Dict of yeast strain to 3-letter string
        self.contacts = np.nan                          # Residue contacts with interaction partners. Should initialize to np.nan
        self.distance_from_center = np.nan              # How close to center of interface
        self.distance_from_edges = np.nan               # How close to the most exposed surface of the interface
        self.substitution = np.nan                      # Evolutionry state of a residue
        self.consurf_score = np.nan                     # ConSurf score of a residue
        self.conf_int_consurf_score = (np.nan,np.nan)   
        
        
    def add_tag(self, tag):
        self.tags.append(tag)

    def n_nonsynonymous(self):
        """Total number of non-synonymous substitutions.

        Number of non-synonymous substitutions between pairs of most closely
        related species' codons in this modeled residue. If there are codons
        from 4 different species then this will return an int from 0-3

        Returns:
            int: number of non-synonymous substitutions.
        """
        return sum([is_non_synonymous(self.codons[i], self.codons[i + 1])
                    for i in range(len(self.codons) - 1)])


class Protein:
    def __init__(self, orfName, pdbID, chainID, aaSeq, resolution, cai=np.nan):
        
        """Model of a s. cer. protein with structural/evolutionary annotations.
        """
        self.residues = []
        self.interfaces = {}
        self.tags = ['AllProteins']
        self.orf = orfName
        self.resolution = resolution
        self.seq = aaSeq  # original aa sequence of s. cer. protein
        self.orthologs = {}
        self.pdb = (pdbID, chainID)   # PDB ID of aligned structure
        self.pdbComplex = ''          # chains of modeled complex
        # aa seq of just the structurally imaged residues of aligned PDB chain
        self.pdbStructureSeq = ''
        self.pdbWholeProteinSeq = ''  # aa seq of entire protein of PDB chain
        self.cai = cai
        self.goTerms = set()

    def add_tag(self, tag):
        self.tags.append(tag)

    def add_interface(self, partnerOrf, partnerChain):
        self.interfaces[partnerChain] = Interface(partnerOrf,
                                                  partnerChain,
                                                  self)

    def add_res(self, res, interfaceChain=None):
        self.residues.append(res)
        if interfaceChain is not None:
            if interfaceChain in self.interfaces:
                self.interfaces[interfaceChain].residues.append(res)
            else:
                raise UserWarning(interfaceChain + ' not in ' +
                                  ','.join(self.interfaces.keys()) +
                                  ' for protein ' + self.orf)

    def n_interfaces(self):
        return len(self.interfaces)

    def mean_rsa(self):
        return np.mean([r.rsa for r in self.residues])

    def coverage(self):
        """Fraction of original protein sequence covered.

        Coverage is lost from gaps in alignment with PDB
        structure and from gaps in alignment with otholgous
        yeast proteins.
        """
        return float(len(self.residues)) / float(len(self.seq))

    def pdb_coverage(self):
        """Fraction of the PDB protein in the PDB structure."""
        return (float(len(self.pdbStructureSeq)) /
                float(len(self.pdbWholeProteinSeq)))

    def __str__(self):
        return ('modeled s. cer. protein ' + self.orf +
                '\naligned with PDB structure ' + '-'.join(self.pdb) +
                '\naligned with yeast orthologs: ' +
                ', '.join(self.orthologs.values()) +
                '\nlength of protein ' + str(len(self.seq)) +
                '\nfraction of s. cer. protein covered by structure ' +
                'and alignment with orthologs %.2f' % self.coverage() +
                '\nlength of PDB structure ' + str(len(self.pdbStructureSeq)) +
                '\nlength of whole PDB protein sequence ' +
                str(len(self.pdbWholeProteinSeq)) +
                '\nresolution of PDB structure ' + str(self.resolution) +
                '\nCAI: ' + str(self.cai) +
                '\ncontains ' + str(len(self.residues)) + ' residues' +
                '\nwith ' + str(len(self.interfaces)) + ' PPI interfaces' +
                '\nGO annotations: ' + str(self.goTerms) +
                '\ntags: ' + str(self.tags))


class Interface:
    def __init__(self, partnerOrfID, chainID, protein):
        self.residues = []
        self.tags = []
        self.partnerProtein = partnerOrfID
        self.partnerChainID = chainID
        self.coexpression = -99.
        self.rankCoexpression = -99.
        self.protein = protein

    def n_residues(self):
        """Number of modeled residues in modeled interface."""
        return len(self.residues)
