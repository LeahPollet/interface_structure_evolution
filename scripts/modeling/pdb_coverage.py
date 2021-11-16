import os
import time
from tqdm import tqdm
from Bio import PDB, SeqIO

"""
Step 2 in the pipeline - 
Generates a aa and residue length summary file for a list of PDBs 

Input:  pdb_reduced.faa (AA sequence and pdbID of all proteins in the analysis in fasta format)
        assemblies (folder with the curated, quality controled PDB files for the analysis)
Output: pdb_lengths.csv (Format: pdb,chain,struct_length,seq_length)

Output file used in build_protein_models.py

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

def load_sequences():
    fH = open('../../data/processed/pdb_reduced.faa', 'r')
    seqs = SeqIO.to_dict(SeqIO.parse(fH, 'fasta'))
    fH.close()
    lengths = {k: len(v) for k, v in seqs.items()}
    return lengths


# read in pdb file, go through each chain and return length of pdb seq
def read_pdb(fPath):
    parser = PDB.PDBParser(PERMISSIVE=1, QUIET=1)
    structure = parser.get_structure('name', (fPath))
    lengths = {}
    for chain in structure[0]:
        lengths[chain.id] = len(set([r.id[1] for r in chain.get_unpacked_list()]))
    return lengths


def main():
    structLengths = {}
    seqLengths = load_sequences()
    opath = '../../data/processed/pdb_lengths.csv'
    if not os.path.exists(opath): 
        fOut = open(opath, 'w')
        fOut.write('pdb,chain,struct_length,seq_length\n')
        pdbDir = '../../data/processed/assemblies'
        pdbs = [i[:4] for i in os.listdir(pdbDir)]
        t1 = time.perf_counter()
        for pdbID in tqdm(pdbs):
            structLengths[pdbID] = read_pdb(pdbDir+'/'+pdbID+'.pdb1')
        for pdbID in pdbs:
            for chainID in structLengths[pdbID]:
                if pdbID+'_'+chainID not in seqLengths:
                    continue
                # Imaged structure shouldn't be longer than whole protein
                # but this happens in cases where substrates are mislabelled 
                # as residues e.g. 5CEV, 1BGV
                # should implement a procedure to remove those resididues
                # and then remove this hack
                if structLengths[pdbID][chainID] > seqLengths[pdbID+'_'+chainID]:
                    structLengths[pdbID][chainID] = seqLengths[pdbID+'_'+chainID]
                fOut.write(','.join([pdbID, chainID, 
                                     str(structLengths[pdbID][chainID]),
                                     str(seqLengths[pdbID+'_'+chainID])])+'\n')
        print ('Took', time.perf_counter() - t1)
        fOut.close()
    else:
        print ('pdb_lengths.csv file already exists. Skipping.')


if __name__ == '__main__':
    main()
