"""
Wrapper for doing SASA and deltaSASA calculations
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os
import sys
import re
import subprocess
import gzip
from collections import OrderedDict


def load_structure(path):
    """Read PDB file. 

    Args:
        path (str): location of PDB file.

    Returns:
        {str: list(str)}: Chain ID --> atom lines from file.

    """
    chains = {}
    regex = re.compile('^ATOM')
    with open(path) as fp:
        for line in fp:
            if regex.search(line):
                chain = line[21]
                chains.setdefault(chain, []).append(line.strip())
    return chains


def write_chains(structure, chainsToWrite, outPath):
    """Write out a PDB file.

    Args:
        structure ({str: list(str)}): Chain ID --> PDB file atom lines.
        chainsToWrite (str): chain IDs to write to file.
        outPath (str): location of output file.

    """
    with open(outPath, 'w') as fh:
        for chain in chainsToWrite:
            if chain not in structure:
                raise UserWarning('Requested chain: ' + chain + ' not loaded.')
            for line in structure[chain]:
                fh.write(line)
                fh.write("\n")

def run_dssp(path, chainIDs, dsspPath):
    """Run dssp on a PDB file and parse the output.

    Args:
        path (str): location of PDB file.
        chainIDs (str): chains to return data for.
        dsspPath (str): path to DSSP executable.

    Returns:
        {str: list(str)/list(int)}: single letter amino acid code
                                    and solvent accessible area of
                                    residue for each chain.
    """
    if not os.path.exists(path):
        raise ValueError('Input file missing for dssp: ' + path)
    strCommand = dsspPath + ' -i ' + path
    cmd = subprocess.run(strCommand, shell=True, stdout=subprocess.PIPE)
    data = {'aminos': {}, 'acc': {}, 'id': {}}
    started = False
    results = cmd.stdout.decode('utf-8')
    for line in results.splitlines():
        # Skip the header
        if not started:
            regex = re.compile('^  #  ')
            if regex.search(line):
                started = True
                continue
        else:
            chainID = line[11]
            if chainID in chainIDs:
                acc = int(line[34:38])
                data['acc'].setdefault(chainID, []).append(acc)
                amino = line[13]
                data['aminos'].setdefault(chainID, []).append(amino)
                residueID = int(line[6:10])
                data['id'].setdefault(chainID, []).append(residueID)
   
    return data


def rescale(values, sequence, maxSasaPerResidue):
    """Scale solvent accessable surface area to relative solvent accessability.

    Args:
        values: solvent accessibility of each protein.
        sequence (list(str)): amino acid sequence of protein.
        maxSasaPerResidue ({str: float}): single letter amino acid code
                                          --> maximum solvent accessibility
    """
    return [min(1, values[i] / maxSasaPerResidue[sequence[i]]) for i in range(len(values))]


def sasa_scan(structure, chainIDs, chainsInComplex, dsspPath, normalize,
              tmpPath='/tmp/atoms.pdb'): 
    """

    Args:
        structure: dict(str: list(str)):
                                Chain ID -> atom lines from parsed PDB file.
        chainIDs (str): chains of interest.
        chainsInComplex (str): chains in complex of interest.
        normalize (bool): normalize to get relative solvent accessability.
        tmpPath (str): path to write out temporary PDB structure.
        dsspPath (str): path to DSSP executable.

    Returns:
        list(OrderedDict(int: dict)):
            PDB residue ID to dssp results of amino acid, solvent accessability
            and difference in solvent accessability between the complexed and
            uncomplexed states
    """
    maxSasa = {'A': 115.,
               'R': 175.,
               'N': 154.,
               'D': 154.,
               'C': 118.,
               'Q': 180.,
               'E': 177.,
               'G': 91.,
               'H': 179.,
               'I': 147.,
               'L': 159.,
               'K': 206.,
               'M': 198.,
               'F': 188.,
               'P': 145.,
               'S': 128.,
               'T': 140.,
               'W': 181.,
               'Y': 178.,
               'V': 135.}

    # compute sasa for each chain in the lone state
    dsspOutput = {}
    for chainID in chainIDs:
        write_chains(structure, chainID, tmpPath)
        dsspOutput[chainID] = run_dssp(tmpPath, chainID, dsspPath)

    # compute sasa for each chain in the complex
    write_chains(structure, chainsInComplex, tmpPath)

    dsspOutput[chainsInComplex] = run_dssp(tmpPath, chainsInComplex, dsspPath)
   
    results = []
    for chainID in chainIDs:
        residueIDs = dsspOutput[chainID]['id'][chainID]
        # ERROR IS HERE, CHECK WHY THE CHAIN IS NOT VALID MAYBE DSSP IS NOT RUNNING PROPERLY WITH NEW INSTALL, MAYBE
        # OUTPUT FORMAT CHANGED AND NOT READING IT PROPERLY NOW
       
        seq = dsspOutput[chainID]['aminos'][chainID]
        sasaMonomer = dsspOutput[chainID]['acc'][chainID]
        sasaComplex = dsspOutput[chainsInComplex]['acc'][chainID]
        if len(sasaMonomer) != len(sasaComplex):
            raise UserWarning('Something went wrong in SASA code.')
        if normalize:
            sasaMonomer = rescale(sasaMonomer, seq, maxSasa)
            sasaComplex = rescale(sasaComplex, seq, maxSasa)
        delta = [chainRSA - complexRSA for chainRSA, complexRSA in
                 zip(sasaMonomer, sasaComplex)]
        resInfo = list(map(lambda a, b, c: {'aa': a, 'sasa': b, 'dsasa': c}, seq, sasaMonomer, delta))
        results.append(OrderedDict(zip(residueIDs, resInfo)))

    return results
