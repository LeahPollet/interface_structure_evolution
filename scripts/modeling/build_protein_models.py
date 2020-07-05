#!/usr/bin/env python

""" Align RSA/dRSA from structural models to yeast protein sequences.

ToDo:
	
Output files:
	modeled_proteins.pkl - structurally modeled proteins

Input files: 
	Scer.aa / Spar.nt etc. - sequences of yeast ORFs
	scer_protein_interactions.txt - BioGrid physical s. cer. interactions
	blast_scer_vs_pdb.txt - blast results of s. cer. proteins against the PDB
	pdb_lengths.csv - length of each PDB structure
	pdb_reduced.faa - AA sequences of PDB chains
	curated_scer_pdb_homologs.txt - map of s. cer. proteins to homologous PDB chains
	pdb_chain_uniprot.csv - mapping of PDB chains to UniProt proteins

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os
import sys
import subprocess
import pickle
from collections import defaultdict, OrderedDict,Counter
import numpy as np
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq 
from Bio.PDB import PDBParser
from pattern_align_local import align
from sasa_scan import load_structure, sasa_scan
from modeled_proteins import Protein, Residue
import tarfile

def download_dssp():
	"""
	Downloads and compiles dssp for macosx

	Returns: 
		path_to_dssp_exec (str): path to the installed dssp executable

	Useful if dssp is not compiling properly, or to install on other OS:
	http://melolab.org/supmat/moma/install.html

	"""
	oPath = "../../data"
	url = "ftp://ftp.cmbi.ru.nl/pub/molbio/software/dssp-2/dssp-2.0.4.tgz"

	zipped = oPath+"/"+os.path.basename(url)
	path_to_dssp = oPath+"/"+os.path.basename(url)[:10]
	path_to_dssp_exec = oPath+"/"+os.path.basename(url)[:10]+"/mkdssp"
	if not os.path.exists(path_to_dssp):
		print("Downloading DSSP")
		os.system('wget -P '+ oPath +' '+url)
		tar = tarfile.open(zipped)
		tar.extractall(oPath)
		tar.close()
		os.remove(zipped)

		subprocess.Popen(["make"], stdout=subprocess.PIPE, cwd=path_to_dssp)

	return path_to_dssp_exec

def load_inputs(externalData,processedData):
	"""
	Load all of the data.
		Args: 
			externalData:  directory with the unprocessed data from previous analysis
			processedData: directory with the processed data from previous analysis
		
		Outputs:
			ints: 			Dictionary of S.cer PPIs (double sided)
			blastResults: 	Stats from the BLAST homology mapping process (between S.cer sequences and PDB sequences)
			pdbFullSeqs:	Dictionary with aa sequences for proteins in PDB
			aaSeq:			Dictionary with aa sequences for proteins in S.cere
			structs:		Dictionary with the mapping between ORF names and PDB structure and chain names
			pdbLengths:		Dictionary with the length of the structure and sequence for each chain
			resolutions:    Dictionary with the resolution of each PDB structure
			orthologs:		Dictionary with the mapping between ORF name in S. cere and ORF name in 3 orthologue species (S. bay, S.mik, S.par)
	"""
	# ints: Double sided dictionary of all S.cer PPIs
	intsList = [[l.split()[0], l.split()[1]] for l in open(os.path.join(processedData,'scer_protein_interactions_multiple_reports.txt'), 'r')]
	ints = defaultdict(set) 
	for y1, y2 in intsList:
		ints[y1].add(y2)
		ints[y2].add(y1)
	
	# blastResults: Stats from homology mapping
	blastResults = read_blast_results(os.path.join(processedData,'blast_scer_vs_pdb.txt'))

	# pdbFullSeqs: Sequences in FASTA format
	pdbFullSeqs = load_whole_protein_pdb_sequences(processedData)
	
	# aaSeq:	dictionary: { yeastORFname: aa sequence for S.cer}
	aaSeq, structs,  = {}, {}
	with open(os.path.join(externalData, 'Scer.aa'), 'r') as fH:
		aaSeq = {(k.split('_')[0]): str(v.seq) for k, v in
				SeqIO.to_dict(SeqIO.parse(fH, 'fasta')).items()}

	# structs: 	dictionary: { yeastORFname : [(struct(e.g: 5fpn), chain(e.g: A)), ...],}
	f_structs = open(os.path.join(processedData,'curated_scer_pdb_homologs.txt'), 'r')
	# count,tot = 0,0
	for l in f_structs:
		# count = count+1
		yeastORFname = l.split()[0][:-5] # strip the trailing "_mRNA"
		struct = l.split()[1]
		chain = l.split()[2]
		if yeastORFname not in structs:
			# tot = tot+1
			structs[yeastORFname] = []
		structs[yeastORFname].append((struct, chain))
	f_structs.close()

	# pdbLengths: dictionary {(pdbID, chainID)] : (length of strcuture, length of sequence)}
	with open(os.path.join(processedData,'pdb_lengths.csv'), 'r') as fCov:
		pdbLengths = {}
		for l in fCov.readlines()[1:]:
			l = l.strip('\n') 	# /!\ NT strip the trailing "\n "at the end of lSeq
			pdbID, chainID, lStruct, lSeq = l.split(',')
			pdbLengths[(pdbID, chainID)] = (lStruct, lSeq)

	# resolutions: resolution of each PDB structure.
	resolutions = {}
	with open(os.path.join(externalData,'resolu.idx'), 'r') as fResolution:
		for line in fResolution.readlines()[6:]:
			pdbID = line.split(';')[0].strip().lower()
			data_resolution = line.split(';')[1].strip()
			if data_resolution == '':
				continue
			res = float(data_resolution)
			if res > 0.:  # file has -1.0 for NMR and other non x-ray methods
				resolutions[pdbID] = res

	# orthologs: distionary: {s. cer ORF : orthologue species (S.bay,S.mik or S.par) : ortholog ORF}
	orthologs = read_yeast_orthologs(externalData)
	
	return (ints, blastResults, pdbFullSeqs, aaSeq,structs, pdbLengths, resolutions, orthologs)


def read_blast_results(blastResultsPath):
	"""Read the file containing the results of a Blast search.
	Check the format of te BLAST result file. Here, trailing "_mRNA" after each ORF name, strip it to match ORF format so far
	Args:
		blastResultsPath (str): path to the file.

	Returns:
		res {(str, str, str): {str: float}}: maps ORF ID, PDB ID, chainID to statistic name to value
	"""
	with open(blastResultsPath, 'r') as file:
		res = {}
		for line in file:
			yeastORF = line.split()[1][:-5] # Strip the trailing "_mRNA"
			pdbID, chainID = line.split()[0].split('_')
			key = (yeastORF, pdbID, chainID)
			res[key] = {}
			res[key]['identity'] = int(line.split()[9])
			res[key]['query_length'] = int(line.split()[2])
	return res

def load_whole_protein_pdb_sequences(dataDir):
	with open(os.path.join(dataDir, 'pdb_reduced.faa'), 'r') as fH:
		seqs = {tuple(k.split('_')): str(v.seq) for k, v in
				SeqIO.to_dict(SeqIO.parse(fH, 'fasta')).items()}
	return seqs

def read_yeast_orthologs(externalData):
	"""
	Maps each S.cer ORF to ortholog ORFs in 8 other yeast species
	Returns:
		{str: {str: str}}: s. cer ORF --> species --> ortholog ORF with the best evalue in the BLAST alignment
	"""
	orthologs = defaultdict(dict)
	for species in ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
			with open(os.path.join(externalData, 'Scer-' + species + '-blast_stats_coverageCutoff.best.txt'),'r') as f:
				for line in f:
					# If want all orthologues in a list:orthologs[line.split()[0]].setdefault(species.lower(),[]).append(line.split()[1]) #/!\ + delete the filtering process below
					# If want only the orthologue with the best evalue:
					orthologs[line.split()[0]].setdefault(species.lower(),[]).append((line.split()[1],float(line.split()[2])))

	# Sort the species list by evalue, smallest first keep only the smallest one
	for scer_orf in orthologs:
		for orth_species in orthologs[scer_orf]:
			(orthologs[scer_orf][orth_species]).sort(key=lambda x: float(x[1]))
			orthologs[scer_orf][orth_species]= orthologs[scer_orf][orth_species][0][0]
	
	return orthologs

def select_best_pdb_matches(matchedStructures, blastResults,
							pdbLengths, resolutions, yeastPPIs):
	"""Pick the best homology matches in the PDB for each yeast PPI

	Decision is made by:
		1) Only consider PPIs for which both partners have:
		- At least one matched PDB structure 
		- with resolution of 3 Angstroms or better
		- and at least 2 chains in the file homology matched to yeast proteins.
	
		2) Construct the composite coverage of each structure as sum of the coverage for each partner protein.
		For each partner protein, the coverage by the structure is:
			the fraction of the yeast protein with identical residues in the 
			blast alignment * coverage of the PDB protein by the PDB structure
		The composite coverage of the PPI by the structure is:
			coverage of partner 1 + coverage of partner 2 by the PDB structure
			
		3) Take all PDBs within 10% of the maximum value of this composite coverage.

		5) Out of those, pick the structure with the best resolution.

	Args:
		matchedStructuress ({str: list((str, str))}): yeast ORF --> PDB chains homology matched to it
		blastResults ({(str, str, str): list(float)}): results of Blast of yeast ORFs against PDB.
			Keys are s. cer. ORF, PDB ID, chain ID.
		resolutions ({str: float}): resolution of each PDB structure.

	Returns:
		{(str, str): (str, str)}: yeast PPI --> best matched structure.
	"""

	bestMatch = {}
	pdbs = [i for j in matchedStructures.values() for i in j] # list of (PDBid,Chainid)s

	chains = {}
	for pdbID, chainID in pdbs:
		chains.setdefault(pdbID, set()).add(chainID)
	
	nChainsPerStruct = {pdbID: len(cIDs) for pdbID, cIDs in chains.items()}
	coverage = {k: float(structureLength) / float(seqLength) for k, (structureLength, seqLength) in pdbLengths.items()}
	
	# PPIs for which both partners have at least one matched structure
	PPIs = list(set([tuple(sorted([orf1,orf2])) for orf1 in yeastPPIs for orf2 in yeastPPIs[orf1] if (orf1 in matchedStructures and orf2 in matchedStructures)]))

	PPI_PDB_match = {}
	for (orf1,orf2) in PPIs:
		
		pdbs_orf1 = list(filter(lambda x: resolutions[x[0]] <= 3.0, matchedStructures[orf1])) 	# Resolution filter
		pdbs_orf2 = list(filter(lambda x: resolutions[x[0]] <= 3.0, matchedStructures[orf2])) 	# Resolution filter
		pdbs_orf1 = list(filter(lambda x: nChainsPerStruct[x[0]] > 1, pdbs_orf1)) 				# Numb of chains filter
		pdbs_orf2 = list(filter(lambda x: nChainsPerStruct[x[0]] > 1, pdbs_orf2)) 				# Numb of chains filter
		
		PPI_matches =[(pdb1,pdb2) for pdb1 in pdbs_orf1 for pdb2 in pdbs_orf2 if pdb1[0]== pdb2[0] and pdb1[1]!= pdb2[1]]
		
		if len(PPI_matches) == 0:
			continue

		compositeCoverage = {(pdb1,pdb2): 
				((float(blastResults[(orf1,) + pdb1]['identity']) /float(blastResults[(orf1,) + pdb1]['query_length'])) * coverage[pdb1]) +
				((float(blastResults[(orf2,) + pdb2]['identity']) /float(blastResults[(orf2,) + pdb2]['query_length'])) * coverage[pdb2]) 
				for (pdb1,pdb2) in PPI_matches}

		threshold = max(compositeCoverage.values()) - 0.01
		filtered = [k for k, v in compositeCoverage.items() if v > threshold]
		best = min(filtered, key=lambda x: resolutions[x[0][0]]+resolutions[x[1][0]])

		PPI_PDB_match.update({(orf1,orf2) : best})

	
	print (len(PPI_PDB_match), "/", len(PPIs) ,"PPIs between ",len(set(list(sum(PPI_PDB_match.keys(), ())))),"S.cere proteins have been matched to a structure")
	return PPI_PDB_match


def align_with_structures(alignmentFile, yeastPPIs, structs,bestMatch, aaSeq, resolutions,
						  pdbWholeProteinSeqs, orthologs, path_to_dssp, outPath):

	""" Build models using yeast PPI names (partner 1, partner 2) as keys (i.e. multiple PDB structures can be used per protein)
	 	If a protein is involved in more than one PPI using more than one PDB structure,
	 	this protein can be modeled multiple times using different PDBs
	 	
	Args:
		alignmentFile (str): file with the alignemnt to 3 other yeast species
		TO DO: add details of input output

	Returns: 
		(orfName, (PDBid,PDB chain)) --> model.
	"""
	if not os.path.exists(outPath):
		# Load alignment with other yeast species
		strains = ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']
		yeastAlignment = load_aligned_yeast_codons(alignmentFile)

		# Get the list of all of the (orfName, (PDBid,PDB chain)) matches
		counts = Counter([(orf[0],bestMatch[orf][0]) for orf in bestMatch ]+[(orf[1],bestMatch[orf][1]) for orf in bestMatch ])

		alignment_keys= list(yeastAlignment.keys())
		count_keys = [y[0] for y in list(counts.keys())]

		matchedChains, intProps = {}, {}
		pdbAligned = {}
		proteins = {}
		print ('\n Building structural models of', len(counts), 'yeast proteins')
		
		count_skipped_no_residue = 0
		count_no_modeled_interface =0

		for y in tqdm(counts.keys()):

			if y[0] not in yeastAlignment:
				continue   # Throw away proteins without yeast orthologs
			
			pdbID = y[1][0]
			chainID = y[1][1]
			print ("\n",y[0]," ", pdbID," Chain: ",chainID) 

			(complex, iP, resInfo) = model_interfaces(y[0], pdbID, chainID, structs,yeastPPIs[y[0]],path_to_dssp)
			
			if len(complex) <= 1:
				count_no_modeled_interface = count_no_modeled_interface + 1

			matchedChains[y] = complex
			if complex != chainID:
				intProps[y] = iP
			
			seqPDB = ''.join([r['aa'] for r in resInfo.values()])

			# Add the ID notation (identifier for the residue in the PDB file - 3LATOM notation. Ex: MET1:A) to the values in resInfo
			conversion = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
						 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
						 'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
						 'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
			for i in resInfo:
				aa = conversion[resInfo[int(i)]['aa']]
				idNotation = aa+str(i)+":"+chainID
				resInfo[int(i)].update({'id':idNotation})

			try:
				pdbAligned[y] = align_yeast_to_pdb(aaSeq[y[0]],
												   resInfo.values())
				
			except UserWarning as e:
				# since pdb structure can only cover a fraction of the protein
				#  alignment working is not guarenteed
				print (e)
				print ('problem aligning', y[0], 'and', pdbID, chainID)
				print ('seq1:', aaSeq[y[0]], '\nseq2:', seqPDB)
				continue

			proteins[y] = Protein(y[0], pdbID, chainID, aaSeq[y[0]], resolutions[pdbID])
			proteins[y].pdbComplex = complex
			proteins[y].pdbStructureSeq = seqPDB
			proteins[y].pdbWholeProteinSeq = pdbWholeProteinSeqs[(pdbID, chainID)]
			proteins[y].orthologs = orthologs[y[0]]

			if y in intProps:
				for intf in intProps[y]:
					proteins[y].add_interface(intf['partner'], intf['chain'])
			if len(list(pdbAligned[y])) != len(list(yeastAlignment[y[0]])):
				raise UserWarning('Structure and yeast alignments differ')

			# Index is the residue in the yeast sequence saved as proteins.seq, save it as a property for the residue when adding it?
			for index, (algnYstCdn, algnPDBRes) in enumerate(zip(yeastAlignment[y[0]],
																 pdbAligned[y])):
				# throw away gaps in any of the alignments
				if any(['-' in algnYstCdn[(i, 'nt')] for i in strains]):
					continue
				if '-' in algnPDBRes:
					continue
				cdns = {k1: v for (k1, k2), v in algnYstCdn.items() if k2 == 'nt'}
				res = Residue(cdns,
							  algnPDBRes['sasa'],
							  algnPDBRes.get('dsasa', 0.),
							  algnYstCdn[('Scer', 'aa')])

				proteins[y].add_res(res, algnPDBRes.get('interface', None))
				
				proteins[y].residues[-1].id = algnPDBRes.get('id', '')		# Identifier notation for the residue to match with dssp output: format examples MET1:A (AaPosition:chain)
				proteins[y].residues[-1].indexInSeq = index 				# Identifier notation for the residue to match with the sequence saved as proteins.seq
				proteins[y].residues[-1].contacts = algnPDBRes.get('contacts', np.nan)
				proteins[y].residues[-1].distance_from_center = algnPDBRes.get('distance_from_center', np.nan)
				proteins[y].residues[-1].distance_from_edges = algnPDBRes.get('distance_from_edges', np.nan)		
				
				if (len(set(proteins[y].residues[-1].codons.values()))== 1) :
					proteins[y].residues[-1].substitution = 0
				else:
					proteins[y].residues[-1].substitution = 1

			if len(proteins[y].residues) == 0:
				# sometimes large gaps in the alignment with the PDB protein
				# and one of the other yeast proteins will cover everything
				count_skipped_no_residue = count_skipped_no_residue +1
				del proteins[y]
				continue 
		
		print ("No residues in model for ",count_skipped_no_residue, "/",len(counts)," PPIs.")
		print ("No interface modeled",count_no_modeled_interface, "/",len(counts)," PPIs. ")	
		pickle.dump(proteins, open(outPath, "wb"))

	else:
		print ("Structural models already built. Skipping.")

	return


def load_aligned_yeast_codons(codonAlignmentPath):
	"""
	Load codon alignment of yeast orthologs.
	Args:
		codonAlignmentPath: input file path

	Returns:
		{str: {(str, str): list(str)}}}: s.cer. ORF -> (strain, aa/nt) -> seq
	"""
	yeastAlignment = {}
	with open(codonAlignmentPath, 'r') as f:
		for l in f:
			scerORFID = l.split()[0]
			if scerORFID not in yeastAlignment:
				yeastAlignment[scerORFID] = {}
			strainName, seqType = l.split()[1].split('.')
			yeastAlignment[scerORFID][(strainName, seqType)] = l.split()[2:]

	# Turn it into a more convenient format
	algn = {}
	(total,inconsistent) = (0,0)
	for orf in yeastAlignment:
		total = total +1
		numAA = len(yeastAlignment[orf][('Scer', 'aa')])
		# Verifying that the lengths are consistent accross the 6 strains
		if any(map(lambda x: len(yeastAlignment[orf][(x, 'nt')]) != numAA,['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb'])):
			raise UserWarning('Inconsistent sequence lengths in:' + codonAlignmentPath)
			inconsistent = inconsistent +1

		keys = [('Scer', 'aa'), ('Scer', 'nt'), ('Spar', 'nt'),
				('Sbay', 'nt'), ('Smik', 'nt'),('Ncas','nt'),('Cgla','nt'),('Agos','nt'),('Klac','nt'),('Calb','nt')]
		algn[orf] = list(map(lambda i: {k: yeastAlignment[orf][k][i] for k in keys}, range(numAA)))

	return algn


def model_interfaces(yeastORF, pdbID, chainID, structs, yeastPPIs,path_to_dssp):
	"""Form a structural model of a yeast protein in complex with others.

	For a yeast protein with structural model, use interaction data to match
	to structural models of partner yeast proteins in the same pdb complex.

	Args:
		yeastORF (str): s. cerevesiae orf ID
		pdbID (str): PDB structure ID of model of y
		chainID (str): PDB chain ID of model of y
		structs (dict): map s. cer. protein ID to list of matched structures
		yeastPPIs (list): list of interactions

	Returns:
		str: IDs of PDB chains that are in complex of interaction partners
		dict(str: str): yeast ORF ID and PDB chain ID of modeled partner
		str: amino acid sequence of PDB structure

	"""
	complex = chainID
	intProps = []

	for interactionPartner in yeastPPIs:
		if yeastORF == interactionPartner:
			continue
		if interactionPartner not in structs:
			continue
		for (pdbIDpartner, chainIDpartner) in structs[interactionPartner]:
			if pdbID == pdbIDpartner and chainID != chainIDpartner:
				if chainIDpartner not in complex:
					complex = complex + chainIDpartner
					intProps.append({'partner': interactionPartner,
									 'chain': chainIDpartner})
	
# # For partner yeast proteins of yeastORF in the same pdb complex:
	# Calc_sasa: 
	# 	 - Assign residue ID to each residue in the complex 
	# 	 		(self.id: Identifier notation for the residue in the PDB file (to match with dssp output: format examples MET1:A (AaPosition:chain)))
	# 	 - Assign aa notation to the residue
	# 	 - Assign monomer solent accessibility value to the residue 
	# 			(rsa, computed with dssp)
	# 	 - Assign delta RSA value to the residue 
	# 			(dRsa, computed with dssp)
	resInfo, redComplex = calc_sasa(pdbID, chainID, complex,path_to_dssp)
	
	# Interface_geometry: Only if len(redComplex) > 1
	# 	 - Assign residue contact number value to each residue in the complex 
	# 	 		(self.contacts: number of residue contacts with interaction partner
	# 	 - Assign distance to the geometric center of the interface to each residue in the complex 
	# 			(distance_from_center: numerical value)
	# 	 - Assign distance from edges of the interface to each residue in the complex 
	#			(distance_from_edges: numerical value)

	if len(redComplex) > 1:
		interface_geometry(pdbID, chainID, redComplex[1:], resInfo)

	redIntProps = []
	for ip in intProps:
		if ip['chain'] in redComplex:
			redIntProps.append(ip)

	return redComplex, redIntProps, resInfo


def calc_sasa(pdbID, chain, complex,path_to_dssp):
	"""Per-residue RSA and deltaRSA for a PDB chain in a complex.
	RSA is the fraction of solvent-accesible surface area relative to maximum.
	deltaRSA is the increase in RSA when the interaction partner is excluded.
	Assigns a residue to a PPI interface based on the partner which causes
	the largest change in RSA.
	Also filters the input complex, retaining only chains that change the RSA
	of at least one residue in the chain of interest.

	Args:
		pdbID (str): PDB ID of structure
		chain (str): chain ID of chain
		complex (str): chains to consider when calculating solvent
					   accessibility. For just the individual chain,
					   complex == chain.

	Returns:
		str: AA sequence of chain in PDB structure (NOTE: not necessarily the
			 whole protein sequence is in the structure)
		list(float): RSA per residue
		list(float): deltaRSA per residue
		list(str): interface number of each residue
				   (0 means not assigned to an interface)
		str: chain IDs of filtered complex
	"""
	# Threshold at which to consider change in RSA non-zero
	minNonZeroDeltaRSA = 0.001
	inDir = '../../data/processed/assemblies/'
	struct = load_structure(inDir + pdbID + '.pdb1')
	# Change in RSA for each interface
	# Run wrapper for DSSP
	redComplex = complex[0]
	if len(complex) > 1:
		dRSAPerPartner = []
		for partner in complex[1:]:
			sa = sasa_scan(struct, chain, chain + partner,path_to_dssp, True)[0]
			if any([r['dsasa'] >= minNonZeroDeltaRSA for r in sa.values()]):
				redComplex += partner
				dRSAPerPartner.append([r['dsasa'] for r in sa.values()])
	
	sasa = sasa_scan(struct, chain, complex,path_to_dssp, True)[0]
	resInfo = OrderedDict([(k, round_delta_RSA(v, minNonZeroDeltaRSA))
						   for k, v in sasa.items()])

	# Assign residue to the interaction which gives largest change in RSA
	if len(redComplex) > 1:
		for i, res in enumerate(resInfo.values()):
			if res['dsasa'] > 0.:
				deltas = [j[i] for j in dRSAPerPartner]
				res['interface'] = redComplex[1:][deltas.index(max(deltas))]
	
	return resInfo, redComplex

def round_delta_RSA(residue, threshold):
	"""Round down values of delta RSA below a threshold to zero.
	Args:
		residue (dict)
		threshold (float): minimum change in RSA considered non-zero

	Returns:
		dict
	"""
	if residue['dsasa'] >= threshold or residue['dsasa'] == 0.:
		return residue
	else:
		return {'aa': residue['aa'], 'sasa': residue['sasa'], 'dsasa': 0.0}

def interface_geometry(pdbID, chainID, partnerChainIDs, resInfo):
	"""Quantities of interest for interfacial residues.
	Quantities:
	- rsa: 	relative solvent accessibility in uncomplexed state
	- dRSA: change in rsa upon complex formation
	- Residue contacts: Number of residues in the partner protein within 10 Angstroms from our residue
	- Distance from center of the interface: distance from the geometric center of the interface. 
	- Distance from edges of the interface: defined as 95th percentile of largest complex RSA
	 
	Args:
		pdbID (str)
		chainID (str)
		partnerChainIDs (str)
		resInfo (OrderdDict(int: dict)): PDB residue ID -> solvent accessibility values of residue.
										Used to define interface residues. 

	"""
	# Protein A & B are partners
	# Protein A
	if len(partnerChainIDs) == 0:
		raise ValueError('Expected at least one partner.')
	parser = PDBParser(PERMISSIVE=1, QUIET=1)
	ipath = '../../data/processed/assemblies/' + pdbID + '.pdb1' 
	structure = parser.get_structure('name', ipath)[0]
	protein = structure[chainID]
	nContacts = {} 
	intf_dsasa = [] # List of all values of complex RSA, just used to plot/ get the 95th percentile of the distribution

	# All partner chains in protein B
	for partnerChainID in partnerChainIDs:
		partner = structure[partnerChainID]
		intf = []
		nContacts[partnerChainID] = np.zeros(len(protein), dtype=int)
		cutoff = 10.0  # Angstroms cutoff to define interface contacts

		# --- Residue contacts ---
		# For all residues in the protein, 
		for i, residue in enumerate(protein):
			resID = residue.id[1]
			if resID not in resInfo:
				continue  # Skip small number of residues not parsed by dssp
			
			resInfo[resID]['contacts'] = 0 # initialize to 0
			dists = list(map(lambda r: residue['CA'] - r['CA'], partner))
			if resInfo[resID]['dsasa'] > 0.0:
				
				resInfo[resID]['contacts'] += np.sum([d < cutoff for d in dists])
				
			if ('interface' in resInfo[resID] and resInfo[resID]['interface'] == partnerChainID):
				intf.append(residue)
				# Append both residue and dsasa to be able to easily transform it to a list of residues (most exposed surface)
				intf_dsasa.append((residue,resInfo[resID]['dsasa']))

		# --- Distance from center/edges of the interface  ---		
		if len(intf) < 10: # Keep looping over all partner chains until we have a big enough interface!
			continue

		# Get the dsasa cutoff to use to define the outer boundary of the interface (95th percentile of largest RSAs complex)
		dsasas = [x[1] for x in intf_dsasa]
		percentile = np.percentile(dsasas, 95)
		# Outer boundaries / edges of the interface
		outerBoundaries = [x[0] for x in intf_dsasa if x[1] >= percentile]

		# Compute the geometrical mean of the interface
		intfCenter = np.mean([r['CA'].coord for r in intf], axis=0) 	# x,y,z coords of geometric center of the interface
		intfDist = [np.sqrt(np.sum((r['CA'].coord - intfCenter) ** 2)) 	# list of average distance of all interfacial residues to the geometric center of the interface
					for r in intf]
		
		# For all residues in the protein, 
		for i, residue in enumerate(protein):
			resID = residue.id[1]
			if resID not in resInfo:
				continue
			if 'interface' not in resInfo[resID]:
				continue
			if resInfo[resID]['interface'] == partnerChainID:
				d = np.sqrt(np.sum((residue['CA'].coord - intfCenter) ** 2))		# Distance of residue to geometric center
				e = min(map(lambda r: residue['CA'] - r['CA'], outerBoundaries))	# Distance of residue to outer boundary
				
				# Assign  distance to geometric center
				resInfo[resID]['distance_from_center'] = d

				# Assign  distance to geometric center
				resInfo[resID]['distance_from_edges'] = e

	return

def align_yeast_to_pdb(scerAASeq, pdbProtein):
	"""Align RSA/deltaRSA of PDB complex to yeast protein.

	Args:
		scerAASeq (str): aa sequence of s. cer. protein
		pdbProtein (list(dict)): residue-level data of homolog in the PDB,
								 including AA sequence, solvent accessibility
								 etc.

	Returns:
		list(dict): alignment of RSA and other values of PDB chain to yeast protein

	"""

	seqPDB = ''.join([r['aa'] for r in pdbProtein])
	iAlnd = align(scerAASeq, seqPDB, [str(i) for i in range(len(seqPDB))])
	
	return [i if '-' in i else list(pdbProtein)[int(i)] for i in iAlnd]


def main():
	
	externalData = '../../data/external'
	processedData = '../../data/processed'
	
	# Download DSSP
	path_to_dssp = download_dssp()

	# Load input files needed
	print ("----------- Loading input -----------")
	(yeastPPIs,matchStats, pdbFullSeqs, aaSeq,structs,pdbLengths, resolutions, orthologs) = load_inputs(externalData,processedData)
	
	# For testing purposes, build models for a subset of the structures (small_structs)
	small_structs ={k: structs[k] for k in list(structs.keys())[:40]}

	# Select best PDB matches
	print ("----------- Selecting best PDB matches -----------")
	bestMatch = select_best_pdb_matches(structs,
										matchStats,
										pdbLengths,
										resolutions,
										yeastPPIs)
	# Align with PBD structures, model with one PDB file per PPI
	print ("----------- Aligning with structures -----------")
	
	alignmentFile =os.path.join(externalData,'codon_alignment.txt')
	align_with_structures(alignmentFile, yeastPPIs, structs, bestMatch, 
						aaSeq, resolutions, pdbFullSeqs, orthologs, path_to_dssp,
						outPath=os.path.join(processedData,'modeled_proteins.pkl'))

if __name__ == '__main__':
	main()


