#!/usr/bin/env python

""" 
Step 3 in the pipeline - 
Align RSA/dRSA from structural models to yeast protein sequences.


Output files:
	modeled_proteins.pkl - structurally modeled proteins

Input files: 
	Spom.aa / Sjap.nt etc. - sequences of yeast ORFs
	scer_protein_interactions_multiple_reports.txt - BioGrid physical s. cer interactions reported more than once
	blast_scer_vs_pdb.txt - blast results of s. sce proteins against the PDB
	pdb_reduced.faa - AA sequences of PDB chains
	curated_scer_pdb_homologs.txt - map of s. scer proteins to homologous PDB chains
	pdb_lengths.csv - length of each PDB structure
	pdb_chain_uniprot.csv - mapping of PDB chains to UniProt proteins

Author Leah Pollet
pollet.leah@gmail.com
"""

import os
import re
import csv
import random
import pickle
from tqdm import tqdm
from  ast import literal_eval
from collections import defaultdict,namedtuple,Counter,OrderedDict
from Bio import SeqIO
from Bio.PDB import PDBParser
import pandas as pd
import tarfile
from modeled_PPIs import PPI,Protein,Interface,Residue
from sasa_scan import load_structure, sasa_scan
import numpy as np
import matplotlib.pyplot as plt
import itertools
from pattern_align_local import align


def load_inputs(externalData,processedData,blastToPDBData,related_species):
	"""
	Load all of the data.
		Args: 
			externalData:  directory with the unprocessed data from previous analysis
			processedData: directory with the processed data from previous analysis
		
		Outputs:
			ints: 			Dictionary of PPIs in this analysis. Format: {orfName: set(orfName)} with orfName all orfs in S.cerevisiae PPIs (not double sided but there could still be duplicates...)
			orfs:			Set of all unique orfs in S. cerevisiae PPIs. Format: set(orf) with orf the name of orfs in S. cerevisae
			blastResults: 	Stats from the BLAST homology mapping process. Format: {(orfID, PDBid, ChainID): {statName: statValue}} with orfID = S.cerevisiae orfs, statName = identity and query length from BLAST with PDB
			pdbFullSeqs:	Dictionary with aa sequences for proteins in PDB. Format {PDBID:PDBsequence}
			aaSeq:			Dictionary with aa sequences for proteins in S. cerevisiae. Format: {yeastORFname: aa sequence}, with yeastORFname orfs in S. cerevisiae PPIs
			structs:		Dictionary with the mapping between ORF names and PDB structure and chain names. Format: { yeastORFname : [(struct, chain), (,)...]} with yeastORFname = orfs in S. cerevisiae PPIs, struct = e.g: 5fpn, chain = e.g: A
			pdbLengths:		Dictionary with the length of the structure and sequence for each chain. Format:  {(pdbID, chainID) : (length of struture, length of sequence)}
			resolutions:    Dictionary with the resolution of each PDB structure. Format: {pdbID: resolution}
			orthologs_spom/orthologs_scer:	Dictionary with the mapping between ORF name in S. pom and ORF name in orthologue species (S. jap, S.oct, S.cry etc.)
	"""
	
	# ints: Double sided dictionary of all S.cer PPIs
	intsList = [[l.split()[0], l.split()[1]] for l in open(os.path.join(processedData,'scer_protein_interactions_multiple_reports.txt'), 'r')]
	ints = defaultdict(set) 
	orfs = set()	
	for y1, y2 in intsList:
		orfs.update([y1,y2])
		ints[y1].add(y2)
		ints[y2].add(y1)
	
	# blastResults: Stats from homology mapping
	blastResults = read_blast_results([os.path.join(blastToPDBData,'blast_scer_vs_pdb.txt')])
	
	# pdbFullSeqs: Sequences in FASTA format
	pdbFullSeqs = load_whole_protein_pdb_sequences(blastToPDBData)
	
	# aaSeq:	
	aaSeq = {}
	with open(os.path.join(blastToPDBData, 'Scer_reduced.aa'), 'r') as fH:
		aaSeq.update({k:str(v.seq) for k,v in SeqIO.to_dict(SeqIO.parse(fH, 'fasta')).items()})  

	# structs: 	
	structs = defaultdict(dict)
	f_structs_scer = open(os.path.join(processedData,'curated_scer_pdb_homologs.txt'), 'r')
	for l in f_structs_scer:
		yeastORFname,struct,chain = l.split()[0],l.split()[1],l.split()[2]
		structs.setdefault(yeastORFname,[]).append((struct, chain))
	
	# pdbLengths: 
	pdbLengths = {}
	with open(os.path.join(processedData,'pdb_lengths.csv'), 'r') as fCov:
		for l in fCov.readlines()[1:]:
			l = l.strip('\n') 	# /!\ NT strip the trailing "\n "at the end of lSeq
			pdbID, chainID, lStruct, lSeq = l.split(',')
			pdbLengths[(pdbID, chainID)] = (lStruct, lSeq)
	
	# resolutions: 
	resolutions = {}
	with open(os.path.join(externalData,'resolu.idx'), 'r') as fResolution:
		for line in fResolution.readlines()[6:]:
			pdbID = line.split(';')[0].strip().lower()
			data_resolution = line.split(';')[1].strip()
			if data_resolution == '':
				continue
			res = float(data_resolution)
			# if res > 0.:  # file has -1.0 for NMR and other non x-ray methods *** TOCHECK
			resolutions[pdbID] = res
	
	# orthologs: dictionary: {s. cer ORF : orthologue species  : ortholog ORF}
	dirList = [os.path.join(externalData,"blastToCloselyRelated"),os.path.join(externalData,"blastToAdditionalCloselyRelated"),os.path.join(externalData,"blastToAdditionalDistantlyRelated")]
	orthologs = read_yeast_orthologs(dirList,related_species,"Scer")

	return (ints,orfs, blastResults, pdbFullSeqs, aaSeq,structs, pdbLengths, resolutions, orthologs)


def read_blast_results(blastResultsPaths):
	"""
	Read the file(s) containing the results of Blast search(es).
	Args:
		blastResultsPath (list): list of path to blast results files (only a single path in the analysis)

	Returns:
		{(str, str, str): {str: float}}: maps (ORF ID, PDB ID, chainID) to statistic name to value
	"""
	res = {}
	for blastResultsPath in blastResultsPaths:
		with open(blastResultsPath, 'r') as file:
			for line in file:
				yeastORF = line.split()[1]
				pdbID, chainID = line.split()[0].split('_')
				key = (yeastORF, pdbID, chainID)
				res[key] = {}
				res[key]['identity'] = int(line.split()[9])
				res[key]['query_length'] = int(line.split()[2])
				res[key]['eval'] = line.split()[12]
				
	return res

def load_whole_protein_pdb_sequences(dataDir):
	with open(os.path.join(dataDir, 'pdb_reduced.faa'), 'r') as fH:
		seqs = {tuple(k.split('_')): str(v.seq) for k, v in
				SeqIO.to_dict(SeqIO.parse(fH, 'fasta')).items()}
	return seqs

def read_yeast_orthologs(dirs_list,related_species,main_species):
	"""
	Maps each main_species ORF to ortholog ORFs in each closely_related_species
	Args:	
		dirs_list (lst(str)): list of directories where the alignments to closely related species are located
		related_species (lst): list of related species 
		main_species(str)	: name of the main species ('Scer')
	Returns:
		{str: {str: str}}: ORF --> species --> ortholog ORF with the best evalue in the BLAST alignment
	"""
	# INSTEAD OF EXTERNAL DATA, HAVE A CHECK IF FILE EXISTS IN ALL POSSIBLE DIRS
	orthologs = defaultdict(dict)
	# Read in all possible orthologs
	for species in related_species:
		for dir in dirs_list:
			path = os.path.join(dir, main_species+'-' + species + '-blast_stats_coverageCutoff.best.txt')
			if os.path.exists(path):
				with open(path,'r') as f:
					for line in f:
						orthologs[line.split()[1]].setdefault(species,[]).append((line.split()[0],float(line.split()[2])))
						
	# Sort the species list by evalue, smallest first keep only the smallest one
	for orf in orthologs:
		for orth_species in orthologs[orf]:
			(orthologs[orf][orth_species]).sort(key=lambda x: float(x[1]))
			orthologs[orf][orth_species]= orthologs[orf][orth_species][0][0]
	
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
		blastResults (dict) Stats from the BLAST homology mapping process. Format: {(orfID, PDBid, ChainID): {statName: statValue}} with orfID = S.cerevisiae orfs, statName = identity and query length from BLAST with PDB
		resolutions ({str: float}): resolution of each PDB structure.

	Returns:
		PPI_PDB_match (dict): best homology mapped structure for each PPI. Format: {(orfA, orfB): ((pdbIDA, pdbChainA),(pdbIDB, pdbChainB))} with orfA/orfB interacting proteins in S. cerevisiae , and PDBIDA,PDBChainA the best homology mapped structure for orfA
	"""
	# CHECK AND REMOVE SINGLE PROTEINS
	# PRINT A SUMMARY FILE WITH THE HOMOLOGY MAPPED STRUCTURE FOR PPIS TO processed/scer_PPIs_PDB_homologs.txt
	# Use the print summary and edit it

	pdbs = [i for j in matchedStructures.values() for i in j] # list of (PDBid,Chainid)s


	chains = {}
	for pdbID, chainID in pdbs:
		chains.setdefault(pdbID, set()).add(chainID)
	
	nChainsPerStruct = {pdbID: len(cIDs) for pdbID, cIDs in chains.items()}
	coverage = {k: float(structureLength) / float(seqLength) for k, (structureLength, seqLength) in pdbLengths.items()}

	# PPIs for which both partners have at least one matched structure (alphabelically sorted + set to remove duplicated)
	PPIs = list(set([tuple(sorted([orf1,orf2])) for orf1 in yeastPPIs for orf2 in yeastPPIs[orf1] if (orf1 in matchedStructures and orf2 in matchedStructures)]))
	
	PPI_PDB_match = {}

	# Finding the best match for PPIs
	for (orf1,orf2) in tqdm(PPIs):
		
		pdbs_orf1 = list(filter(lambda x: resolutions[x[0]] <= 3.0, matchedStructures[orf1])) 	# Resolution filter
		pdbs_orf2 = list(filter(lambda x: resolutions[x[0]] <= 3.0, matchedStructures[orf2])) 	# Resolution filter
		pdbs_orf1 = list(filter(lambda x: nChainsPerStruct[x[0]] > 1, pdbs_orf1)) 				# Numb of chains filter
		pdbs_orf2 = list(filter(lambda x: nChainsPerStruct[x[0]] > 1, pdbs_orf2)) 				# Numb of chains filter

		# PDB structures that can be used to model the PPI (have both partners structured, and good enough resolution)
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

		
	ppis_matched = {k:v for (k, v) in PPI_PDB_match.items() if None not in k} # This was for 1gle prots
	print (str(len(PPI_PDB_match)), "/", str(len(PPIs)) ,"PPIs between",str(len(set(list(sum(PPI_PDB_match.keys(), ())))))," proteins have been matched to a structure")
	
	ofile= "../data/processed/PPIs_to_model.txt"
	PPIs_to_model,orfs_in_PPIs_to_model=print_summary(PPI_PDB_match,ofile,blastResults,coverage,compositeCoverage)

	
	return PPI_PDB_match,PPIs_to_model,orfs_in_PPIs_to_model

def print_summary(bestMatch,ofile,blastResults,coverage,compositeCoverage):
	"""
		Print summary of the numbers we have so far
		# Edit so we have: orf1, orf2, number of reports, best match pdb, eval, compositecoverage, orf1 chain, orf2 chain, coverage, coverage
		# 
	"""
	PPIs_to_model = []
	PPI_summary = namedtuple('PPI', ["ORF_A","ORF_B","NumReports","PDB", "chain_A","chain_B","e_value_A","e_value_B","coverage_A","coverage_B"])
	
	# Read the number of reports for each PPI
	reportsFile = os.path.split(ofile)[0]+"/scer_protein_interactions.txt"
	reportsNum = {(line.split()[0],line.split()[1]):int(line.split()[2]) for line in open(reportsFile) }

	for (orfA,orfB) in bestMatch:

		pdb_ID = bestMatch[(orfA,orfB)][0][0] #= bestMatch[(orfA,orfB)][1][0]
		chainA = bestMatch[(orfA,orfB)][0][1]
		chainB = bestMatch[(orfA,orfB)][1][1]
		ppi = PPI_summary(str(orfA),str(orfB),str(reportsNum[(orfA,orfB)]),str(pdb_ID),str(chainA),str(chainB),str(blastResults[(orfA,pdb_ID,chainA)]['eval']),str(blastResults[(orfB,pdb_ID,chainB)]['eval']),str(coverage[(pdb_ID,chainA)]),str(coverage[(pdb_ID,chainB)]))
		PPIs_to_model.append(ppi)


	# Write to the output file
	if not os.path.exists(ofile):
		with open(ofile, "w" ) as f:
			f.write('\t'.join(["ORF_A","ORF_B","NumReports","PDB", "chain_A","chain_B","e_value_A","e_value_B","coverage_A","coverage_B"])+"\n")
			for ppi in set(PPIs_to_model):
				f.write('\t'.join([i for i in ppi])+"\n")
	else:
		print(ofile,"already exists. Skipping.")

	ORFs = set([ppi.ORF_A for ppi in PPIs_to_model]+[ppi.ORF_B for ppi in PPIs_to_model])
	
	return PPIs_to_model,ORFs


def download_dssp():
	# useful if dssp is not compling properly:
	# http://melolab.org/supmat/moma/install.html
	oPath = "../data"
	url = "https://github.com/cmbi/dssp/archive/2.0.4.tar.gz"
	# url = "ftp://ftp.cmbi.ru.nl/pub/molbio/software/dssp-2/dssp-2.0.4.tgz"

	zipped = oPath+"/"+os.path.basename(url)
	
	path_to_dssp = oPath+"/dssp-"+os.path.basename(url)[:5]
	
	path_to_dssp_exec = oPath+"/dssp-"+os.path.basename(url)[:5]+"/mkdssp"
	
	if not os.path.exists(path_to_dssp):
		print("Downloading DSSP")
		os.system('wget -P '+ oPath +' '+url)
		tar = tarfile.open(zipped)
		tar.extractall(oPath)
		tar.close()
		os.remove(zipped)	
		
	# subprocess.Popen(["./autogen.sh"], stdout=subprocess.PIPE, cwd=path_to_dssp)
	# subprocess.Popen(["./configure"], stdout=subprocess.PIPE, cwd=path_to_dssp)
	# subprocess.Popen(["make mkdssp"], stdout=subprocess.PIPE, cwd=path_to_dssp)
	
	print(path_to_dssp_exec)
	return path_to_dssp_exec

def align_with_structures_multiple(PPIsToModel,alignmentFile,yeastPPIs,structs,bestMatch, aaSeq, resolutions,
								   pdbWholeProteinSeqs, orthologs_mapping,path_to_dssp, outPath,ofile):

	""" Build models using yeast PPI names (partner 1, partner 2) as keys (i.e. multiple PDB structures can be used per protein)
	 	If a protein is involved in more than one PPI using more than one PDB structure,
	 	this protein can be modeled multiple times using different PDBs

	 														PPI(ORF_A, ORF_B, NumReports, PDB, chain_A, chain_B, e_value_A, e_value_B, coverage_A, coverage_B)
	Args:													PPI(ORF_A='YER148W', ORF_B='YER148W', PDB='1tbp', chain_A='B', chain_B='A', e_value_A='3.76e-130', e_value_B='3.76e-130', coverage_A='1.0', coverage_B='1.0')
		PPIsToModel (str): 								List of PPI objects with format namedtuple('PPI', ["ppi", "species","interolog","interologspecies","reportsInSpom","reportsInScer","HaiywanData","multipleReports","ppiType","PDB","PDBinterolog"])
		alignmentFile :								 	Path to the codon alignment file to closely related species for S. cerevisiae
		yeastPPIs (dict): 								NOT NEEDED??? Dictionary of PPIs in this analysis. Format: {orfName: set(orfName)} with orfName all orfs in S. pombe and S.cerevisiae PPIs (double sided)
		structs (dict): 								Dictionary with the mapping between ORF names and PDB structure and chain names. Format: { yeastORFname : [(struct, chain), (,)...]} with yeastORFname = orfs in S. pombe and S. cerevisiae PPIs, struct = e.g: 5fpn, chain = e.g: A
		bestMatch (dict): 								Best homology mapped structure for each PPI. Format: {(orfA, orfB): ((pdbIDA, pdbChainA),(pdbIDB, pdbChainB))} with orfA/orfB interacting proteins in S. cerevisiae or S. pombe, and PDBIDA,PDBChainA the best homology maped structure for orfA
		aaSeq (dict):									Dictionary with aa sequences for proteins in S.pombe and S. cerevisiae. Format: {yeastORFname: aa sequence}, with yeastORFname  orfs in S. pombe and S. cerevisiae PPIs
		pdbWholeProteinSeqs (dict): 					Dictionary with aa sequences for proteins in PDB. Format {PDBIS:PDBsequence}
		resolutions:    								Dictionary with the resolution of each PDB structure. Format: {pdbID: resolution}
		orthologs_mapping (dict):						Dictionary with the mapping between ORF name in S. cere and ORF name in related species (S. jap, S.oct, S.cry etc.) 
		path_to_dssp (str):								Path to the dssp executable
		outPath (str):									Path to the pickled final models
	
		
	"""
	if not os.path.exists(outPath):
		# (1) Load alignment to closely related species.
		yeastAlignment = load_aligned_yeast_codons(alignmentFile)
		

		modeledPPIs = {}
		count_no_modeled_interface =0
		error1,error2,error3,error4,error5,error6,count = 0,0,0,0,0,0,0

		# (2) Create PPI models 
		for ppiToModel in tqdm(PPIsToModel):
			
			PDBs = ((ppiToModel.PDB,ppiToModel.chain_A) , (ppiToModel.PDB,ppiToModel.chain_B))
			evalues = (ppiToModel.e_value_A,ppiToModel.e_value_B)
			coverages = (ppiToModel.coverage_A,ppiToModel.coverage_B)
			key = (ppiToModel.ORF_A,ppiToModel.ORF_B)
			
			# (3) Create a PPI object and add it to models. 
			modeledPPIs[key]= PPI(ppiToModel.ORF_A,ppiToModel.ORF_B,ppiToModel.NumReports)

			
			# (4) Loop through both proteins in the PPI to model (orf1 and orf2) 
			for orf1,orf2,pdbOrf1,pdbOrf2,evalue,coverage in zip([key[0],key[1]],[key[1],key[0]],[PDBs[0],PDBs[1]],[PDBs[1],PDBs[0]],[evalues[0],evalues[1]],[coverages[0],coverages[1]]):
				# Can't model a "None" (missing/lost) protein
				if orf1 == None: 
					continue

				orthologs = orthologs_mapping[orf1]
				
				# Get the structural info for the protein to model
				pdbSeqs,resolution = None,None
				if pdbOrf1 in pdbWholeProteinSeqs:
					pdbSeqs = pdbWholeProteinSeqs[pdbOrf1]
				else:
					print("can't model", orf1, "with matched PDB" ,pdbOrf1, "is not in the pdbWholeProteinSeqs dict")
					
				if pdbOrf1[0] in resolutions:
					resolution = resolutions[pdbOrf1[0]]
				else:
					print("can't model", orf1, "with matched PDB" ,pdbOrf1, "is not in the resolutions dict")
					
				# Create protein object, with modeled interface and residues
				
				modeled_protein= model_protein(orf1,orf2,'Scer',pdbOrf1,pdbOrf2,
												evalue,coverage,
												aaSeq[orf1],pdbSeqs,
												orthologs,resolution,
												path_to_dssp,yeastAlignment)
				
				# Check for errors in a protein model
				# Error messages:
				# error1: No sequence for orthologs in closely related species (needed for evolutionary rate calculations)
				# error2: Coverage is too small on the PDB structure (cannot align sequence and structure)
				# error3: Alignment of the protein sequence with the PDB structure and alignemnt of the protein sequence with orthologs in closely related species have different lengths
				# error4: The protein has no modeled residues

				if not (modeled_protein in ["error1","error2","error3","error4","error5","error6"]):
					if not modeled_protein in modeledPPIs[key].proteins:
						modeledPPIs[key].add_protein(modeled_protein)
						count = count+1
					else:
						continue
				else:
					for error,cnt in zip(["error1","error2","error3","error4","error5","error6"],[error1,error2,error3,error4,error5,error6]):
						if modeled_protein == error:
							cnt=cnt+1
		
		pickle.dump(modeledPPIs, open(outPath, "wb"))	

		ppi_model_summary(modeledPPIs,ofile,outPath)
		
	else:

		modeledPPIs = pickle.load(open(outPath,"rb"))
		ppi_model_summary(modeledPPIs,ofile,outPath)
		
		
		
	return

def ppi_model_summary(modeledPPIs,ofile,outPath):

	# Check  number of proteins modeled for each PPI 
	# PPIs with no/1 modeled protein happened if pdb structure has too small coverage alignment & structural annotation transfer will fail
	count_noProtein,count_oneProtein = 0,0
	for key, ppi in list(modeledPPIs.items()):
		if len(ppi.proteins) == 0:
			count_noProtein=count_noProtein+1
			del modeledPPIs[key]
		elif len(ppi.proteins) == 1: # delete PPIs with only 1 modeled protein -for now- 
			count_oneProtein = count_oneProtein+1
			del modeledPPIs[key]
		elif len(ppi.proteins) == 2:
			continue
		else:
			print("error, too many proteins")
	
	print("PPIs with only 1 protein modeled (deleted):",count_oneProtein)
	print("PPIs with no modeled proteins (deleted): ",count_noProtein)
	print("PPIs with 2 proteins modeled :",len([ppi for ppi in modeledPPIs.values() if len(ppi.proteins)==2]))
	print("Total number of PPIs in the models:",len(modeledPPIs))
	print("Total number of residues in the models:",len([res for ppi in modeledPPIs.values() for protein in ppi.proteins for res in protein.residues]))
	print("Total number of interfacial residues in the models:",len([res for ppi in modeledPPIs.values() for protein in ppi.proteins for res in protein.residues if res.interfacial == True]))
	

	print("Total number of PPIs with at least 1 interface modeled (not size 0):", 
		len([ppi for ppi in modeledPPIs.values() if (ppi.proteins[0].interface[0].get_size()+ppi.proteins[1].interface[0].get_size())>0]))
	print("Total number of PPIs with at least 1 interface modeled and large enough for structural calculations (interface size >=10):", 
		len([ppi for ppi in modeledPPIs.values() if ( (ppi.proteins[0].interface[0].get_size()>=10) or (ppi.proteins[1].interface[0].get_size()>=10) ) ]))
	print("Total number of PPIs with  2 interfaces modeled (not size 0):", 
			len([ppi for ppi in modeledPPIs.values() if ( (ppi.proteins[0].interface[0].get_size()>0) and (ppi.proteins[1].interface[0].get_size()>0) ) ]))
	print("Total number of PPIs with 2 interfaces modeled and large enough for structural calculations (interface size >=10):", 
		len([ppi for ppi in modeledPPIs.values() if ( (ppi.proteins[0].interface[0].get_size()>=10) and (ppi.proteins[1].interface[0].get_size()>=10) ) ]))
	print("Total number of interfaces modeled (size >0)",
		len([intf for ppi in modeledPPIs.values() for protein in ppi.proteins for intf in protein.interface if (intf.get_size()>0) ])
		)
	print("Total number of interfaces modeled (size >=10)",
		len([intf for ppi in modeledPPIs.values() for protein in ppi.proteins for intf in protein.interface if (intf.get_size()>=10) ])
		)
	print("Total number of residues in interfaces modeled (size >0)",
		len([res for ppi in modeledPPIs.values() for protein in ppi.proteins for intf in protein.interface for res in intf.residues if (intf.get_size()>0) ])
		)
	print("Total number of residues in interfaces modeled (size >=10)",
		len([res for ppi in modeledPPIs.values() for protein in ppi.proteins for intf in protein.interface for res in intf.residues if (intf.get_size()>=10) ])
		)
	print("\n\n")
	# Re-pickle with the single prots and no prots deleted
	if count_noProtein != 0 or count_oneProtein !=0:
		pickle.dump(modeledPPIs, open(outPath, "wb"))	

	# intf_sizes = [intf.get_size() for ppi in modeledPPIs.values() for protein in ppi.proteins for intf in protein.interface]
	# plt.hist(intf_sizes)
	# plt.savefig('Interface sizes')
	 
	# Add check for if doesn;t exit
	with open(ofile,"w") as f:
		f.write('\t'.join(["ORF_A","ORF_B","NumReports","PDB", "chain_A","chain_B","e_value_A","e_value_B","coverage_A","coverage_B","interface_A_size","interface_B_size","interface_A_residues","interface_B_residues"])+"\n")

		for key, ppi in list(modeledPPIs.items()):
			orfA = ppi.orf1
			orfB = ppi.orf2
			prot_A = ppi.proteins[0] 
			prot_B = ppi.proteins[1] 
			if prot_A.orf != orfA:
				print("error1")
			if prot_B.orf != orfB:
				print("error2")
			NumReports = str(ppi.numberOfReports)
			PDB = prot_A.pdb[0] #=prot_B.pdb[0]
			chain_A = prot_A.pdb[1]
			chain_B = prot_B.pdb[1]
			e_value_A = prot_A.evalue
			e_value_B = prot_B.evalue
			coverage_A = prot_A.coverage
			coverage_B = prot_B.coverage
			if len(prot_A.interface) >1:
				print("too many interfaces A")
			if len(prot_B.interface) >1:
				print("too many interfaces B")
			if len(prot_A.interface) == 0:
				print("no interface A")
			if len(prot_B.interface) == 0:
				print("no interface B")
			interface_A_size = str(prot_A.interface[0].get_size())
			interface_B_size = str(prot_B.interface[0].get_size())
			interface_A_residues = ';'.join([res.id for interface in prot_A.interface for res in interface.residues])
			interface_B_residues = ';'.join([res.id for interface in prot_B.interface for res in interface.residues])

			f.write('\t'.join([orfA,orfB,NumReports,PDB,chain_A,chain_B,e_value_A,e_value_B,coverage_A,coverage_B,interface_A_size,interface_B_size,interface_A_residues,interface_B_residues])+"\n")

	return

	
def model_protein(orf,partner,species,pdb,pdbPartner,evalue,coverage,seq,pdbWholeProteinSeqs,orthologs,resolution,path_to_dssp,yeastAlignment):
	'''
	Models protein 'orf' involved in a PPI with the protein 'partner'
	Args:
		orf (str): 					Orf name of the protein to model
		partner (str): 				Orf name of the partner of the protein to model
		species (str):				'Scer' 
		pdb (str):					(PDB ID,PDB chain) for the structure best homology mapped to the protein to model
		pdbPartner (str):			(PDB ID,PDB chain) for the structure best homology mapped to the partner of the protein to model
		evalue (str):				e-value of the BLAST alignment of orf to PDB structure
		coverage (str):				coverage of orf by the PDB structure
		seq (str):					Original aa sequence of the protein in the species
		pdbWholeProteinSeqs (str):	Sequence (aa) of entire protein modeled in the PDB structure
		orthologs (dict):			Dictionary with the orf names of orthologs of 'orf' in closely related species. Format: {'Smik': 'ORFP:7528', 'Spar': 'ORFP:15866', 'Sbay': 'ORFP:16877'}
		resolution (float):			Resolution of the PDB structure best homology mapped to the protein to model
		path_to_dssp (str):			Path to the dssp executable
	Returns:
		Modeled protein with residues and interface modeled
		Error messages:
			error1: No sequence for orthologs in closely related species (needed for evolutionary rate calculations)
			error2: Coverage is too small on the PDB structure (cannot align sequence and structure)
			error3: Alignment of the protein sequence with the PDB structure and alignemnt of the protein sequence with orthologs in closely related species have different lengths
			error4: The protein has no modeled residues
			error5: The protein has no modeled interface
			error6: The protein has a modeled interface containing no residue
	'''

	# (1) Throw away proteins without yeast orthologs
	if (orf not in yeastAlignment): 
		return "error1"	

	# (2) Compute structural properties for residues in the protein
	pdbAligned = {}
	pdbID,chainID,partnerChain = pdb[0],pdb[1],pdbPartner[1]
	
	resInfo = calculate_structural_properties(pdbID,chainID,partnerChain,path_to_dssp)

	# Add the ID notation (identifier for the residue in the PDB file - 3LATOM notation. Ex: MET1:A) to the values in resInfo
	conversion = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
				 'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN', 
				 'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP', 
				 'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}
	for i in resInfo:
		aa = conversion[resInfo[int(i)]['aa']]
		idNotation = aa+str(i)+":"+chainID
		resInfo[int(i)].update({'id':idNotation})

	# (3) Align and transfers annotations back to PDB seq
	seqPDB = ''.join([r['aa'] for r in resInfo.values()])
	try:
		pdbAligned[orf] = align_yeast_to_pdb(seq,resInfo.values())

	except UserWarning as e: # If pdb structure has too small coverage alignment will fail
		print (e)
		print ('problem aligning', orf, 'and', pdbID, chainID)
		print ('seq1:', seq, '\nseq2:', seqPDB)
		return "error2"

	# (4) Create the protein model
	modeled_protein = Protein(orf,partner,species,pdb,evalue,coverage,pdbPartner,seq,seqPDB,pdbWholeProteinSeqs,orthologs,resolution)
	

	# (5) Create the interface model 
	modeled_protein.add_interface(Interface(orf,partner,pdb,pdbPartner))
	
	if len(list(pdbAligned[orf])) != len(list(yeastAlignment[orf])):
		raise UserWarning('Structure and yeast alignments differ')
		return "error3"

	# (6) Add the residues to our models
	for index, (algnYstCdn, algnPDBRes) in enumerate(zip(yeastAlignment[orf],pdbAligned[orf])):
		# Throw away gaps in any of the alignments
		if any(['-' in algnYstCdn[(i, 'nt')] for i in list(modeled_protein.orthologs.keys())]):
			continue
		if '-' in algnPDBRes:
			continue
		
		# Create the residue
		cdns = {k1: v for (k1, k2), v in algnYstCdn.items() if k2 == 'nt'}
		ID = algnPDBRes.get('id', '') 
		num_closely_related = len(list(modeled_protein.orthologs.keys()))
		
		res=Residue(ID,index,algnYstCdn[(species, 'aa')],cdns,num_closely_related, algnPDBRes['rsaMonomer'],algnPDBRes['rsaComplex'],algnPDBRes['drsa'])
		res.interRRC = algnPDBRes.get('interRRC', np.nan)
		res.dCenter = algnPDBRes.get('dCenter', np.nan)
		res.dEdges = algnPDBRes.get('dEdges', np.nan)
		if(len(set(res.codons.values()))==1):
			res.substitution = 0
		else:
			res.substitution = 1

		# Add residue to the protein model (will also be added to the interface model if interfacial)
		modeled_protein.add_residue(res)


	if len(modeled_protein.residues) == 0:
		return "error4"
	# if len(modeled_protein.interface) ==0:
	# 	return "error5"
	# if modeled_protein.interface[0].get_size() ==0:
	# 	return "error6"
	
	return modeled_protein

def align_yeast_to_pdb(sequence, pdbProtein):
	"""Align RSA/deltaRSA of PDB complex to yeast protein.

	Args:
		spomAASeq (str): aa sequence of s. pom protein
		pdbProtein (list(dict)): residue-level data of homolog in the PDB,
								 including AA sequence, solvent accessibility
								 etc.

	Returns:
		list(dict): alignment of RSA and other values of PDB chain to yeast protein

	"""
	seqPDB = ''.join([r['aa'] for r in pdbProtein])
	iAlnd = align(sequence, seqPDB, [str(i) for i in range(len(seqPDB))])

	return [i if '-' in i else list(pdbProtein)[int(i)] for i in iAlnd]


def calculate_structural_properties(pdbID,chainID,partnerChain,path_to_dssp):
	'''
	Compute structural properties for residues in the protein modeled as chainID, 
	which interacts with the protein modeled as partnerChain in the structure pdbID
	- For all residues in the protein we compute 
		rsaMonomer (float)			: 	The solvent accessibility of the residue when the protein is in monomer state
		rsaComplex (float)			: 	The solvent accessibility of the residue when the protein is co-complexed with partnerChain
		drsa (float)      			:	The change in solvent accessibility upon complex formation (rsaMonomer-rsaComplex)
		interface (bool)  			:	Whether or not the residue is interfacial (has a change in solvent accessibility upon complex formation i.e. drsa !=0)
	- For interfacial residues in the protein we compute
		distance_from_center (float):	The distance between the residue and the geometric center of the interface
		distance_from_edges (float) :	The distance between the residue and the edges of the interface (90th percentile of rsaComplex)
		contacts (float)			:	The number of contacts that the residue makes with residues in the partner protein
	Args:
		pdbID (str):		ID for the PDB structure that best models the PPI			
		chainID (str):		Chain ID for the protein of interest to the calculations
		partnerChain (str):	Chain ID for the partner protein in the PPI ('' if we are modeling a single protein)
		path_to_dssp (str):	Path to the dssp executable
	Returns:
		resInfo (OrderedDict): Format: { residueChainID : {'aa':'M','rsaMonomer': float, 'rsaComplex': float, 'drsa': float, 'interface': bool,'distance_from_center': float (only if 'interface'=True), 'distance_from_edges': float (only if 'interface'=True), 'contacts': int (only if 'interface'=True)}}
	'''

	minNonZeroDeltaRSA = 0.001
	inDir = '../data/processed/assemblies/'
	struct = load_structure(inDir + pdbID + '.pdb1')
	
	# (1) Compute solvent accessibilities (rsaMonomer,rsaComplex,drsa) for all residues in the protein
	sasa = sasa_scan(struct, chainID, chainID + partnerChain,path_to_dssp, True)[0]

	resInfo = OrderedDict([(k, round_delta_RSA(v, minNonZeroDeltaRSA))
					   for k, v in sasa.items()])
	
	# (2) Label interfacial residues
	for i, res in enumerate(resInfo.values()):
		if res['drsa'] > 0.0:
			res['interface'] = True
		else:
			res['interface'] = False

	# (3) Compute properties for interfacial residues
	parser = PDBParser(PERMISSIVE=1, QUIET=1)
	ipath = '../data/processed/assemblies/' + pdbID + '.pdb1' 
	structure = parser.get_structure('name', ipath)[0]
	protein = structure[chainID]
	
	intf = []				# List of interfacial residues
	intf_complexRSA = [] 	# List of all values of complex RSA, to plot/ get the 90th percentile of the distribution

	# Get interfacial residues list (intf) and rsaComplex for all interfacial residues (intf_complexRSA)
	for i, residue in enumerate(protein):
		resID = residue.id[1]
		# Skip small number of residues not parsed by dssp
		if resID not in resInfo:
			continue  
		if (resInfo[resID]['interface']): 
			intf.append(residue)
			intf_complexRSA.append((residue,resInfo[resID]['rsaComplex']))

	# Compute structural properties of interfacial residues
	if (len(intf)>=10): # Only for interfaces that are large enough
		partner = structure[partnerChain]
		# Define the outer boundaries of the interface (edges) (90th percetile largest complex RSA values)
		rsacomplexes = [x[1] for x in intf_complexRSA]
		percentile = np.percentile(rsacomplexes, 90)
		outerBoundaries = [x[0] for x in intf_complexRSA if x[1] >= percentile]
		
		# Compute the geometrical mean of the interface
		intfCenter = np.mean([r['CA'].coord for r in intf], axis=0) 	# x,y,z coords of geometric center of the interface
	
		# For interfacial residues in the protein, 
		for i, residue in enumerate(protein):
			resID = residue.id[1]
			if resID not in resInfo: 				# Skip small number of residues not parsed by dssp
				continue
			if not resInfo[resID]['interface']: 	# Only for interfacial residues
				continue

			cutoff = 10.0  # Angstroms cutoff to define interface contacts
			dists = list(map(lambda r: residue['CA'] - r['CA'], partner))		# Number of residue contacts with the partner protein
			d = np.sqrt(np.sum((residue['CA'].coord - intfCenter) ** 2))		# Distance of residue to geometric center
			e = min(map(lambda r: residue['CA'] - r['CA'], outerBoundaries))	# Distance of residue to outer boundary
			
			# Assign distance to geometric center
			resInfo[resID]['dCenter'] = d
			# Assign distance to geometric center
			resInfo[resID]['dEdges'] = e
			# Assign number of residue contacts
			resInfo[resID]['interRRC'] = 0 # initialize to 0
			resInfo[resID]['interRRC'] += np.sum([d < cutoff for d in dists])

	return resInfo

def round_delta_RSA(residue, threshold):
	"""Round down values of delta RSA below a threshold to zero.
	Args:
		residue (dict)
		threshold (float): minimum change in RSA considered non-zero

	Returns:
		dict
	"""
	if residue['drsa'] >= threshold or residue['drsa'] == 0.:
		return residue
	else:
		return {'aa': residue['aa'], 'rsaMonomer': residue['rsaMonomer'],'rsaComplex':residue['rsaComplex'] ,'drsa': 0.0}


def interface_geometry(pdbID, chainID, partnerChainID, resInfo):
	"""Quantities of interest for interfacial residues.
	Quantities:
	- rsa: 	relative solvent accessibility in uncomplexed state
	- dRSA: change in rsa upon complex formation
	Both computed in calc_sasa. Check if saved under correct model portion
	- Residue contacts: Number of residues in the partner protein within 10 Angstroms from our residue
	- Distance from center of the interface: distance from the geometric center of the interface. 
	- Distance from edges of the interface: defined as 90th percentile of largest complex RSA
	Args:
		pdbID (str)
		chainID (str)
		partnerChainIDs (str)
		resInfo (OrderdDict(int: dict)): PDB residue ID -> solvent accessibility values of residue.
										Used to define interface residues. 

	"""
	
	parser = PDBParser(PERMISSIVE=1, QUIET=1)
	ipath = '../data/processed/assemblies/' + pdbID + '.pdb1' 
	structure = parser.get_structure('name', ipath)[0]
	protein = structure[chainID]
	partner = structure[partnerChainID]
	intf = []
	intf_complexRSA = [] # List of all values of complex RSA, just used to plot/ get the 90th percentile of the distribution

	# Get  interfacial residues
	for i, residue in enumerate(protein):
		resID = residue.id[1]
		# Skip small number of residues not parsed by dssp
		if resID not in resInfo:
			continue  

		# Find the interfacial residues (dRSA > 0)
		if ('interface' in resInfo[resID]): 
			intf.append(residue)
			intf_complexRSA.append((residue,resInfo[resID]['rsaComplex']))
		
	
	# Only for interfaces that are big enough
	if (len(intf)<10):
		print("Too small, interface length= ",len(intf))
	else:	
		# Define the outer boundaries of the interface (edges) (95th percetile largest complex RSA values)
		rsacomplexes = [x[1] for x in intf_complexRSA]
		percentile = np.percentile(rsacomplexes, 95)
		outerBoundaries = [x[0] for x in intf_complexRSA if x[1] >= percentile]
		outerBoundariesDist = []
		# Get the average distance to the outer boundary 
		for residue in intf:
			dists = map(lambda r: residue['CA'] - r['CA'], outerBoundaries) # Distances of the residue to all residues in the outer boundary
			outerBoundariesDist.append(min(dists)) 							# Smallest distance to a residue in the outer boundary
		
		# Compute the geometrical mean of the interface
		intfCenter = np.mean([r['CA'].coord for r in intf], axis=0) 	# x,y,z coords of geometric center of the interface
		intfDist = [np.sqrt(np.sum((r['CA'].coord - intfCenter) ** 2)) 	# list of average distance of all interfacial residues to the geometric center of the interface
					for r in intf]
	
		# For all residues in the protein, 
		for i, residue in enumerate(protein):
			resID = residue.id[1]
			if resID not in resInfo: # Skip small number of residues not parsed by dssp
				continue
			if 'interface' not in resInfo[resID]: # Only for interfacial residues
				continue

			cutoff = 10.0  # Angstroms cutoff to define interface contacts
			dists = list(map(lambda r: residue['CA'] - r['CA'], partner))		# Number of residue contacts with the partner protein
			d = np.sqrt(np.sum((residue['CA'].coord - intfCenter) ** 2))		# Distance of residue to geometric center
			e = min(map(lambda r: residue['CA'] - r['CA'], outerBoundaries))	# Distance of residue to outer boundary
			
			# Assign numerical distance to geometric center
			resInfo[resID]['distance_from_center'] = d
			# Assign numerical distance to geometric center
			resInfo[resID]['distance_from_edges'] = e
			# Assign number of residue contacts
			resInfo[resID]['contacts'] = 0 # initialize to 0
			resInfo[resID]['contacts'] += np.sum([d < cutoff for d in dists])
			

	return

def load_aligned_yeast_codons(codonAlignmentPath):
	"""
	Load codon alignment file between all orfs in a species and closely related species (generated in process_data.py) into a dictionary
	Format: yeastAlignment = { orf : [{codon1},{codon2},...]} with orf each orf in species and {codon} a dictionary with the sequence at the codon position in each of the closely related species
			example: yeastAlignment["SPBC106.05c"] = [{('Spom', 'aa'): 'M',('Spom', 'nt'): 'ATG',...},{},...]
	Args:
		codonAlignmentPath (str): path to the codon alignment file
	"""
	yeastAlignment = {}
	with open(codonAlignmentPath, 'r') as f:
		for l in f:
			scerORFID = l.split()[0]
			if scerORFID not in yeastAlignment:
				yeastAlignment[scerORFID] = {}
			strainName, seqType = l.split()[1].split('.')
			yeastAlignment[scerORFID][(strainName, seqType)] = l.split()[2:]

	# Switch to a more convenient format
	algn = {}
	for orf in yeastAlignment:
		# Verifying that the lengths are consistent accross closely related species
		lengths =  set([len(yeastAlignment[orf][key]) for key in yeastAlignment[orf].keys()])
		numAA = list(lengths)[0]
		if(len(lengths) != 1):
			raise UserWarning('Inconsistent sequence lengths in:' + codonAlignmentPath)
		# Format change
		algn[orf] = list(map(lambda i: {k: yeastAlignment[orf][k][i] for k in yeastAlignment[orf].keys()}, range(numAA)))

	return algn


def main():
	externalData = '../data/external'
	processedData = '../data/processed'
	blastToPDBData = os.path.join(externalData, 'blastToPDB')

	# Running the models 3 times:
	
	scer_closely_related = ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']
	scer_additional_closely_related = ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Ctro', 'Tpha','Cfab','Ndai','Zrou']
	scer_additional_distantly_related = ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Fgra','Ncra','Anid','Spom','Nirr','Plac','Pjir','Pmur']
	
	alignmentFile_closely_related = os.path.join(externalData,'codon_alignment.txt')
	alignmentFile_additional_closely_related = os.path.join(externalData,'codon_alignment_AdditionalCloseSpecies.txt')
	alignmentFile_additional_distantly_related = os.path.join(externalData,'codon_alignment_AdditionalDistantSpecies.txt')

	pickled_path_closely_related =os.path.join(processedData,'modeled_ppis.pkl')
	pickled_path_additional_closely_related =os.path.join(processedData,'modeled_ppis_additional_closely_related.pkl')
	pickled_path_additional_distantly_related =os.path.join(processedData,'modeled_ppis_additional_distantly_related.pkl')

	summary_path_closely_related =os.path.join(processedData,'model_summary.txt')
	summary_path_additional_closely_related =os.path.join(processedData,'model_summary_additional_closely_related.txt')
	summary_path_additional_distantly_related =os.path.join(processedData,'model_summary_additional_distantly_related.txt')

	for related_species,alignmentFile,ofile1,ofile2 in zip( [scer_closely_related,scer_additional_closely_related,scer_additional_distantly_related],
															[alignmentFile_closely_related,alignmentFile_additional_closely_related,alignmentFile_additional_distantly_related],
															[pickled_path_closely_related,pickled_path_additional_closely_related,pickled_path_additional_distantly_related],
															[summary_path_closely_related,summary_path_additional_closely_related,summary_path_additional_distantly_related]):
	

		print ("----------- Loading input -----------")
		(yeastPPIs,orfs, matchStats, pdbFullSeqs, aaSeq,structs, pdbLengths, resolutions, orthologs) = load_inputs(externalData,processedData,blastToPDBData,related_species)
		
		# # For testing purposes, build models for a subset of the structures (small_structs)
		# structs = dict(random.sample(structs.items(), 50))
		
		print ("----------- Selecting best PDB matches -----------")
		bestMatch,PPIs_to_model,orfs_in_PPIs_to_model = select_best_pdb_matches(structs,
											matchStats,
											pdbLengths,
											resolutions,
											yeastPPIs)

		## Testing which species have orthologs info in related species for future dN/dS calculations
		# orfs_with_all_related_species = [orf for (orf,spec) in orthologs.items() if len(spec) == len(related_species)]
		# print(len(orfs_with_all_related_species), len(set(orfs_with_all_related_species)))
		# ppis_with_all_related_species = [ppi for ppi in PPIs_to_model if ppi.ORF_A in orfs_with_all_related_species and ppi.ORF_B in orfs_with_all_related_species]
		# print(len(ppis_with_all_related_species),"/",len(PPIs_to_model),"have MSA to perform evolutionary rate calculations")
		
		# Download DSSP
		path_to_dssp = download_dssp()
		path_to_dssp ="../data/dssp-2.0.4/mkdssp"
		
		print ("----------- Aligning with structures -----------")
		align_with_structures_multiple(PPIs_to_model,alignmentFile,yeastPPIs,structs,bestMatch,
										aaSeq, resolutions, pdbFullSeqs, orthologs, path_to_dssp,
										ofile1,ofile2)
		

		
	
if __name__ == '__main__':
	main()


