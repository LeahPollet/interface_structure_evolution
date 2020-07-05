#!/usr/bin/env python

"""
Bin codons according to structural features of the corresponding residues and into random bin.
Preliminary step before computing dN/dS

File outputs:
	- many fasta files for each binning of the codons
	- many files containing average RSA values in different bins
	- csv files of protein and interface properties

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os
from collections import Counter
import pickle
import numpy as np
import pandas as pd
from amino_acid import AminoAcid
from modeled_proteins import Protein, Residue
from pathlib import Path

import math
import matplotlib.pyplot as plt
import sys


MODELS_PATH = '../../data/processed/modeled_proteins.pkl'
OUT_DIR = '../../data/evolutionary_rates/'+os.path.basename(MODELS_PATH)
FASTA_DIR = os.path.join(OUT_DIR, 'binned_codons')
RSA_DIR = os.path.join(OUT_DIR, 'RSA_distributions')


def write_readme_file(subdir,binEdges,readmeTxT):
	""" Write a readme file detailing the bin edges used """

	# Format the bin edges used as a readable string
	listEdges = [(dn, up) for i,(dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:]))]
	
	if listEdges[0][0] == listEdges[0][1]: 	# First bin has lower bound included
		listEdges_string = "["
	else:									# Lower bound of first bin not included	
		listEdges_string = "]"

	for i in range(len(listEdges)):
		listEdges_string = listEdges_string + str(listEdges[i][0])+","+str(listEdges[i][1])+"]\n]"
	listEdges_string = listEdges_string[:-1] # Remove the last open backet

	# Write the file
	ofile = os.path.join(subdir,'README.txt')
	with open(ofile, 'w') as f:
		f.write(readmeTxT)
		f.write("\n"+str(len(listEdges))+" bins used: \n")
		f.write(listEdges_string)

	return
def save_binned_aligned_codons_random(buckets, binEdges,varName, subdir, names, readmeTxT, dir=FASTA_DIR):
	"""Save binned aligned codons as fasta files in a  directory."""

	# Create the directory to store fasta files
	if not os.path.exists(subdir):
		os.makedirs(subdir)
		# Add a README.txt file to describe the binning
		if not readmeTxT==None:
			write_readme_file(subdir,binEdges,readmeTxT)
			
	# Write the fasta file for each bin
	for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:])):
		if len(buckets[i]) == 0:
			raise UserWarning('Empty bin: ',varName + '_' + str(dn) + '_to_' + str(up)+".fasta")
			# Should add a min number of residue needed for dnds calculations
			return

		if names is None:
			name = varName + '_' + str(dn) + '_to_' + str(up)+".fasta"
		else:
			name = str(names[i])+".fasta"

		if not os.path.exists(os.path.join(subdir,name)):
			with open(os.path.join(subdir,name), 'w') as f:
				for strain in ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
					f.write('>' + strain + '\n')
					f.write(''.join([r.codons[strain] for r in buckets[i]]) + '\n')
	return

def save_binned_aligned_codons(buckets, binEdges,varName, subdir, names, readmeTxT, dir=FASTA_DIR):
	"""Save binned aligned codons as fasta files in a  directory."""

	# Create the directory to store fasta files
	if not os.path.exists(subdir):
		os.makedirs(subdir)
		# Add a README.txt file to describe the binning
		if not readmeTxT==None:
			write_readme_file(subdir,binEdges,readmeTxT)
			
	# Write the fasta file for each bin
	for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:])):
		if len(buckets[i]) == 0:
			raise UserWarning('Empty bin: ',varName + '_' + str(dn) + '_to_' + str(up)+".fasta")
			# Should add a min number of residue needed for dnds calculations
			return

		if names is None:
			name = varName + '_' + str(dn) + '_to_' + str(up)+".fasta"
		else:
			name = str(names[i])+".fasta"

		if not os.path.exists(os.path.join(subdir,name)):
			with open(os.path.join(subdir,name), 'w') as f:
				for strain in ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
					f.write('>' + strain + '\n')
					f.write(''.join([r.codons[strain] for r in buckets[i]]) + '\n')
	return

def save_average_RSA_for_bins(buckets,binEdges,saveDir,varName):
	""" Saves a file with the average monomer and complex RSA value for each bin in binEdges"""

	rsaPath = os.path.join(saveDir,'rsa_by_' + varName + '.csv')
	if not os.path.exists(rsaPath):
		with open(rsaPath, 'w') as f:
			f.write('bin_low,bin_high,number_residues,''mean_rsa_monomer,mean_rsa_complex\n')
	
			for res, dn, up in zip(buckets, binEdges[:-1], binEdges[1:]):
				nResidues = len(res)
				if nResidues > 0:
					meanRSAMonomer = np.mean([r.rsa for r in res])
					meanRSAComplex = np.mean([r.rsa - r.dRsa for r in res])
				else:
					meanRSAMonomer = np.nan
					meanRSAComplex = np.nan
				
				f.write(','.join(map(str, [dn,
										   up,
										   nResidues,
										   meanRSAMonomer,
										   meanRSAComplex])))
				f.write('\n')

	return

def barPlots(buckets,binEdges,varName,saveDir):
	'''Plot histograms for the number of codons in each bin'''

	counts = [len(bucket) for bucket in buckets]
	y_pos = np.arange(len(counts))
	labels = [str((round(dn,2),round(up,2))) for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:]))]#[1:]#Don't plot the zero bin

	plt.bar(y_pos, counts, align='center', alpha=0.5)
	plt.xticks(y_pos, labels,rotation=45)
	plt.ylabel('Number of residues in bin')
	plt.xlabel(str(varName)+" bin")
	plt.title("Binning residues by "+str(varName))
	plt.tight_layout()
	plt.savefig(os.path.join(saveDir,str(varName)+'_binning_graph.png'),dpi=300)
	plt.show()

	# If we need the plot without the [0,0] bin:
	if binEdges[0] == binEdges[1]:
		counts = [len(bucket) for bucket in buckets][1:] #Don't plot the zero bin
		y_pos = np.arange(len(counts))
		labels = [str((round(dn,2),round(up,2))) for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:]))][1:]#Don't plot the zero bin

		plt.bar(y_pos, counts, align='center', alpha=0.5)
		plt.xticks(y_pos, labels,rotation=45)
		plt.ylabel('Number of residues in bin')
		plt.xlabel(str(varName)+" bin")
		plt.title("Binning residues by "+str(varName))
		plt.tight_layout()
		plt.savefig(os.path.join(saveDir,str(varName)+'_binning_graph_noFirstBin.png'),dpi=300)
		plt.show()

	return

def split(a, n):
	
	k, m = divmod(len(a), n)
	return list(a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n))


def bin_codons_random(proteins, binNum,saveDir):
	"""Bin codons into random groups and save them.

	Args:
		proteins (dict(proteins))
		binNum (int): number of bins wanted
	"""
	res_interfacial= [r for p in proteins.values() for r in p.residues if r.dRsa != 0]
	binned_residues = split(res_interfacial,binNum)
	
	properties = pd.DataFrame(columns=['bin_number','num_of_residues','RSA_monomer', 'RSA_complex', 'deltaRSA',
							   'residueContacts', 'distance_from_center',
							   'distance_from_edges'])

	# Write the summary text file for all bins
	for idx,bin in enumerate(binned_residues):
		row={}
		row['bin_number']			= idx
		row['num_of_residues']		= len(bin)
		row['RSA_monomer']			= float((np.nanmean([float(r.rsa) for r in bin])))
		row['RSA_complex']			= float((np.nanmean([float(r.rsa - r.dRsa) for r in bin])))
		row['deltaRSA']				= float((np.nanmean([float(r.dRsa) for r in bin])))
		row['residueContacts']		= float((np.nanmean([float(r.contacts) for r in bin])))
		row['distance_from_center']	= float((np.nanmean([float(r.distance_from_center) for r in bin])))
		row['distance_from_edges']	= float((np.nanmean([float(r.distance_from_edges) for r in bin])))
		properties = properties.append(row, ignore_index=True)

	properties = properties.astype({"bin_number": int, "num_of_residues": int})
	readmeTxT="All interfacial residues split into "+str(binNum)+"equal size bins (~"+str(len(binned_residues[0]))+" residues per bin) \n"+str(properties.to_string())
	if not readmeTxT==None:
		ofile = os.path.join(saveDir,'README.txt')
		with open(ofile, 'w') as f:
			f.write(readmeTxT)

	# Write the fasta file for each bin
	for bin_num in range(binNum):
		if len(binned_residues[bin_num]) == 0:
			raise UserWarning('Empty bin: ',varName + '_' + str(dn) + '_to_' + str(up)+".fasta")
			return
		fasta_name = os.path.join(saveDir,str(bin_num)+".fasta")	
		if not os.path.exists(fasta_name):
			with open(fasta_name, 'w') as f:
				for strain in ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
					f.write('>' + strain + '\n')
					f.write(''.join([r.codons[strain] for r in binned_residues[bin_num]]) + '\n')
		
	return

def bin_codons(residueVar, binEdges, proteins, varName, saveDir, 
			   fNames=None, saveRSA=False,rsaDir=RSA_DIR,readmeTxT=None,plot=False):
	"""Bin codons by a single continuous residue-level variable and save them.

	Saves the bins for each of the protein and residue level categories.

	Args:
		residueVar (function): call to access value of residue variable.
		binEdges (list(float)): edges of the bins.
		proteins (dict(proteins))
		varName (str): name of variable that is binned on.
		saveDir(str): path to the directory where to save the binning
		fNames (str): overwrite file name.
		saveRSA (bool): also save the average RSA of each bin
		readmeTxT (str): text to include as a README.txt file
		plot (bool): whether or not to plot the number of residues in each bin as a histogram

	"""
	# Bin the codons
	if len(binEdges) < 2:
		raise ValueError('Expected at least 2 bin edges')
	
	buckets = {}
	
	buckets = [[] for __ in range(len(binEdges) - 1)]
	
	for p in proteins.values():
		for r in p.residues:
			val = residueVar(r)
			if val == np.nan or type(val)== str:
				continue
			if val > binEdges[-1]:
				raise UserWarning('Overflow from bin.')
			for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:])):
				if (val <= up and val > dn) or (dn == up and val == dn):
					buckets[i].append(r)
	# Save the binning
	save_binned_aligned_codons(buckets,binEdges,varName, saveDir, fNames, readmeTxT)
	
	# Save RSA
	if saveRSA:
		save_average_RSA_for_bins(buckets,binEdges,rsaDir,varName)	

	# Visualize the binning
	if plot:
		barPlots(buckets,binEdges,varName,saveDir)
		
	return


def bin_all(proteins):
	"""Set the bins that residues are grouped into.

	Bins by the protein-level and residue-level categories defined by the tags
	and by RSA, dRSA, CAI.

	Args:
		proteins ({str: Protein}): structurally modeled proteins.

	"""

# 1) Define the edges of the bins for all categories (Value val is included in bin ]dn,up] if (val <= up and val > dn) or (dn == up and val == dn))
	binEdgesRSA_interface = np.arange(0.00, 1.01, 0.05)			# [0,0] bin is empty for interfaces
	binEdgesRSA = np.insert(binEdgesRSA_interface, 0, 0.00) 	# adding a [0,0] bin for all residues
	binEdgesDeltaRSA= np.arange(0.00, 1.01, 0.1)
	binEdgesDeltaRSA= np.insert(binEdgesDeltaRSA, 0, 0.00) 		# adding a [0,0] bin
	binEdgesContacts = [0, 0, 2, 4, 6, 8, 10, 12, 14, 30]
	binEdgesDistance_from_center = [0,10,15,20,25,30,40,50,65,93]
	binEdgesDistance_from_edges = [0,10,15,20,25,30,35,45,65,100]
	
	proteinsWithInterface = {k: v for k, v in proteins.items() if len(v.interfaces) > 0}
	print(len(proteinsWithInterface),"/" ,len(proteins),"proteins do not have a modeled interface. Discarded from the interface involvement analysis")

# 2) For all proteins 
	if not os.path.exists(FASTA_DIR):
		os.makedirs(FASTA_DIR)
	if not os.path.exists(RSA_DIR):
		os.makedirs(RSA_DIR)

# - RSA (monomer) 
	# All residues
	bin_codons(lambda r: r.rsa, binEdgesRSA, proteins, 'rsa_monomer_all',  
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_all"),saveRSA = True, 
				rsaDir = RSA_DIR,readmeTxT = "Monomer RSA, for all residues in all proteins.", plot=False) 
	# Only interfacial residues (drsa != 0)
	bin_codons(lambda r: r.rsa if (r.dRsa > 0.0) else (np.nan), binEdgesRSA_interface, proteins, 'rsa_monomer_interfacial',  
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_interfacial"),saveRSA = True, 
				rsaDir = RSA_DIR,readmeTxT = "Monomer RSA, for all interfacial residues in all proteins.", plot=False) 

	# Only non-interfacial residues (drsa == 0)
	bin_codons(lambda r: r.rsa if (r.dRsa == 0.0) else (np.nan), binEdgesRSA, proteins, 'rsa_monomer_non_interfacial',  
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_non_interfacial"),saveRSA = True, 
				rsaDir = RSA_DIR,readmeTxT = "Monomer RSA, for all non-interfacial residues in all proteins.", plot=False) 
	
# - RSA (complex) 
	# All residues
	bin_codons(lambda r: r.rsa - r.dRsa, binEdgesRSA, proteins, 'rsa_complex_all', 
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_all"),saveRSA = True, 
				rsaDir = RSA_DIR, readmeTxT="Complex RSA, for all residues in all proteins.",plot=False) 
	# Only interfacial residues (drsa != 0)
	bin_codons(lambda r: r.rsa - r.dRsa if (r.dRsa > 0.0) else (np.nan), binEdgesRSA, proteins, 'rsa_complex_interfacial', 
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_interfacial"),saveRSA = True, 
				rsaDir = RSA_DIR, readmeTxT="Complex RSA, for all interfacial residues in all proteins.",plot=False) 
	# # Only non-interfacial residues (drsa == 0)
	bin_codons(lambda r: r.rsa - r.dRsa if (r.dRsa == 0.0) else (np.nan), binEdgesRSA, proteins, 'rsa_complex_non_interfacial', 
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_non_interfacial"),saveRSA = True, 
				rsaDir = RSA_DIR, readmeTxT="Complex RSA, for all non-interfacial residues in all proteins.",plot=False) 
	
		
# - Measures of interface involvement 	
	# - Delta RSA 
	bin_codons(lambda r: r.dRsa, binEdgesDeltaRSA, proteins,'drsa', 
				saveDir= os.path.join(FASTA_DIR,"drsa"), saveRSA=True, 
				rsaDir = RSA_DIR, readmeTxT='deltaRSA (Difference in solvent accessible surface between complex and monomer state?) for all residues in all proteins. [0,0] bin represents non-interfacial residues.', 
				plot=False) 
	# - Residue contacts
	bin_codons(lambda r: r.contacts, binEdgesContacts, proteins,'residueContacts',
				saveDir= os.path.join(FASTA_DIR,"residueContacts"),saveRSA=True, 
				rsaDir = RSA_DIR, readmeTxT='Residue contacts (Number of contacts with residues in partner protein) for all residues in all proteins. We only have residue contact values for interfacial residues in proteins with a modeled interface', 
				plot=False) 
	# - Distance from center
	bin_codons(lambda r: r.distance_from_center, binEdgesDistance_from_center, 
				proteins,'distance_from_center',saveDir= os.path.join(FASTA_DIR,"distance_from_center"), 
				saveRSA=True, rsaDir = RSA_DIR, readmeTxT='Binning by distance to geometric center of the interface for all residues in all proteins. We only have distance from center values for interfacial residues in proteins with a modeled interface', 
				plot=False) 
	# - distance from edges
	bin_codons(lambda r: r.distance_from_edges, binEdgesDistance_from_edges, 
				proteins,'distance_from_edges',saveDir= os.path.join(FASTA_DIR,"distance_from_edges"), 
				saveRSA=True, rsaDir = RSA_DIR, readmeTxT='Binning by distance to edges (most exposed surface in complex state) of the interface for all residues in all proteins. We only have distance from edges values for interfacial residues in proteins with a modeled interface', 
				plot=False) 

	# Create 'dir_list.txt' file, listing the directories for which to calculate dnds
	# Modify manually if calculations are not needed for everything
	# Modify to parallelize
	ofile = os.path.join(FASTA_DIR,'dir_list.txt') 
	subfolders = [root for root,dirs,files in os.walk(FASTA_DIR) if not dirs] # Only innermost directories ("leafs")

	with open(ofile, 'w') as f:
		for path in subfolders:
			f.write(path)
			f.write("\n")
	
	return

def bin_random(proteins):
	# Create the directorty
	fasta_path100 = os.path.join(FASTA_DIR,"100_random_bins")
	if not os.path.exists(fasta_path100):
		os.makedirs(fasta_path100)
	
	res = [r for p in proteins.values() for r in p.residues ]
	res_interfacial= [r for p in proteins.values() for r in p.residues if r.dRsa != 0]
	print(len(res_interfacial),"/" ,len(res),"residues are interfacial")
	
	# 100 bins of ~ 435 residues
	bin_codons_random(proteins, 100,fasta_path100)

	return	

def main():

	if not os.path.exists(FASTA_DIR):
	    os.makedirs(FASTA_DIR)
	
	print ("loading protein models")
	with open(MODELS_PATH, 'rb') as f:
		proteins = pickle.load(f)
	print ("protein models loaded")
	#Bin by structural properties
	bin_all(proteins)
	#Random bins of residues for possible multiple linear regression analysis
	bin_random(proteins)

if __name__ == '__main__':
	main()
