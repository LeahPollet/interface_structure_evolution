#!/usr/bin/env python

"""
Bin codons according to structural features of the corresponding residues.

File outputs:
	- Multiple fasta files for each binning of the codons
	- Multiple files containing average RSA values in different bins
	- csv files of protein and interface properties
"""
import os
from collections import Counter
import pickle
import numpy as np
import pandas as pd
from modeled_PPIs import PPI,Protein,Interface,Residue
from pathlib import Path
import warnings
from tqdm import tqdm


pd.set_option("display.max_rows", 20, "display.max_columns", 20)

def save_averages_for_bins(df,varName,saveDir):
	""" Saves a file with the average monomer and complex RSA value for each bin in binEdges"""
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)

	rsaPath = os.path.join(saveDir,'averages_for_' + varName + '.csv')
	
	if not os.path.exists(rsaPath):
		with open(rsaPath,"w") as f:
			df.to_csv(rsaPath,index= False,na_rep='NaN')
		
	return

def bin_codons_two_bins(bin1,name1, bin2, name2, saveDir, strains,RSA_DIR,saveRSA = False,plot=False):
	''' Bin codons into 2 bins, writes a README.txt file with details of the binning, creates a fasta file for each bin
	bin1,name1 	: codons in the first bin and name of first bin
	bin2,name2 	: codons in the second bin and name of second bin
	strains 	: closely related species to use for the dN/dS calculations
	saveDir 	: directory where to save the fasta a readme files
	'''
	bins=[(name1,bin1),(name2,bin2)]
	
	properties = pd.DataFrame(columns=['bin_name','num_of_residues','RSA_monomer', 'RSA_complex', 'deltaRSA',
							   'interRRC', 'dCenter','dEdges','consurfDB_score','consurfDB_std','consurfDB_resNum','consurfRate4Site_score','consurfRate4Site_std','consurfRate4Site_resNum'])

	# Write the summary text file for all bins
	for name,bin in bins:

		row={}
		row['bin_name']			= name
		row['num_of_residues']		= len(bin)
		row['RSA_monomer']			= float((np.nanmean([float(r.rsaMonomer) for r in bin])))
		row['RSA_complex']			= float((np.nanmean([float(r.rsaComplex) for r in bin])))
		row['deltaRSA']				= float((np.nanmean([float(r.dRsa) for r in bin])))
		row['interRRC']				= float((np.nanmean([float(r.interRRC) for r in bin])))
		row['dCenter']				= float((np.nanmean([float(r.dCenter) for r in bin])))
		row['dEdges']				= float((np.nanmean([float(r.dEdges) for r in bin])))
		row['consurfDB_score']		= float((np.nanmean([float(r.consurfDB_score) for r in bin])))
		row['consurfDB_std']		= float((np.nanstd([float(r.consurfDB_score) for r in bin])))
		row['consurfDB_resNum']		= len([float(r.consurfDB_score) for r in bin if not np.isnan(r.consurfDB_score)])
		row['consurfRate4Site_score']	= float((np.nanmean([float(r.consurf_score) for r in bin])))
		row['consurfRate4Site_std']		= float((np.nanstd([float(r.consurf_score) for r in bin])))
		row['consurfRate4Site_resNum']	= len([float(r.consurf_score) for r in bin if not np.isnan(r.consurf_score)])
		properties = properties.append(row, ignore_index=True)

	properties = properties.astype({"num_of_residues": int})
	properties = properties.astype({"consurfDB_resNum": int})
	properties = properties.astype({"consurfRate4Site_resNum": int})

	df= properties[['consurfDB_score','consurfDB_std','consurfDB_resNum','consurfRate4Site_score','consurfRate4Site_std','consurfRate4Site_resNum']]
	
	
	readmeTxT="#All residues in modeled S.cerevisiae PPIs split into two bins: \n#"+name1 +" bin (" +str(len(bin1))+" residues) \n#"+name2+" bin (" +str(len(bin2))+" residues)\n#\nSummary:\n"+str(properties.to_string())
	if not readmeTxT==None:
		ofile = os.path.join(saveDir,'README.txt')
		if not os.path.exists(ofile):
			with open(ofile, 'w') as f:
				f.write(readmeTxT)
	
	# Write the fasta file for each bin
	for name,bin in bins:
		if len(bin) == 0:
			raise UserWarning('Empty bin: ',varName + '_' + str(dn) + '_to_' + str(up)+".fasta")
			return
		fasta_name = os.path.join(saveDir,str(name)+".fasta")	
		if not os.path.exists(fasta_name):
			with open(fasta_name, 'w') as f:
				for strain in strains:
					f.write('>' + strain + '\n')
					f.write(''.join([r.codons[strain] for r in bin]) + '\n')

	# Save RSA
	if saveRSA:
		save_averages_for_bins(properties,"interfacial_vs_nonInterfacial",RSA_DIR)	

	# # Visualize the binning
	# if plot:
	# 	barPlots(buckets,binEdges,varName,saveDir)

	return df



def bin_interf_vs_nonInterf(modeledPPIs,FASTA_DIR,strains,RSA_DIR,saveRSA = False,plot=False):
	''' Creates directory structure and run the binning for interfacial vs non-interfacial residues in multiple report PPIs
		
	'''
	fasta_path = os.path.join(FASTA_DIR,"interfacial_vs_non_interfacial")
	if not os.path.exists(fasta_path):
		os.makedirs(fasta_path)

	interf= [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True ]
	noninterf =[res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False ]
	
	df = bin_codons_two_bins(interf,"interfacial",noninterf,"non_interfacial",fasta_path,strains,RSA_DIR,saveRSA = saveRSA,plot=plot)
	
	return df

def bin_codons(residueVar, binEdges, modeledPPIs, varName, strains, saveDir,RSA_DIR,residueType,
			   saveRSA=False,readmeTxT=None,plot=False):
	"""Bin codons by a single continuous residue-level variable and save them.

	Saves the bins for each of the protein and residue level categories.

	Args:
		residueVar (function): call to access value of residue variable.
		binEdges (list(float)): edges of the bins.
		modeledPPIs (dict(modeledPPIs))
		varName (str): name of variable that is binned on.
		strains: yeast strains to use in evolutionary rate calculations
		saveDir
		fNames (str): overwrite file name.
		saveRSA (bool): also save the average RSA of each bin
		rsaDir
		residueType (str): type of residues binned (ex: all, interfacial...)

	"""
	print("> Binning by: ",varName)
	if not os.path.exists(saveDir):
		os.makedirs(saveDir)

	if len(binEdges) < 2:
		raise ValueError('Expected at least 2 bin edges')

	residues = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues]
	bins = [(varName+"_"+str(dn) + '_to_' + str(up),[]) for (dn, up) in zip(binEdges[:-1], binEdges[1:])]
	# Put residues into respective bins
	for r in tqdm(residues):
		val = residueVar(r)
		if val == np.nan or type(val)== str:
			continue
		if val > binEdges[-1]:
			raise UserWarning('Overflow from bin')
		for i, (dn, up) in enumerate(zip(binEdges[:-1], binEdges[1:])):
			if (val <= up and val > dn) or (dn == up and val == dn):
				bins[i][1].append(r)

	properties = pd.DataFrame(columns=['bin_name','num_of_residues','RSA_monomer', 'RSA_complex', 'deltaRSA',
				'interRRC', 'dCenter','dEdges','consurfDB_score','consurfDB_std','consurfDB_resNum','consurfRate4Site_score','consurfRate4Site_std','consurfRate4Site_resNum'])

	# Write the summary text file for all bins
	for name,bin in bins:

		row={}
		row['bin_name']			= name
		row['num_of_residues']		= len(bin)
		row['RSA_monomer']			= float((np.nanmean([float(r.rsaMonomer) for r in bin])))
		row['RSA_complex']			= float((np.nanmean([float(r.rsaComplex) for r in bin])))
		row['deltaRSA']				= float((np.nanmean([float(r.dRsa) for r in bin])))
		row['interRRC']				= float((np.nanmean([float(r.interRRC) for r in bin])))
		row['dCenter']				= float((np.nanmean([float(r.dCenter) for r in bin])))
		row['dEdges']				= float((np.nanmean([float(r.dEdges) for r in bin])))
		row['consurfDB_score']		= float((np.nanmean([float(r.consurfDB_score) for r in bin])))
		row['consurfDB_std']		= float((np.nanstd([float(r.consurfDB_score) for r in bin])))
		row['consurfDB_resNum']		= len([float(r.consurfDB_score) for r in bin if not np.isnan(r.consurfDB_score)])
		row['consurfRate4Site_score']	= float((np.nanmean([float(r.consurf_score) for r in bin])))
		row['consurfRate4Site_std']		= float((np.nanstd([float(r.consurf_score) for r in bin])))
		row['consurfRate4Site_resNum']	= len([float(r.consurf_score) for r in bin if not np.isnan(r.consurf_score)])
		properties = properties.append(row, ignore_index=True)

	properties = properties.astype({"num_of_residues": int})
	properties = properties.astype({"consurfDB_resNum": int})
	properties = properties.astype({"consurfRate4Site_resNum": int})

	df= properties[['consurfDB_score','consurfDB_std','consurfDB_resNum','consurfRate4Site_score','consurfRate4Site_std','consurfRate4Site_resNum']]
	
	
	readmeTxT="#"+residueType+" residues in modeled S.cerevisiae PPIs split into "+str(len(bins))+" bins of "+varName+"\n#\n#Summary:\n"+str(properties.to_string())
	
	
	if not readmeTxT==None:
		ofile = os.path.join(saveDir,'README.txt')
		if not os.path.exists(ofile):
			with open(ofile, 'w') as f:
				f.write(readmeTxT)
	
	# Write the fasta file for each bin
	for name,bin in tqdm(bins):
		if len(bin) == 0:
			raise UserWarning('Empty bin: ',varName + '_' + str(dn) + '_to_' + str(up)+".fasta")
			return
		fasta_name = os.path.join(saveDir,str(name)+".fasta")	
		if not os.path.exists(fasta_name):
			with open(fasta_name, 'w') as f:
				for strain in strains:
					f.write('>' + strain + '\n')
					f.write(''.join([r.codons[strain] for r in bin]) + '\n')

	# Save RSA
	if saveRSA:
		save_averages_for_bins(properties,varName+"_bins",RSA_DIR)	
	
	# # Visualize the binning
	# if plot:
	# 	barPlots(buckets,binEdges,varName,saveDir)

	return df


def bin_structure(modeledPPIs,FASTA_DIR,strains,RSA_DIR,saveRSA = False,plot=False):
	'''Creates directory structure and run the binning for structural properties '''
	# 1) Directory set up
	if not os.path.exists(FASTA_DIR):
		os.makedirs(FASTA_DIR)
	# fasta_path = os.path.join(FASTA_DIR,"interfacial_vs_non_interfacial")
	# if not os.path.exists(fasta_path):
	# 	os.makedirs(fasta_path)


	# 2) Define the edges of the bins for all categories 
	# The value val is included in bin ]dn,up] if (val <= up and val > dn) or (dn == up and val == dn)
	binEdgesRSA_interface = np.arange(0.00, 1.01, 0.05)			# [0,0] bin is empty for interfaces
	binEdgesRSA = np.insert(binEdgesRSA_interface, 0, 0.00) 	# adding a [0,0] bin for all residues
	binEdgesDeltaRSA= np.arange(0.00, 1.01, 0.1)
	binEdgesDeltaRSA= np.insert(binEdgesDeltaRSA, 0, 0.00) 		# adding a [0,0] bin
	binEdgesContacts = [0, 0, 2, 4, 6, 8, 10, 12, 14, 30]
	binEdgesDistance_from_center = [0,10,15,20,25,30,40,50,65,93]
	binEdgesDistance_from_edges = [0,10,15,20,25,30,35,45,65,100]
	
	# 3) Check numbers and print disclaimer
	allProteins = [prot for ppi in modeledPPIs.values() for prot in ppi.proteins]
	proteinsWithInterface = [prot for prot in allProteins if len(prot.interface)>0]
	proteinsWithoutInterface = [prot for prot in allProteins if len(prot.interface)==0]
	print(len(proteinsWithInterface),"/" ,len(allProteins),"proteins have a modeled interface.")
	print(len(proteinsWithoutInterface),"/" ,len(allProteins),"proteins do not have a modeled interface. Discarded from the interface involvement analysis")

	# 4) Start binning
# - RSA (monomer) 
	# All residues
	bin_codons(lambda r: r.rsaMonomer, binEdgesRSA, modeledPPIs, 'rsa_monomer_all', strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_all"),RSA_DIR=RSA_DIR,residueType="All",saveRSA = True, 
				readmeTxT = "Monomer RSA, for all residues in all proteins.", plot=False) 
	# Only interfacial residues (drsa != 0)
	bin_codons(lambda r: r.rsaMonomer if (r.dRsa > 0.0) else (np.nan), binEdgesRSA_interface, modeledPPIs, 'rsa_monomer_interfacial', strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_interfacial"),RSA_DIR=RSA_DIR,residueType="Interfacial",saveRSA = True, 
				readmeTxT = "Monomer RSA, for all interfacial residues in all proteins.", plot=False) 
	# Only non-interfacial residues (drsa == 0)
	bin_codons(lambda r: r.rsaMonomer if (r.dRsa == 0.0) else (np.nan), binEdgesRSA, modeledPPIs, 'rsa_monomer_non_interfacial',strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_monomer","rsa_monomer_non_interfacial"),RSA_DIR=RSA_DIR,residueType="Non-interfacial",saveRSA = True, 
				readmeTxT = "Monomer RSA, for all non-interfacial residues in all proteins.", plot=False) 

# - RSA (complex) 
	# All residues
	bin_codons(lambda r: r.rsaComplex, binEdgesRSA, modeledPPIs, 'rsa_complex_all', strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_all"),RSA_DIR=RSA_DIR,residueType="All",saveRSA = True, 
				readmeTxT="Complex RSA, for all residues in all proteins.",plot=False) 
	# Only interfacial residues (drsa != 0)
	bin_codons(lambda r: r.rsaComplex if (r.dRsa > 0.0) else (np.nan), binEdgesRSA, modeledPPIs, 'rsa_complex_interfacial' , strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_interfacial"),RSA_DIR=RSA_DIR,residueType="Interfacial",saveRSA = True, 
				readmeTxT="Complex RSA, for all interfacial residues in all proteins.",plot=False) 
	# # Only non-interfacial residues (drsa == 0)
	bin_codons(lambda r: r.rsaComplex if (r.dRsa == 0.0) else (np.nan), binEdgesRSA, modeledPPIs, 'rsa_complex_non_interfacial', strains,
				saveDir= os.path.join(FASTA_DIR,"rsa_complex","rsa_complex_non_interfacial"),RSA_DIR=RSA_DIR,residueType="Non-interfacial",saveRSA = True, 
				readmeTxT="Complex RSA, for all non-interfacial residues in all proteins.",plot=False) 

# - Measures of interface involvement 	
	# - Delta RSA 
	bin_codons(lambda r: r.dRsa, binEdgesDeltaRSA, modeledPPIs,'drsa', strains,
				saveDir= os.path.join(FASTA_DIR,"drsa"),RSA_DIR=RSA_DIR,residueType="All",saveRSA = True, 
				readmeTxT='deltaRSA (Difference in solvent accessible surface between complex and monomer state) for all residues in all proteins. [0,0] bin represents non-interfacial residues.', 
				plot=False) 
	# - Residue contacts
	bin_codons(lambda r: r.interRRC, binEdgesContacts, modeledPPIs,'residueContacts', strains,
				saveDir= os.path.join(FASTA_DIR,"residueContacts"),RSA_DIR=RSA_DIR,residueType="All",saveRSA = True,
				readmeTxT='Residue contacts (Number of contacts with residues in partner protein) for all residues in all proteins. We only have residue contact values for interfacial residues in proteins with a modeled interface', 
				plot=False) 
	# - Distance from center
	bin_codons(lambda r: r.dCenter, binEdgesDistance_from_center, modeledPPIs,'distance_from_center', strains,
				saveDir= os.path.join(FASTA_DIR,"distance_from_center"), RSA_DIR=RSA_DIR,residueType="All",saveRSA = True,
				readmeTxT='Binning by distance to geometric center of the interface for all residues in all proteins. We only have distance from center values for interfacial residues in proteins with a modeled interface', 
				plot=False) 
	# - distance from edges
	bin_codons(lambda r: r.dEdges, binEdgesDistance_from_edges, modeledPPIs,'distance_from_edges',strains,
				saveDir= os.path.join(FASTA_DIR,"distance_from_edges"), RSA_DIR=RSA_DIR,residueType="All",saveRSA = True,
				readmeTxT='Binning by distance to edges (most exposed surface in complex state) of the interface for all residues in all proteins. We only have distance from edges values for interfacial residues in proteins with a modeled interface', 
				plot=False) 		
	
	return

def write_dir_list(FASTA_DIR,parallelize = False):
	'''
	FASTA_DIR: binned_codons directory
	TO EDIT
	NOT WORKING
	'''

	# TODO: Add and check NOT parallelized
	dirs = []
	subdirs = [os.path.join(FASTA_DIR,dir) for dir in os.listdir(FASTA_DIR) if os.path.isdir(os.path.join(FASTA_DIR,dir))]
	subsubdirs = [os.path.join(subdir,subsubdir) for subdir in subdirs for subsubdir in os.listdir(subdir) if os.path.isdir(os.path.join(subdir,subsubdir))]
	subsubdirs = subdirs + subsubdirs
	
	files= []
	for dir in subsubdirs:
		files = files+[os.path.join(dir,file) for file in os.listdir(dir) if file[0:6] != "README"]
	print(files)
	print(len(files))
	exit()
	if parallelize:
		for i,k,fileNum in zip(files[0::2], files[1::2],range(len(files))):
			fileName = os.path.join(FASTA_DIR,"dir_list"+str(fileNum)+".txt")	
			with open(fileName, 'w') as f:
				f.write(i+"\n")
				f.write(k+"\n")
	else:
		# Check?
		for dir,filenum in zip(subsubdirs,range(len(subsubdirs))):
			fileName = os.path.join(FASTA_DIR,"dir_list"+str(fileNum)+".txt")	
			with open(fileName, 'w') as f:
				f.write(dir+"\n")
				
	return

def main():

	path_closely_related = '../data/processed/modeled_ppis.pkl'
	path_additional_closely_related = '../data/processed/modeled_ppis_additional_closely_related.pkl'
	path_additional_distantly_related = '../data/processed/modeled_ppis_additional_distantly_related.pkl'

	closely_related = ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']
	additional_closely_related = ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Ctro', 'Tpha','Cfab','Ndai','Zrou']
	additional_distantly_related = ['Scer','Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Fgra','Ncra','Anid','Spom','Nirr','Plac','Pjir','Pmur']
	
	names= ["","_AdditionalCloseSpecies","_AdditionalDistantSpecies"]
	# Run binning once per model
	for name,modelPath,strains in zip(names,[path_closely_related,path_additional_closely_related,path_additional_distantly_related],[closely_related,additional_closely_related,additional_distantly_related]):
		# Set up
		OUT_DIR = '../data/evolutionary_rates'+name
		FASTA_DIR = os.path.join(OUT_DIR, 'binned_codons')
		RSA_DIR = os.path.join(OUT_DIR, 'RSA_distributions')
		if not os.path.exists(FASTA_DIR):
			os.makedirs(FASTA_DIR)
		# # Load model
		# print ("loading protein models")
		# with open(modelPath, 'rb') as f:
		# 	modeledPPIs = pickle.load(f)
		# print ("protein models loaded")

		# Ignore "RuntimeWarning: Mean of empty slice" due to non-interfacial residues not having some structural prooperties computed
		with warnings.catch_warnings(): 
			warnings.simplefilter("ignore")

			# # 1) Bin interfacial vs. non interfacial residues
			# df1 = bin_interf_vs_nonInterf(modeledPPIs,FASTA_DIR,strains,RSA_DIR,saveRSA = True,plot=True)
			
			# # 2) Bin by structural properties
			# bin_structure(modeledPPIs,FASTA_DIR,strains,RSA_DIR,saveRSA = True,plot=True)
			# # bin_all(modeledPPIs,FASTA_DIR,RSA_DIR,strains)
			
			# 3) Write the directory list
			write_dir_list(FASTA_DIR,parallelize = False)


	# TODO:
	# Visualize the binning?
	# Summary?
	

	# 	# # 4) Write dir_list files to run dN/dS calculations
	# 	# # 1 file / folder if parallelize = True
	# 	# write_dir_list(parallelize = True)

	# 	# 5) Save summary
	# 	save_summary(scerMult1,scerMult2,scerMult3,scerMult4,scerMult5,spomMult1,spomMult2,spomMult3,spomMult4,spomMult5,scer1,scer2,scer3,scer4,scer5,spom1,spom2,spom3,spom4,spom5)
		
if __name__ == '__main__':
	main()
