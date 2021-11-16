"""
Manually compute consurf scores using the rate4site program

"""
import os, sys, stat
from modeled_PPIs import PPI,Protein,Interface,Residue
import pickle
from collections import Counter
import re
from tqdm import tqdm
import zipfile
import subprocess
from shutil import copyfile
from Bio import SeqIO
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
sns.set(style="whitegrid")


def download_rate4sites(): 
	# useful if rate4site is not compling properly:
	# https://www.tau.ac.il/~itaymay/cp/rate4site.html
	oPath = "../data"
	url = "https://www.tau.ac.il/~itaymay/cp/rate4site.3.2.source.zip"

	# Download the source code
	zipped = oPath+"/"+os.path.basename(url)
	path_to_rate4site = oPath+"/"+os.path.basename(url)[:-4]
	if not os.path.exists(path_to_rate4site):
		print("Downloading rate4site")
		os.system('wget -P '+ oPath +' '+url)
		with zipfile.ZipFile(zipped, 'r') as zip_ref:
			zip_ref.extractall(oPath)
		os.remove(zipped)

	# Compiling the source code
	path_to_rate4site_exec = oPath+"/"+os.path.basename(url)[:-4]+"/sourceMar09/rate4site"
	if not os.path.exists(path_to_rate4site_exec):	
		path = path_to_rate4site_exec[:-10]
		subprocess.call('make', shell=True,cwd=path)
	
	# print(path_to_rate4site_exec)
	return path_to_rate4site+"/sourceMar09"

def set_up_for_rate4sites(name,modelPath,rate4site_dir,inDir,alignment_folder,tree_file): 
	'''Create the consurf directory structure and the input files for consurf analysis (MSA and phylo tree for each protein)
	Returns: 
		inputs (list): list of ORF names with input available for running rate4site
		proteins (dict): PPI models from the build_protein_models.py script
	'''
	print("------ Set up for rate4Site ------ ")
	# Set up the consurf directory
	if not os.path.exists(rate4site_dir):
		os.makedirs(rate4site_dir)
	if not os.path.exists(inDir):
		os.makedirs(inDir)
	
	# Load protein models
	if not os.path.exists(modelPath):
		raise UserWarning('Protein models not created')
	with open(modelPath, 'rb') as f:
		modeledPPIs = pickle.load(f)
		print ("Loaded protein models")
	
	# Create 1 input dir per ORF in the model
	# Uses same tree and alignment file  as in dnds calculations
	# Alignment file seq.id format modifed to match the tree format
	inputs = list(set([protein.orf for ppi in modeledPPIs.values() for protein in ppi.proteins]))
	proteins = [protein for ppi in modeledPPIs.values() for protein in ppi.proteins]
	
	for protein in tqdm(proteins):
		ORF_folder = inDir+"/"+protein.orf
		if not os.path.exists(ORF_folder):
			os.makedirs(ORF_folder)
			# Copy the tree file  and alignment file into the new directory
			copyfile(tree_file, ORF_folder+"/yeasts.tree")
			alignment_file_original = alignment_folder+"/"+protein.orf+"/"+protein.orf+".aa.aln"

			# Copy the alignment file into the new directory, modifying its sequence names to match the tree file
			alignment_file_corrected = ORF_folder+"/"+protein.orf+".aln"
			with open(alignment_file_original) as original, open(alignment_file_corrected, 'w') as corrected:
					records = SeqIO.parse(original, 'fasta')
					for record in records:
						record.id = record.description.split()[1][:4]
						record.description = record.id 
						SeqIO.write(record, corrected, 'fasta')
		
	return inputs, modeledPPIs

def run_rate4sites(name,path_to_rate4site,input_list,inDir,outDir):
	print("-------- Run rate4Site -------- ")
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	for orf in tqdm(input_list):
		
		tree_file = "../../"+inDir+"/"+orf+"/yeasts.tree"
		aln_file = "../../"+inDir+"/"+orf+"/"+orf+".aln"
		outfile = "../../"+outDir+"/"+orf+".consurf.txt"

		if not os.path.exists(os.path.join(path_to_rate4site,outfile)):
			cmd = "./rate4site -s " + aln_file+" -t "+tree_file+" -o "+outfile
			# subprocess.Popen(["./autogen.sh"], stdout=subprocess.PIPE, cwd=path_to_dssp)
			# stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,cwd=path_to_rate4site,
			# subprocess.call("./data/rate4site.3.2.source/sourceMar09/rate4site -s "+aln_file+" -t "+tree_file+" -o "+outfile, shell=True)
			subprocess.call(cmd, shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,cwd=path_to_rate4site)
			
	return 

def add_to_protein_models(inputs,modeledPPIs,modelPath,resultsDir):
	
	'''
	Add the consurf results to proteins in our protein models and prints a log file
	Args:
		inputs (list): list of yeast ORFs with ConSurf scores computed
		modeledPPIs (dict):  PPI models from the build_protein_models.py script
		TODO, add logfile?
	'''
	print("-------- Add to models -------- ")
	# --- For all proteins modeled
	noConsurfScore = 0
	proteins = [protein for ppi in modeledPPIs.values() for protein in ppi.proteins]
	
	for protein in tqdm(proteins):
		if protein.orf in inputs:
			# print(">" ,protein.orf )
			consurf_file = resultsDir+"/"+protein.orf+".consurf.txt"
			consurf_summary = get_results_from_consurf_file(consurf_file)
			for res in protein.residues:
				if (res.indexInSeq in consurf_summary):
					if (res.aa.letter != consurf_summary[res.indexInSeq][0]):
						print("Problem, Sequence mismatch in the mapping of Consurf scores")

					res.consurf_score = consurf_summary[res.indexInSeq][1]
					res.conf_int_consurf_score=consurf_summary[res.indexInSeq][2]
					res.consurf_std =consurf_summary[res.indexInSeq][3]
					
				else:
					noConsurfScore += 1

	# --- Check counts in the models:
	total_res 		= len([r for ppi in modeledPPIs.values() for p in ppi.proteins for r in p.residues])
	score_for_res 	= len([r for ppi in modeledPPIs.values() for p in ppi.proteins for r in p.residues if not np.isnan(r.consurf_score)]) 
	print(score_for_res,"/",total_res," residues in the PPI models have a rate4site consurf score.")
	
	
	# --- Save the new protein models updated with consurf (rate4site) scores
	pickle.dump(modeledPPIs, open(modelPath, 'wb'))

	return

def get_results_from_consurf_file(consurf_file):
	# Read out the results from a consurf file to a dictionary

	consurf_summary = {} # Dictionary to save the consurf results

	with open(consurf_file, 'r') as h:
		rows = list((line.split() for line in h))
		rows = rows[11:-2] # skip the header and footer

		for row in rows:
			if len(row )== 5:
				row=[row[0],row[1],row[2],row[3].split("]")[0]+"]",row[3].split("]")[1],row[4]]
				
			index = float(row[0])-1
			residue = row[1]
			score = float(row[2])
			confidence_interval= eval(row[3])
			std = float(row[4])
			data = eval(row[5])


			consurf_summary.update({index: (residue,score,confidence_interval,std,data)})
					
	return consurf_summary

def print_summary(modeledPPIs,txt):
	print("\n>",txt,":")

	totalRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues]
	interfacialRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True ]
	nonInterfacialRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False ]

	print("Total number of modeled residues:", len(totalRes))
	print("(",len(nonInterfacialRes),"non-interfacial +",len(interfacialRes),"interfacial )")
	total_resWithConSurf = [res for res in totalRes if not np.isnan(res.consurf_score)]
	intf_resWithConSurf = [res for res in interfacialRes if not np.isnan(res.consurf_score)]
	nonInt_resWithConSurf = [res for res in nonInterfacialRes if not np.isnan(res.consurf_score)]

	print(" >Model has ConSurf (rate4site) info for: ")
	print(len(total_resWithConSurf),"/",len(totalRes),"residues overall")
	print(len(nonInt_resWithConSurf),"/",len(nonInterfacialRes),"non-interfacial residues")
	print(len(intf_resWithConSurf),"/",len(interfacialRes),"interfacial residues")

	print(" >Average ConSurf (rate4site) score:")
	print("Overall average ConSurf score:",np.mean([res.consurf_score for res in total_resWithConSurf]))
	print("Non-interfacial residues average ConSurf score:",np.mean([res.consurf_score for res in nonInt_resWithConSurf]))
	print("Interfacial residues average ConSurf score:",np.mean([res.consurf_score for res in intf_resWithConSurf]))
	print("\n")
	
	return

def main():
	path_to_rate4site_exec= download_rate4sites()

	path_closely_related = '../data/processed/modeled_ppis.pkl'
	path_additional_closely_related = '../data/processed/modeled_ppis_additional_closely_related.pkl'
	path_additional_distantly_related = '../data/processed/modeled_ppis_additional_distantly_related.pkl'

	names= ["","_AdditionalCloseSpecies","_AdditionalDistantSpecies"]
	for name,modelPath in zip(names,[path_closely_related,path_additional_closely_related,path_additional_distantly_related]):
		rate4site_dir = '../data/consurf'+name
		inDir = os.path.join(rate4site_dir, 'input_files')
		outDir = os.path.join(rate4site_dir, 'output_files')
		alignment_folder = "../data/external/MSA"+name
		tree_file = "../data/external/yeasts"+name+".tree"

		
		inputs,modeledPPIs = set_up_for_rate4sites(name,modelPath,rate4site_dir,inDir,alignment_folder,tree_file)
		run_rate4sites(name,path_to_rate4site_exec,inputs,inDir,outDir)
		add_to_protein_models(inputs,modeledPPIs,modelPath,outDir)
		print_summary(modeledPPIs,os.path.basename(modelPath))

		# TODO: add a skip if already done for adding to add_to_protein_models (log file check)

if __name__ == '__main__':
    main()




