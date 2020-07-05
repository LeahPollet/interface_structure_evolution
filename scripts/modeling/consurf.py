"""
Download consurfDB results for a list of PDBs / chain IDs and adds them to PPI models

Input   = 
	curated_scer_pdb_homologs.txt - map of s. cer. proteins to homologous PDB chains
	perl scrip from http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip
Output  = consurfDB_input.txt right format for batch consurf DB query

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os, sys, stat
from modeled_proteins import Protein, Residue
import pickle
from collections import Counter
import re
from tqdm import tqdm


OUT_DIR = '../../data/consurf' 
INFILES_DIR = os.path.join(OUT_DIR, 'input_file')

PERL_SCRIPT_DIR = os.path.join(OUT_DIR, 'perl_script')
CONSURF_DB_DIR = os.path.join(OUT_DIR, 'consurfDB_files')
PROTEIN_MODELS_PATH = '../../data/processed/modeled_proteins.pkl'

url = "http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip"

def prep_input(): 
	'''Create the consurf directory structure and the input file for consurf analysis (consurfDB_input.txt: list of pdbID and pdbChain matching to S. cere ORFs)
	Returns: 
		idsList (list): list of the combined pdbID+chainID in the right format to query the consurf database
		proteins (dict): PPI models from the build_protein_models.py script
	'''
	print("--------Loading protein models-------- ")
	modelName = os.path.basename(PROTEIN_MODELS_PATH)
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	if not os.path.exists(INFILES_DIR):
		os.makedirs(INFILES_DIR)
	idsList = []
	infile = os.path.join(INFILES_DIR, 'consurfDB_input_'+str(modelName)+'.txt')

	if not os.path.exists(PROTEIN_MODELS_PATH):
		raise UserWarning('Protein models not created')
		
	with open(PROTEIN_MODELS_PATH, 'rb') as f:
		proteins = pickle.load(f)
		print ("Loaded protein models")

	fOut = open(infile, 'w')
	for protein in proteins.values():
		pdbID = protein.pdb[0]
		chainID = protein.pdb[1]
		DB_query_ID = pdbID.upper()+chainID
		idsList.append(DB_query_ID)
		fOut.write(' '.join([pdbID, chainID])+'\n')
	fOut.close()
		
	return idsList, proteins

def download_consurf_results(idsList):
	'''
	Download consurf result files from the ConSurf database and prints a log file
	Args:
		idsList (list):list of the combined pdbID+chainID in the right format to query the consurf database

	Returns:
		 hasConsurfResults (list): list of all of the (pdbID,chainID)s with ConSurf scores downloaded
	'''
	print("--------Downloading ConSurf files-------- ")
	if not os.path.exists(CONSURF_DB_DIR):
		os.makedirs(CONSURF_DB_DIR)

	modelName = os.path.basename(PROTEIN_MODELS_PATH)
	already_downloaded = [(filename[:4],filename[5]) for filename in os.listdir(CONSURF_DB_DIR)]
	logfile = os.path.join(OUT_DIR, 'ConSurfDB_download_'+str(modelName)+'.log')
	failed_downloads = []
	if os.path.exists(logfile):
		failed_downloads = [(line.split("_")[0],line.split("_")[1]) for line in open(logfile, 'r')]
	log = open(logfile, 'a')
	hasConsurfResults = []

	for queryId in tqdm(idsList):
		pdbID = queryId[:4]
		chainID = queryId[4]
		
		if (pdbID,chainID) in failed_downloads:
			continue
		elif (pdbID,chainID) not in already_downloaded+hasConsurfResults:
			
			queryURL = "https://consurfdb.tau.ac.il/DB/"+queryId+"/consurf_summary.txt"
			downloaded_file_Name = pdbID+"_"+chainID+"_consurf_summary.txt"
			downloaded_file_path = os.path.join(CONSURF_DB_DIR, downloaded_file_Name)
			if not os.path.exists(downloaded_file_path):
				os.system('wget -q -O'+downloaded_file_path+' '+queryURL)
				if os.stat(downloaded_file_path).st_size == 0: # If download failed
					log.write(downloaded_file_Name+" -> Unable to find the grades page on ConSurfDB site\n")
					os.remove(downloaded_file_path)
				else:
					hasConsurfResults.append((pdbID,chainID))
		else:
			hasConsurfResults.append((pdbID,chainID))

	print("Consurf files downloaded.")
	
	print(len(hasConsurfResults),"/",len(idsList)," proteins in the model have consurf results ")
	return hasConsurfResults

def add_to_protein_models(has_consurf_results,proteins):
	'''
	Add the consurf results to proteins in our protein models and prints a log file
	Args:
		has_consurf_results (list): list of all of the (pdbID,chainID)s with ConSurf scores downloaded
		proteins (dict):  PPI models from the build_protein_models.py script
	'''

	print("--------Adding to PPI models-------- ")
	modelName = os.path.basename(PROTEIN_MODELS_PATH)
	logfile = os.path.join(OUT_DIR, 'match_with_models_'+str(modelName)+'.log')
	if not os.path.exists(logfile):
		
		proteins_with_errors = {}
		total_residue_numbers = {}
		double_counted = {}      # number of time a PDB is used, to avoid double counting errors in the log file

		# --- For all proteins modeled
		for protein in tqdm(proteins.values()):
			if (protein.pdb[0].upper(),protein.pdb[1]) in has_consurf_results:

				ID = protein.pdb[0].upper()+"_"+protein.pdb[1] # Format: PDBID_chainID
				# print ("> Curent protein for consurf update:", protein.pdb[0].upper(),"Chain:",protein.pdb[1])

				# 1) Updates the stats dictionaries
				consurf_file = os.path.join(CONSURF_DB_DIR,ID +'_consurf_summary.txt')
				Prot_and_PDB = str(protein.orf+", "+ID) # ORFNAME,PDBID_chainID (string)
				total_residue_numbers.update({Prot_and_PDB:len(protein.residues)})
				double_counted[Prot_and_PDB] = double_counted.get(Prot_and_PDB, 0) + 1 # Add the PDBid to the double_counted dict if first time, add 1 to its count if already encountred the id
				
				# 2) Load consurf results
				consurf_summary = get_results_from_consurf_file(consurf_file)
				
				# 3) Update protein models and error log dictionary
				# -> Check if all the chains in the Consurf result file are the same (i.e. Consurf results from only one PDB file chain)
				chains = set([val[0] for val in consurf_summary.values()])
				if len(chains)!=1:
					proteins_with_errors.setdefault(Prot_and_PDB,[]).append("Different chains used in PDB")
					continue
				chain= list(chains)[0] # Get the chain used

				# -> For each residue
				for residue in protein.residues:
					id = residue.id

					# - If another chain was used in the consurf calculations
					# Remove the current chain used from the id and replace it with the one used by consurf
					if chain != id.split(':', 1)[-1]:
						id = id.split(':', 1)[0]+':'+ chain

					# - If we have consurf results for this residue in the dictionary
					if id in consurf_summary:

						# Below confidence cutoff error
						if consurf_summary[id][5] == "yes": # Below the confidence cutoff
							proteins_with_errors.setdefault(Prot_and_PDB,[]).append("Score below ConSurf confidence cutoff")
							continue
						# Add to models
						residue.consurf_score = consurf_summary[id][3]	
						residue.conf_int_consurf_score = consurf_summary[id][4]	
						

					# - No results
					else:
						proteins_with_errors.setdefault(Prot_and_PDB,[]).append("Consurf file format error")

		# --- Write the log file
		write_log_file(proteins_with_errors,double_counted,total_residue_numbers,logfile)

		# --- Check counts in the models:
		total_prot 		= len(proteins.values())
		score_for_prot 	= len([p for p in proteins.values() if len(p.residues) != len([r for r in p.residues if r.consurf_score != r.consurf_score])])
		total_res 		= len([r for p in proteins.values() for r in p.residues])
		score_for_res 	= len([r for p in proteins.values() for r in p.residues if r.consurf_score == r.consurf_score]) # Only nan is not equal to itself
		
		print(score_for_prot,"/",total_prot," proteins in the PPI models have no consurf score at all.")
		print(score_for_res,"/",total_res," residues in the PPI models have no consurf score.")

		# --- Save the new protein models updated with consurf score
		pickle.dump(proteins, open(PROTEIN_MODELS_PATH, 'wb'))


	else:
		print("Consurf score already added to models. Skipping")

def get_results_from_consurf_file(consurf_file):
	# Read out the results from a consurf file to a dictionary

	consurf_summary = {} # Dictionary to save the consurf results

	with open(consurf_file, 'r') as h:
		rows = list((line.split() for line in h))
		rows = rows[15:-4] # skip the header and footer

		for row in rows:
			atom = row[2] # Format: MET2:B, need to split into residue, number and chain
			# Amino acids that are not in the pdb file but somehow still have a score. Skipping them for now.
			if atom == "-":
				continue
			chain = atom [-1:]
			number = float(re.findall(r'\d+',atom[:-1])[0])
			residue = atom[:3]
			score = float(row[3])
			# File format sometimes splits the confidence interval 
			if len(row)== 10: 
				confidence_interval = eval(row[5]+row[6])
			else: # len(row ==9)
				confidence_interval = eval(row[5])

			below_cutoff = "no"
			if '*' in row[4]:
				below_cutoff ="yes"

			consurf_summary.update({atom: (chain,number,residue,score,confidence_interval,below_cutoff)})
					
	return consurf_summary
	
def write_log_file(proteins_with_errors,double_counted,total_residue_numbers,logfile):
	# Write the log file for adding to models

	with open(logfile, 'a') as f:
		f.write("ORF_name, PDB_id_used ---> error_message_1: number_of_error_messages_1/total_num_res_in_prot ---> error_message_2: number_of_error_messages_2/total_num_res_in_prot ...\n" )
		for Prot_and_PDB in proteins_with_errors:
			error_counts = list(Counter(proteins_with_errors[Prot_and_PDB]).items())
			
			f.write(Prot_and_PDB)
			for i in range(len(error_counts)):
				number_or_times_used = double_counted[Prot_and_PDB]
				error_message = error_counts[i][0]
				number_of_error_messages = error_counts[i][1]
				corrected_number_of_error_messages =  number_of_error_messages/number_or_times_used
				total_num_res_in_prot = total_residue_numbers[Prot_and_PDB]
				f.write(" ---> " + error_message+": "+str(corrected_number_of_error_messages)+"/"+str(total_num_res_in_prot))
			f.write("\n")
		
	return


def main():
	idsList,proteins = prep_input()
	hasConsurfResults = download_consurf_results(idsList)
	add_to_protein_models(hasConsurfResults,proteins)

if __name__ == '__main__':
    main()

