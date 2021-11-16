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
from modeled_PPIs import PPI,Protein,Interface,Residue
import pickle
from collections import Counter
import re
from tqdm import tqdm
import numpy as np


OUT_DIR = '../data/consurfDB' 
INFILES_DIR = os.path.join(OUT_DIR, 'input_file')

# PERL_SCRIPT_DIR = os.path.join(OUT_DIR, 'perl_script')
CONSURF_DB_DIR = os.path.join(OUT_DIR, 'consurfDB_files')

url = "http://bental.tau.ac.il/new_ConSurfDB/ConSurfDB_list_feature.zip"

def prep_input(models): 
	'''Create the consurf directory structure and the input file for consurf analysis (consurfDB_input.txt: list of pdbID and pdbChain matching to S. cer ORFs)
	Args:
		models(str): 	list of paths to pickled PPI model to download ConSurf results for
	Returns: 
		idsList (list): list of the combined pdbID+chainID in the right format to query the consurf database
		proteins (dict): PPI models from the build_protein_models.py script
	'''
	print("--------Preping input-------- ")
	if not os.path.exists(OUT_DIR):
		os.makedirs(OUT_DIR)
	if not os.path.exists(INFILES_DIR):
		os.makedirs(INFILES_DIR)

	loaded_models = []
	pdbID_chainID = []
	idsList = []

	for model_path in tqdm(models):
		model = pickle.load(open(model_path,"rb"))
		loaded_models.append(model)
		pdbID_chainID = pdbID_chainID+[protein.pdb for ppi in model.values() for protein in ppi.proteins]
		idsList = idsList + [key[0].upper()+key[1] for key in pdbID_chainID]

	# Remove duplicates
	pdbID_chainID = list(set(pdbID_chainID))
	idsList= list(set(idsList))

	infile = os.path.join(INFILES_DIR, 'consurfDB_input.txt')
	fOut = open(infile, 'w')
	for pdbID,chainID in pdbID_chainID:
		fOut.write(' '.join([pdbID, chainID])+'\n')
	fOut.close()
	
	return idsList, loaded_models


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


	already_downloaded = [(filename[:4],filename[5]) for filename in os.listdir(CONSURF_DB_DIR)]
	print("Already downloaded:",len(already_downloaded))
	# print(already_downloaded[0])
	test_idsList = [(queryId[:4],queryId[4]) for queryId in idsList]
	print("To download:",len(test_idsList))
	# print(test_idsList[0])
	logfile = os.path.join(OUT_DIR, 'ConSurfDB_download.log')
	failed_downloads = [(line.split("_")[0],line.split("_")[1]) for line in open(logfile, 'r')]
	print("Failed downloads:",len(failed_downloads))
	# print(failed_downloads[0])
	#
	print("No data:",len([x for x in test_idsList if x not in already_downloaded+failed_downloads]))
	#

	# Doing with what i have already downloaded for now TO edit...
	hasConsurfResults = [x for x in test_idsList if x in already_downloaded]
	return hasConsurfResults
	exit()



	logfile = os.path.join(OUT_DIR, 'ConSurfDB_download.log')
	failed_downloads = []
	if os.path.exists(logfile):
		failed_downloads = [(line.split("_")[0],line.split("_")[1]) for line in open(logfile, 'r')]
	log = open(logfile, 'a')
	print(failed_downloads)
	hasConsurfResults = []

	for queryId in tqdm(idsList):
		pdbID = queryId[:4]
		chainID = queryId[4]
		
		if (pdbID,chainID) in failed_downloads:
			continue
		if (pdbID,chainID) in hasConsurfResults:
			continue
		if (pdbID,chainID) in already_downloaded:
			hasConsurfResults.append((pdbID,chainID))
			continue
		
		# Need to download:
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
		

	print("Consurf files downloaded.")
	
	print(len(hasConsurfResults),"/",len(idsList)," proteins in checked have consurf results ")
	return hasConsurfResults

def add_to_protein_models(has_consurf_results,models,paths):
	'''
	Add the consurf results to proteins in our PPI models and prints a log file with errors
	Args:
		has_consurf_results (list): list of all of the (pdbID,chainID)s with ConSurf scores downloaded
		models (list):  		list of loaded PPI models 
		paths(list):			lits of paths to the models
	'''

	for i,model in zip([0,1,2], models):
		name = os.path.basename(paths[i])
		logfile = os.path.join(OUT_DIR, 'errorsAddingTo_'+name+'.log')

		if not os.path.exists(logfile):
			proteins_with_errors = {} 	# List of error messages
			double_counted = {}      	# number of time a PDB is used, to avoid double counting errors in the log file

			modeledProteins = [protein for ppi in model.values() for protein in ppi.proteins]
			
			# --- For all proteins modeled
			for protein in tqdm(modeledProteins):
				ID = protein.pdb[0].upper()+"_"+protein.pdb[1] # Format: PDBID_chainID
				Prot_and_PDB = str(protein.orf+", "+ID) # ORFNAME,PDBID_chainID (string)
				# If we have a ConsurfDB file downloaded
				if (protein.pdb[0].upper(),protein.pdb[1]) in has_consurf_results:
					# Get Consurf results
					consurf_file = os.path.join(CONSURF_DB_DIR,ID +'_consurf_summary.txt')
					double_counted[Prot_and_PDB] = double_counted.get(Prot_and_PDB, 0) + 1 # Add the PDBid to the double_counted dict if first time, add 1 to its count if already encountred the id
					consurf_summary = get_results_from_consurf_file(consurf_file)

					# Update PPI models and error log dictionary
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
							residue.consurfDB_score = consurf_summary[id][3]	
							residue.consurfDB_conf_int = consurf_summary[id][4]	
							
						# - No results
						else:
							proteins_with_errors.setdefault(Prot_and_PDB,[]).append("Consurf file format error")

				# No file on the Consurf DB
				else:
					proteins_with_errors.setdefault(Prot_and_PDB,[]).append("No consurfDB file")
					double_counted[Prot_and_PDB] = double_counted.get(Prot_and_PDB, 0) + 1 # Add the PDBid to the double_counted dict if first time, add 1 to its count if already encountred the id
			 
			 # --- Write the log file
			write_log_file(proteins_with_errors,double_counted,logfile,modeledProteins)

			# --- Check counts in the models:
			total_res 		= len([res for prot in modeledProteins for res in prot.residues])
			score_for_res 	= len([res for prot in modeledProteins for res in prot.residues if res.consurfDB_score != res.consurfDB_score]) 
		
			print(score_for_res,"/",total_res," residues in the PPI models have a consurf DB score.")
			
			# --- Save the new protein models updated with consurf score
			pickle.dump(model, open(paths[i], 'wb'))	
		

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
	
def write_log_file(proteins_with_errors,double_counted,logfile,modeledProteins):
	
	# Get the number of residues for every modeled protein
	resNum = {str(protein.orf+", "+protein.pdb[0].upper()+"_"+protein.pdb[1]) : len(protein.residues) for protein in modeledProteins}
	
	# Write the log file for adding to models
	with open(logfile, 'a') as f:
		f.write("ORF_name, PDB_id_used: error_message (number_of_error_messages_1/total_num_res_in_prot) \n" )
		
		for Prot_and_PDB in proteins_with_errors:
			error_counts = list(Counter(proteins_with_errors[Prot_and_PDB]).items())
			
			f.write(Prot_and_PDB+": ")
			for i in range(len(error_counts)):
				error_message = error_counts[i][0]
				number_or_times_used = double_counted[Prot_and_PDB]
				number_of_error_messages = error_counts[i][1]
				total_num_res_in_prot = resNum[Prot_and_PDB]
				if error_message == "No consurfDB file":
					corrected_number_of_error_messages = resNum[Prot_and_PDB]
				else:
					corrected_number_of_error_messages =  int(number_of_error_messages/number_or_times_used)

				f.write(error_message+" ("+str(int(corrected_number_of_error_messages))+"/"+str(total_num_res_in_prot)+" residues). ")
			f.write("\n")
		
	return

def print_summary(modeledPPIs,txt):
	print("\n>",txt,":")

	totalRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues]
	interfacialRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == True ]
	nonInterfacialRes = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues if res.interfacial == False ]

	print("Total number of modeled residues:", len(totalRes))
	print("(",len(nonInterfacialRes),"non-interfacial +",len(interfacialRes),"interfacial )")
	total_resWithConSurfDB = [res for res in totalRes if not np.isnan(res.consurfDB_score)]
	intf_resWithConSurfDB = [res for res in interfacialRes if not np.isnan(res.consurfDB_score)]
	nonInt_resWithConSurfDB = [res for res in nonInterfacialRes if not np.isnan(res.consurfDB_score)]

	print(" >Model has ConSurf DB info for: ")
	print(len(total_resWithConSurfDB),"/",len(totalRes),"residues overall")
	print(len(nonInt_resWithConSurfDB),"/",len(nonInterfacialRes),"non-interfacial residues")
	print(len(intf_resWithConSurfDB),"/",len(interfacialRes),"interfacial residues")

	print(" >Average ConSurf DB score:")
	print("Overall average ConSurf score:",np.mean([res.consurfDB_score for res in total_resWithConSurfDB]))
	print("Non-interfacial residues average ConSurf score:",np.mean([res.consurfDB_score for res in nonInt_resWithConSurfDB]))
	print("Interfacial residues average ConSurf score:",np.mean([res.consurfDB_score for res in intf_resWithConSurfDB]))
	print("\n")
	
	return


def main():
	path_closely_related = '../data/processed/modeled_ppis.pkl'
	path_additional_closely_related = '../data/processed/modeled_ppis_additional_closely_related.pkl'
	path_additional_distantly_related = '../data/processed/modeled_ppis_additional_distantly_related.pkl'

	

	idsList,[model_closely_related,model_additional_closely_related,model_additional_distantly_related] = prep_input([path_closely_related,path_additional_closely_related,path_additional_distantly_related])
	
	hasConsurfResults = download_consurf_results(idsList)
	add_to_protein_models(hasConsurfResults,[model_closely_related,model_additional_closely_related,model_additional_distantly_related],[path_closely_related,path_additional_closely_related,path_additional_distantly_related])
	
	print_summary(model_closely_related,"PPI models for closely related species")
	print_summary(model_additional_closely_related,"PPI models for additional closely related species")
	print_summary(model_additional_distantly_related,"PPI models for additional distantly related species")

if __name__ == '__main__':
    main()


