"""
Script to read and format results from the analysis for plotting

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""
import os,re, sys
import feather
import pickle
import pandas as pd
from modeled_PPIs import PPI,Protein,Interface,Residue
from tqdm import tqdm
import csv
import numpy as np
from scipy import optimize as opt
# from collections import Counter

def linear_fit(x, y, y_err, guess=[1, 1]):
	"""Linear regression for points with uncertainties on y measurement.

	Args:
		x (list(float)): x values.
		y (list(float)): y values.
		y_err (list(float): unceratinty on each y value.
		guess (list(float)): initial value of parameters for optimization.

	Returns:
		dict(float): regression results.

	"""
	if 0. in y_err:
		raise ValueError('can\'t have 0 values in the uncertainties')
	if len(x) != len(y) or len(y) != len(y_err):
		raise ValueError('inconsistent length of lists')
	if len(x) < 3:
		raise ValueError('need more than 2 points')

	def f(x, m, c):
		return (x * m) + c

	coeffs, covMatrix = opt.curve_fit(f, np.array(x), np.array(y),
									  p0=guess,
									  sigma=np.array(y_err),
									  absolute_sigma=True)
	slope, yInt = coeffs
	slopeErr, yIntErr = np.sqrt(np.diag(covMatrix))
	
	return {'slope': slope,
			'intercept': yInt,
			'slope_error': slopeErr,
			'intercept_error': yIntErr
			}

def get_expected_dnds(RSA_DIR,rsaPath,varName,dNdS_DIR):
	'''
	Get the expected fit if constrained only by non-interfacial burrial
	Args:
		RSA_DIR(str): path to folder with average for bins of all structural properties
		rsaPath(str): path to averages for bins of the structural property varName
		varName(str): name of the structural property
		dNdS-DIR(str): path to the folder with dNdS results
	'''
	dNdS_DIR = "/".join(dNdS_DIR.split("/")[0:4])

	# 1) Get the average RSA value of each bin
	RSAs=  pd.read_csv(rsaPath, header=0,usecols=["RSA_complex"])

	# 2) Get expected fit if constrainted only by non-interfacial burrial (non-interfacial complexRSA-dNdS relationship)
	expected_dNdS = read_dnds_values(dNdS_DIR+"/rsa_complex/rsa_complex_non_interfacial", 'rsa_complex_non_interfacial')
	expected_binCenters, expected_yVals, expected_yErr = expected_dNdS['bin_center'].tolist(), expected_dNdS['dnds'].tolist(), expected_dNdS['dnds_err'].tolist()
	expected_fit = linear_fit(expected_binCenters, expected_yVals, expected_yErr)
	
	# 3) Compute the expected values from the fits
	expected_vals = [expected_fit['slope'] * x + expected_fit['intercept'] for x in RSAs["RSA_complex"]]
	
	return expected_vals


def read_dnds_values(inDir, varName):
	"""
	Read dnds values from 1D binning of a variable.

	Args:
		inDir (str): directory with the files containing the dN/dS values
		varName (str): the name of the variable that is used to bin,
					   exactly as appears in the file names

	Returns:
		pandas DataFrame: bin values and dN/dS values

	"""
	df = pd.DataFrame(columns=['bin_name','bin_center','bin_low', 'bin_high','dnds', 'dnds_err','dnds_resNum'])
	for resultsName in os.listdir(inDir):

		resultsPath = os.path.join(inDir, resultsName)
		if os.path.isdir(resultsPath): 			# skip the directories
			continue

		if resultsName.startswith(varName):
			if resultsName.replace('.dnds', '').split('_')[-2] != 'to':
				raise UserWarning(resultsName +' does not follow expected file name spec')
			r = {}
			r['bin_name'] = resultsName.replace('.dnds', '')

			r['bin_high'] = float(resultsName.replace('.dnds', '').split('_')[-1])
			r['bin_low'] = float(resultsName.replace('.dnds', '').split('_')[-3])

			with open(resultsPath, 'r') as f:
				lines = f.readlines()
				r['dnds'] = float(lines[4].split()[1])
				r['dnds_err'] = float(lines[4].split()[4])
				r['dnds_resNum'] = int(lines[0][10:])
				
			df = df.append(r, ignore_index=True)

		df = df.assign(bin_center=(df['bin_high'] - df['bin_low']) / 2. + df['bin_low'])
		df = df.sort_values('bin_center')
		df.reset_index(drop=True, inplace=True) 

	
	return df

def read_dnds_values_2bins(inDir, varName):
	"""
	Read dnds values from 2 bins (interfacial vs non-interfacial)

	Args:
		inDir (str): directory with the files containing the dN/dS values
		varName (str): the name of the variable that is used to bin,
					   exactly as appears in the file names

	Returns:
		pandas DataFrame: bin values and dN/dS values

	"""
	df = pd.DataFrame(columns=['bin_name','dnds', 'dnds_err','dnds_resNum'])

	for resultsName in os.listdir(inDir):

		resultsPath = os.path.join(inDir, resultsName)
		if os.path.isdir(resultsPath): 			# skip the directories
			continue

		if resultsName.startswith('interfacial') or resultsName.startswith('non_interfacial') :
			
			r = {}
			r['bin_name'] = resultsName.replace('.dnds', '')
			
			with open(resultsPath, 'r') as f:
				lines = f.readlines()
				r['dnds'] = float(lines[4].split()[1])
				r['dnds_err'] = float(lines[4].split()[4])
				r['dnds_resNum'] = int(lines[0][10:])

			df = df.append(r, ignore_index=True)
		
	
	return df

def load_consurf(inPath,name,plotDir):
	'''
	Load data into a df and write a summary feather file for R containing
	(1) the value of each structural property (monomer RSA, complex RSA, dRSA, interRRC, dCenter,dEdges)
	(2) the ConSurf (DB) score
	(3) the ConSurf (rate4site) score
	for each residue in the models
	Args:
		inPath (str): 	path to the pickled model
		name (str):		name of the extension to add to the feather file ('' if standard analysis, '_AdditionalDistantSpecies' if adding species for supplementary)
		plotDir (str):	path to the plotting directory
	'''

	path = os.path.join(plotDir,'data_from_models'+name+'.feather')
	if not os.path.exists(path):
		# Load the data from protein models
		print ("loading protein models")
		with open(inPath, 'rb') as f:
			modeledPPIs = pickle.load(f)
		print ("protein models loaded")

		df = pd.DataFrame(columns=['RSA_monomer', 'RSA_complex', 'deltaRSA','interRRC', 
								   'dCenter','dEdges', 'consurfDB_score',
								   'consurfRate4Site_score','consurfRate4Site_std'])
		
		# Create the data frame
		residues = [res for ppi in modeledPPIs.values() for prot in ppi.proteins for res in prot.residues]
		print(len(residues), "total residues in model")

		for r in tqdm(residues):
			res = {}
			res['RSA_monomer']			= float(r.rsaMonomer)
			res['RSA_complex']			= float(r.rsaComplex)
			res['deltaRSA']				= float(r.dRsa)
			res['interRRC']		= float(r.interRRC)
			res['dCenter']	= float(r.dCenter)
			res['dEdges']	= float(r.dEdges)
			res['consurfDB_score']		= float(r.consurfDB_score)
			res['consurfRate4Site_score'] = float(r.consurf_score)
			res['consurfRate4Site_std'] = float(r.consurf_std)
			
			
			df = df.append(res, ignore_index=True)
		
		# Save the dataframe
		feather.write_dataframe(df, path)
		
	# Feather file already exists
	else:
		
		df = pd.read_feather(path)
		
	return df

def load_interfacial_vs_nonInterfacial(inDir,name,plotDir,RSA_DIR,addExpected=False):
	'''
	Load data into a df and write a summary feather file for R containing
	(1) dNds score, error, number of residues, expected dNdS if constrained only by burrial
	(2) ConSurf rate4site score, error, number of residues
	(3) ConSurf rate4site score, error, number of residues
	for two bins: interfacial vs non-interfacial residues
	Args:
		inDir (str): 	path to the folder with dNdS results
		name (str):		name of the extension to add to the feather file ('' if standard analysis, '_AdditionalDistantSpecies' if adding species for supplementary)
		plotDir (str):	path to the plotting directory
	'''
	path = os.path.join(plotDir,'interfacial_vs_non_interfacial.feather')
	if not os.path.exists(path):
		# 1) dN/dS
		dnds_dir = os.path.join(inDir, 'interfacial_vs_non_interfacial')	
		dnds_df = read_dnds_values_2bins(dnds_dir, 'interfacial_vs_non_interfacial') 
		
		# 2) ConSurf
		rsaPath = os.path.join(RSA_DIR, 'averages_for_interfacial_vs_nonInterfacial.csv')
		columns = ["bin_name","consurfDB_score","consurfDB_std","consurfDB_resNum","consurfRate4Site_score","consurfRate4Site_std","consurfRate4Site_resNum"]
		consurf_df = pd.read_csv(rsaPath, header=0,usecols=columns)
		
		# 3) Combined
		intf_vs_nnintf_df = pd.merge(dnds_df,consurf_df, how='left',on=["bin_name"])

		# 4) Add expected fit of constrained only by non-interfacial burrial 
		if addExpected:
			expected = get_expected_dnds(RSA_DIR,rsaPath,'interfacial_vs_non_interfacial',inDir)
			intf_vs_nnintf_df['dnds_expected'] = expected 

		# 5) Write the file
		feather.write_dataframe(intf_vs_nnintf_df, path)

	else:
		intf_vs_nnintf_df = pd.read_feather(path)
		
	return intf_vs_nnintf_df

def load_structural_prop_bins (inDir,name,plotDir,varName,RSA_DIR,addExpected=False):
	'''
	Load data into a df and write a summary feather file for R containing
	(1) dNds score, error, number of residues, expected dNdS if constrained only by burrial
	(2) ConSurf rate4site score, error, number of residues
	(3) ConSurf rate4site score, error, number of residues
	for each bin of the structural property
	Args:
		inDir (str): 	path to the folder with dNds results
		name (str):		name of the extension to add to the feather file ('' if standard analysis, '_AdditionalDistantSpecies' if adding species for supplementary)
		plotDir (str):	path to the plotting directory
		varName(str): 	name of the structural property
	'''
	path = os.path.join(plotDir,varName+'.feather')
	if not os.path.exists(path):
		# 1) dN/dS:
		dnds_dir = os.path.join(inDir, varName)	
		dnds_df = read_dnds_values(dnds_dir, varName) 
		
		# 2) ConSurf
		rsaPath = os.path.join(RSA_DIR, 'averages_for_'+varName+'_bins.csv')
		columns = ["bin_name","consurfDB_score","consurfDB_std","consurfDB_resNum","consurfRate4Site_score","consurfRate4Site_std","consurfRate4Site_resNum"]
		consurf_df = pd.read_csv(rsaPath, header=0,usecols=columns)

		# 3) Combined
		df = pd.merge(dnds_df,consurf_df, how='left',on=["bin_name"])

		# 4) Add expected fit of constrained only by non-interfacial burrial 
		if addExpected:
			expected = get_expected_dnds(RSA_DIR,rsaPath,varName,inDir)
			df['dnds_expected'] = expected 

		# 5) Write the feather file
		feather.write_dataframe(df, path)
	else:
		df = pd.read_feather(path)
	
	return df 

def main():

	# 1) Set up
	path_closely_related = '../data/processed/modeled_ppis.pkl'
	path_additional_closely_related = '../data/processed/modeled_ppis_additional_closely_related.pkl'
	path_additional_distantly_related = '../data/processed/modeled_ppis_additional_distantly_related.pkl'

	name_closely_related= ""
	name_additional_closely_related ="_AdditionalCloseSpecies"
	name_additional_distantly_related ="_AdditionalDistantSpecies"

	# 2) Select the model we are running calculations for
	# modelPath = path_additional_distantly_related
	# name = name_additional_distantly_related
	modelPath = path_closely_related
	name = name_closely_related
	
	# # Run binning once per model
	# for name,modelPath in zip(names,[path_closely_related,path_additional_closely_related,path_additional_distantly_related]):
	# Run binning once per model
	for name,modelPath in zip([name],[modelPath]):

		plotDir = '../plots'+name+'/dataFromPython'
		if not os.path.exists(plotDir):
			os.makedirs(plotDir)
		dataDir = '../data/evolutionary_rates'+name
		inDir = os.path.join(dataDir, 'dnds')
		RSA_DIR = os.path.join(dataDir, 'RSA_distributions')
		FASTA_DIR = os.path.join(dataDir, 'binned_codons')
	
	# Un-binned data (structural properties, ConSurf DB and ConSurf rate4site only - for all residues)

		# 0) Load data and write feather file for ConSurf (rate4site and DB) in the model (Figure 3, Figure 4,Figure 5)
		df_unbinned = load_consurf(modelPath,name,plotDir)
		
	# Binned data dN/dS + average consurfs / RSA for bins of structural properties
		
		# 1) interfacial vs non-interfacial bins (Figure 3) /!\ NOT DONE RUNNING!
		df_interf_nnInterf_bins = load_interfacial_vs_nonInterfacial(inDir,name,plotDir,RSA_DIR,addExpected=True)
		
		# 2) complex RSA interfacial bins (Figure 3, Figure 4)
		path= os.path.join(inDir,"rsa_complex")
		df_complexRSA_interfacial_bins = load_structural_prop_bins(path,name,plotDir,'rsa_complex_interfacial',RSA_DIR,addExpected=True)
		
		# 3) complex RSA non-interfacial bins /!\ NOT DONE RUNNING! (Figure 3)
		path= os.path.join(inDir,"rsa_complex")
		df_complexRSA_non_interfacial_bins = load_structural_prop_bins(path,name,plotDir,'rsa_complex_non_interfacial',RSA_DIR,addExpected=True)
		
		# 4) monomer RSA interfacial bins (Figure 4)
		path= os.path.join(inDir,"rsa_monomer")
		df_monomerRSA_interfacial_bins = load_structural_prop_bins(path,name,plotDir,'rsa_monomer_interfacial',RSA_DIR,addExpected=True)
		
		# 5) deltaRSA (Figure 6) /!\ NOT DONE RUNNING!
		path= inDir
		df_drsa_bins = load_structural_prop_bins(path,name,plotDir,'drsa',RSA_DIR,addExpected=True)

		# 6) interRRC (Figure 6)
		path= inDir
		df_interRRC_bins = load_structural_prop_bins(path,name,plotDir,'residueContacts',RSA_DIR,addExpected=True)
		
		# 7) dCenter (Figure 6)
		path= inDir
		df_dCenter_bins = load_structural_prop_bins(path,name,plotDir,'distance_from_center',RSA_DIR,addExpected=True)
		
		# 8) dEdges (Figure 6)
		path= inDir
		df_dEdges_bins = load_structural_prop_bins(path,name,plotDir,'distance_from_edges',RSA_DIR,addExpected=True)
		

		exit() # TODO: only remove when we have run more distantly related species.
	
	
	return


if __name__ == '__main__':
	main()








