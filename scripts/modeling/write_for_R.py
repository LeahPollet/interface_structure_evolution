"""
Script to read and format results from the analysis for plotting

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os,re, sys
import feather
import pickle
import pandas as pd
from modeled_proteins import Protein, Residue
from tqdm import tqdm
import csv
from scipy import optimize as opt
import numpy as np
from collections import Counter


modeled_proteins = '../../data/processed/modeled_proteins.pkl'


plotDir = '../../plots/dataFromPython'
dataDir = '../../data/evolutionary_rates/modeled_proteins.pkl'
inDir = os.path.join(dataDir, 'dnds')
RSA_DIR = os.path.join(dataDir, 'RSA_distributions')
FASTA_DIR = os.path.join(dataDir, 'binned_codons')

	
def read_dnds_values_random_bins(inDir,avg_struct_prop_file):
	"""Reads the dN/dS values and average value of each structural property from a random binning of residues
	Args:
		inDir (str): directory with the files containing the dN/dS values
		avg_struct_prop_file (str): File with the average values of structural properties for each bin

	Returns:
		df (pandas dataframe): dataframe summarizing results
	"""
	# READ DNDS RESULTS
	df = pd.DataFrame(columns=['bin_num', 'num_nt',
							   'dnds', 'dnds_err',
							   'dn', 'dn_err',
							   'ds', 'ds_err'])
	for resultsName in os.listdir(inDir):

		resultsPath = os.path.join(inDir, resultsName)
		if os.path.isdir(resultsPath): 			# skip the directories
			continue

		r = {}

		r['bin_num'] = float(resultsName.replace('.dnds', ''))
		

		with open(resultsPath, 'r') as f:
			lines = f.readlines()
			r['num_nt'] = int(lines[0][10:])
			r['dnds'] = float(lines[4].split()[1])
			r['dnds_err'] = float(lines[4].split()[4])
			r['dn'] = float(lines[2].split()[1])
			r['dn_err'] = float(lines[2].split()[4])
			r['ds'] = float(lines[3].split()[1])
			r['ds_err'] = float(lines[3].split()[4])

		df = df.append(r, ignore_index=True)

		df = df.sort_values('bin_num')
		df.reset_index(drop=True, inplace=True) 
	
	# READ STRUCTURAL PROPERTIES
	prop = pd.read_csv(avg_struct_prop_file,skiprows=[0], header=0,delim_whitespace=True)
	prop = prop.drop(['bin_number', 'num_of_residues'],axis=1)
	df = pd.concat([df.reset_index(drop=True), prop], axis=1)
	
	return df

def read_dnds_values(inDir, varName, addAvrgRSA=False):
	"""Read dnds values from 1D binning of a variable.

	Args:
		inDir (str): directory with the files containing the dN/dS values
		varName (str): the name of the variable that is used to bin,
					   exactly as appears in the file names

	Returns:
		df (pandas dataframe): dataframe summarizing results
	"""
	df = pd.DataFrame(columns=['bin_low', 'bin_high', 'num_nt',
							   'dnds', 'dnds_err',
							   'dn', 'dn_err',
							   'ds', 'ds_err'])
	for resultsName in os.listdir(inDir):

		resultsPath = os.path.join(inDir, resultsName)
		if os.path.isdir(resultsPath): 			# skip the directories
			continue

		if resultsName.startswith(varName):

			if resultsName.replace('.dnds', '').split('_')[-2] != 'to':
				raise UserWarning(resultsName +' does not follow expected file name spec')
			r = {}

			r['bin_high'] = float(resultsName.replace('.dnds', '').split('_')[-1])
			r['bin_low'] = float(resultsName.replace('.dnds', '').split('_')[-3])

			with open(resultsPath, 'r') as f:
				lines = f.readlines()
				r['num_nt'] = int(lines[0][10:])
				r['dnds'] = float(lines[4].split()[1])
				r['dnds_err'] = float(lines[4].split()[4])
				r['dn'] = float(lines[2].split()[1])
				r['dn_err'] = float(lines[2].split()[4])
				r['ds'] = float(lines[3].split()[1])
				r['ds_err'] = float(lines[3].split()[4])

			df = df.append(r, ignore_index=True)

		df = df.assign(bin_center=(df['bin_high'] - df['bin_low']) / 2. + df['bin_low'])
		df = df.sort_values('bin_center')
		df.reset_index(drop=True, inplace=True) 

	if addAvrgRSA:
		rsaPath = os.path.join(RSA_DIR, 'rsa_by_' + varName + '.csv')
		df = pd.merge(df,
					  pd.read_csv(rsaPath, header=0),
					  how='left',
					  on=['bin_low', 'bin_high'])
	return df

def read_RSA(RSA_path):
	"""Read average RSA (monomer and complex values) from 1D binning of a variable.

	Args:
		RSA_path (str): path to the file containing average RSAs

	Returns:
		df (pandas dataframe): dataframe summarizing results
	"""
	df = pd.DataFrame(columns=['bin_low', 'bin_high', 'num_res','meanRSAMonomer', 'meanRSAComplex'])
	with open(RSA_path) as file:
		reader = csv.reader(file)
		for line in reader:
			# skip header
			if line[0] =='bin_low':
				continue
			# fill dataframe one row at a time
			row = {}
			row['bin_low'] = float(line[0])
			row['bin_high'] = float(line[1])
			row['num_res'] = float(line[2])
			row['meanRSAMonomer'] = float(line[3])
			row['meanRSAComplex'] = float(line[4])
			df = df.append(row, ignore_index=True)
			df = df.sort_values('bin_low')
			df = df.assign(bin_center=(df['bin_high'] - df['bin_low']) / 2. + df['bin_low'])
	return df

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
	test = np.sqrt(np.diag(covMatrix))
	r = calc_r(x, y, y_err)
	rsq = calc_r_squared(x, y, y_err, slope, yInt)
	chi2, ndof = calc_chi2_old(x, y, y_err, slope, yInt)
	return {'slope': slope,
			'intercept': yInt,
			'slope_error': slopeErr,
			'intercept_error': yIntErr,
			'PCC': r,
			'r_squared': rsq,
			'chi2': chi2 ,
			'ndof': ndof}

def weighted_mean(x, xErr):
	"""Mean of set of values with uncertainties.

	Args:
		x (list(float)): values.
		xErr (list(float)): uncertainties.

	Returns:
		float

	"""
	if len(x) != len(xErr):
		raise UserWarning('Each value should have its own uncertainty.')
	if len(x) == 0:
		raise UserWarning('Empty list of values')
	w = [1 / i for i in xErr]
	return sum([xi * wi for xi, wi in zip(x, w)]) / sum(w)


def calc_r(x, y, y_err):
	"""Pearson correlation coeff, with uncertainty on y points.

	Equation adapted from wikipedia (questionable).

	Args:
		x (list(float)): x values.
		y (list(float)): y values.
		y_err (list(float): uncertinaty on y values.

	Returns:
		float: PCC value taking into account uncertainties.

	"""
	w = [1 / i for i in y_err]
	m_x = np.mean(x)
	m_y = weighted_mean(y, y_err)
	cov_xy = sum([wi * (xi - m_x) * (yi - m_y) for xi, yi, wi in
				  zip(x, y, w)]) / sum(w)
	var_y = sum([wi * (yi - m_y) * (yi - m_y) for yi, wi in
				 zip(y, w)]) / sum(w)
	var_x = np.var(x)
	return cov_xy / np.sqrt(var_y * var_x)

def calc_r_squared(x, y, y_err, slp, intercept):
	"""read_dnds_valuesuared value of data given a linear fit.

	Args:
		x (list(float)): x values.
		y (list(float)): y values.
		y_err (list(float)): uncertainty on y values.
		slp (float): slope of linear fit.
		int (float): y-intercept of linear fit.

	Returns:
		float

	"""
	ssres, sstot = 0., 0.
	ybar = weighted_mean(y, y_err)
	for xi, yi, ei in zip(x, y, y_err):
		sstot = sstot + ((yi - ybar) / ei)**2
		ssres = ssres + ((yi - (slp * xi + intercept)) / ei)**2
	return 1 - ssres / sstot

def calc_chi2_old(x, y, y_err, slope, yInt):
	"""chi2 and degrees of freedom for single-variable linear fit and data."""
	chi2 = sum([(yi - (xi * slope + yInt))**2 / yi_err**2
				for xi, yi, yi_err in zip(x, y, y_err)])
	ndof = len(y) - 2  # assuming 2-parameter linear fit
	return chi2, ndof

def calc_chi2(x, y, yErr, coeffs):
	"""chi2 and number of degrees of freedom for linear regression and data.

	Args:
		x (array): X values.
		y (array): y values.
		yErr (array): uncertainties on y values.
		coeffs(array): linear regression coefficients.

	Returns:
		float: chi^2
		int: ndof

	"""
	# need to add checks on inputs
	chi2 = sum((y - ((x.transpose() * coeffs[:-1]).sum(axis=1) +
			   coeffs[-1]))**2 / yErr**2)
	ndof = len(y) - (len(x) + 1)  # assuming linear fit
	return chi2, ndof

def load_data(inPath):
	""" Loads structural properties, evolutionary state and consurf data from our PPI models
	Writes fether files to summarize the results
	Args: 
		inPath (str): path to the PPI models
	Returns:
		df (pandas dataframe): summary of results

	"""
	path = os.path.join(plotDir,'data_from_models.feather')
	if not os.path.exists(path):
		# Load the data from protein models
		print ("loading protein models")
		with open(inPath, 'rb') as f:
			proteins = pickle.load(f)
		print ("protein models loaded")

		df = pd.DataFrame(columns=['RSA_monomer', 'RSA_complex', 'deltaRSA',
								   'residueContacts', 'distance_from_center',
								   'distance_from_edges', 'consurf_score',
								   'substitution'])
		# Create the data frame
		for p in tqdm(proteins.values()):
			for r in p.residues:
				res = {}
				res['RSA_monomer']			= float(r.rsa)
				res['RSA_complex']			= float(r.rsa-r.dRsa)
				res['deltaRSA']				= float(r.dRsa)
				res['residueContacts']		= float(r.contacts)
				res['distance_from_center']	= float(r.distance_from_center)
				res['distance_from_edges']	= float(r.distance_from_edges)
				res['consurf_score']		= float(r.consurf_score)
				res['substitution']			= float(r.substitution)

				df = df.append(res, ignore_index=True)
		
		# Save the dataframe
		feather.write_dataframe(df, path)

	# Feather file already exists
	else:
		df = pd.read_feather(path)

	return (df)

def load_RSA_dnds_interfacial():
	""" Loads dN/dS values for the binning of interfacial residues by their values of RSA (monomer and complex)
	Writes fether files to summarize the results
	
	Returns:
		df1/df2 (pandas dataframes): summary of results for the monomer and complex RSA binnings respectively

	"""
	# 1) Monomer RSA for interfacial residues 
	monomer_RSA_interfacial = os.path.join(inDir, 'rsa_monomer/rsa_monomer_interfacial')
	df1 = read_dnds_values(monomer_RSA_interfacial, 'rsa_monomer_interfacial') 
	path = os.path.join(plotDir,'monomerRSA_interfacial.feather')
	# Drop the last row (error bars are too large)
	df1.drop(df1.tail(1).index,inplace=True)
	feather.write_dataframe(df1, path)

	# 2) Complex RSA, for interfacial residues
	complex_RSA_interfacial = os.path.join(inDir, 'rsa_complex/rsa_complex_interfacial')
	df2 = read_dnds_values(complex_RSA_interfacial, 'rsa_complex_interfacial') 
	path = os.path.join(plotDir,'complexRSA_interfacial.feather')
	# Drop the last row (error bars are too large)
	df2.drop(df2.tail(1).index,inplace=True)
	feather.write_dataframe(df2, path)

	return df1,df2

def load_RSA_dnds():
	""" Loads dN/dS values for the binning of all residues by their values of RSA (monomer and complex)
	Writes fether files to summarize the results
	
	Returns:
		df1/df2 (pandas dataframes): summary of results for the monomer and complex RSA binnings respectively

	"""
	# RSA complex all, RSA monomer all, RSA complex interfacial, RSA complex non-interfacial
	# 1) Monomer RSA for all residues 
	monomer_RSA_all = os.path.join(inDir, 'rsa_monomer/rsa_monomer_all')
	df1 = read_dnds_values(monomer_RSA_all, 'rsa_monomer_all') 
	path = os.path.join(plotDir,'monomerRSA.feather')
	feather.write_dataframe(df1, path)

	# 2) Complex RSA, for all residues
	complex_RSA_all = os.path.join(inDir, 'rsa_complex/rsa_complex_all')
	df2 = read_dnds_values(complex_RSA_all, 'rsa_complex_all') 
	path = os.path.join(plotDir,'complexRSA.feather')
	feather.write_dataframe(df2, path)

	# 3) Complex RSA for interfacial residues 
	complex_RSA_interfacial = os.path.join(inDir, 'rsa_complex/rsa_complex_interfacial')
	df3 = read_dnds_values(complex_RSA_interfacial, 'rsa_complex_interfacial') 
	path = os.path.join(plotDir,'complexRSA_interfacial_full.feather')
	# Drop the last row (error bars are too large)
	df3.drop(df3.tail(1).index,inplace=True)
	feather.write_dataframe(df3, path)

	# 4) Complex RSA, for non-interfacial residues
	complex_RSA_non_interfacial = os.path.join(inDir, 'rsa_complex/rsa_complex_non_interfacial')
	df4 = read_dnds_values(complex_RSA_non_interfacial, 'rsa_complex_non_interfacial') 
	path = os.path.join(plotDir,'complexRSA_non_interfacial.feather')
	# Drop the last row (error bars are too large)
	df4.drop(df4.tail(1).index,inplace=True)
	feather.write_dataframe(df4, path)
	
	return df1,df2

def load_interf_inv_dnds():
	""" Loads dN/dS values for the binning of all residues by their values of interface involvments (dRSA, residue contacts, dCenter, dEdges)
	Writes fether files to summarize the results
	
	Returns:
		df1/df2/df3/df4 (pandas dataframes): summary of results for the monomer and complex RSA binnings respectively

	"""
	# 1) Delta RSA
	dRSA_dir = os.path.join(inDir, 'drsa')	
	df1 = read_dnds_values(dRSA_dir, 'drsa') 
	# Get expected dnds:
	RSA_distribution_file = 'rsa_by_drsa.csv'
	expected_dnds_all,expected_dnds_interfacial,expected_dnds_non_interfacial = get_expected_dnds(RSA_distribution_file)
	# Add to the dataframe
	df1['dnds_expected_all'] = expected_dnds_all 
	df1['dnds_expected_interfacial'] = expected_dnds_interfacial 
	df1['dnds_expected_non_interfacial'] = expected_dnds_non_interfacial 
	path = os.path.join(plotDir,'deltaRSA.feather')
	feather.write_dataframe(df1, path)
	
	# 2) Residue contacts
	residueContact_dir = os.path.join(inDir, 'residueContacts')	
	df2 = read_dnds_values(residueContact_dir, 'residueContacts') 
	# Get expected dnds:
	RSA_distribution_file = 'rsa_by_residueContacts.csv'
	expected_dnds_all,expected_dnds_interfacial,expected_dnds_non_interfacial = get_expected_dnds(RSA_distribution_file)
	# Add to the dataframe
	df2['dnds_expected_all'] = expected_dnds_all 
	df2['dnds_expected_interfacial'] = expected_dnds_interfacial 
	df2['dnds_expected_non_interfacial'] = expected_dnds_non_interfacial 
	path = os.path.join(plotDir,'residueContacts.feather')
	feather.write_dataframe(df2, path)

	# 3) Distance from center
	distFromCenter_dir = os.path.join(inDir, 'distance_from_center')	
	df3 = read_dnds_values(distFromCenter_dir, 'distance_from_center')
	# Get expected dnds:
	RSA_distribution_file = 'rsa_by_distance_from_center.csv'
	expected_dnds_all,expected_dnds_interfacial,expected_dnds_non_interfacial = get_expected_dnds(RSA_distribution_file)
	# Add to the dataframe
	df3['dnds_expected_all'] = expected_dnds_all 
	df3['dnds_expected_interfacial'] = expected_dnds_interfacial 
	df3['dnds_expected_non_interfacial'] = expected_dnds_non_interfacial 
	# Drop the two last rows error bars are too large
	df3.drop(df3.tail(2).index,inplace=True) 
	path = os.path.join(plotDir,'distance_from_center.feather')
	feather.write_dataframe(df3, path)

	# 1) Distance from edges
	distFromEdges_dir = os.path.join(inDir, 'distance_from_edges')	
	df4 = read_dnds_values(distFromEdges_dir, 'distance_from_edges') 
	# Get expected dnds:
	RSA_distribution_file = 'rsa_by_distance_from_edges.csv'
	expected_dnds_all,expected_dnds_interfacial,expected_dnds_non_interfacial = get_expected_dnds(RSA_distribution_file)
	# Add to the dataframe
	df4['dnds_expected_all'] = expected_dnds_all 
	df4['dnds_expected_interfacial'] = expected_dnds_interfacial 
	df4['dnds_expected_non_interfacial'] = expected_dnds_non_interfacial  
	# Drop the last row error bars are too large
	df4.drop(df4.tail(2).index,inplace=True)
	path = os.path.join(plotDir,'distance_from_edges.feather')
	feather.write_dataframe(df4, path)
	return df1,df2,df3,df4

def get_expected_dnds(path_to_file):
	""" Get the expected value of dN/dS for a bin of residues if constrained only by burrial (general burrial, interfacial burrial, non-interfacial burrial)
	Returns:
		expected values according to the 3 possible fits (general burrial, interfacial burrial, non-interfacial burrial)

	"""
	# Get the average RSA value of each bin
	RSA_distribution_path = os.path.join(RSA_DIR,path_to_file)
	df = read_RSA(RSA_distribution_path) # Average RSA (monomer and complex per bin) 
	RSAs_complex = df['meanRSAComplex'].tolist()

	# Get expected fit if constrainted only by burrial (complex RSA)
	# Get the RSA-dn/ds relationship for all residues: 
	df_expected_all = read_dnds_values(inDir+"/rsa_complex/rsa_complex_all", 'rsa_complex_all')
	binCenters_all, yVals_all, yErr_all = df_expected_all['bin_center'].tolist(), df_expected_all['dnds'].tolist(), df_expected_all['dnds_err'].tolist()
	fit_all = linear_fit(binCenters_all, yVals_all, yErr_all)
	# Get the RSA-dn/ds relationship for interfacial residues: 
	df_expected_interfacial = read_dnds_values(inDir+"/rsa_complex/rsa_complex_interfacial", 'rsa_complex_interfacial')
	binCenters_interfacial, yVals_interfacial, yErr_interfacial = df_expected_interfacial['bin_center'].tolist(), df_expected_interfacial['dnds'].tolist(), df_expected_interfacial['dnds_err'].tolist()
	fit_interfacial = linear_fit(binCenters_interfacial, yVals_interfacial, yErr_interfacial)
	# Get the RSA-dn/ds relationship for non interfacial residues: 
	df_expected_non_interfacial = read_dnds_values(inDir+"/rsa_complex/rsa_complex_non_interfacial", 'rsa_complex_non_interfacial')
	binCenters_non_interfacial, yVals_non_interfacial, yErr_non_interfacial = df_expected_non_interfacial['bin_center'].tolist(), df_expected_non_interfacial['dnds'].tolist(), df_expected_non_interfacial['dnds_err'].tolist()
	fit_non_interfacial = linear_fit(binCenters_non_interfacial, yVals_non_interfacial, yErr_non_interfacial)

	# Compute the expected values from the fits
	expected_dnds_all = [fit_all['slope'] * x + fit_all['intercept'] for x in RSAs_complex]
	expected_dnds_interfacial = [fit_interfacial['slope'] * x + fit_interfacial['intercept'] for x in RSAs_complex]
	expected_dnds_non_interfacial = [fit_non_interfacial['slope'] * x + fit_non_interfacial['intercept'] for x in RSAs_complex]
	
	return expected_dnds_all,expected_dnds_interfacial,expected_dnds_non_interfacial

def get_python_linear_fits_interface():
	# For interfacial residues only
	df = pd.DataFrame(columns=['slope', 'intercept', 'slope_error',
							   'intercept_error', 'PCC',
							   'r_squared', 'chi2',
							   'ndof'],index= ["rsa_monomer","rsa_complex","delta_rsa","residue_contacts","dCenter","dEdges"])
	# 1) Monomer RSA
	df_temp = read_dnds_values(inDir+"/rsa_monomer/rsa_monomer_interfacial", 'rsa_monomer_interfacial')
	# Drop the last row error bars are too large
	df_temp.drop(df_temp.tail(1).index,inplace=True)
	binCenters, yVals, yErr = df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_monomer'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	# 2) Complex RSA
	df_temp = read_dnds_values(inDir+"/rsa_complex/rsa_complex_interfacial", 'rsa_complex_interfacial')
	# Drop the last row error bars are too large
	df_temp.drop(df_temp.tail(1).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_complex'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	# 3) Delta RSA
	df_temp = read_dnds_values(inDir+"/drsa", 'drsa')
	# Drop the 0 bin, to get only interfacial residues
	df_temp.drop(df_temp.head(1).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['delta_rsa'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	# 4) Residue contacts
	df_temp = read_dnds_values(inDir+"/residueContacts", 'residueContacts')
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['residue_contacts'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	# 5) dCenter
	df_temp = read_dnds_values(inDir+"/distance_from_center", 'distance_from_center')
	# Drop the two last rows error bars are too large
	df_temp.drop(df_temp.tail(2).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['dCenter'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	# 6) dEdges
	df_temp = read_dnds_values(inDir+"/distance_from_edges", 'distance_from_edges')
	# Drop the last row error bars are too large
	df_temp.drop(df_temp.tail(2).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['dEdges'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	
	path = os.path.join(plotDir,'linearFits_interfaces.feather')
	feather.write_dataframe(df, path)
	
	return df

def get_python_linear_fits():
	# For RSA monomer all, RSA complex all, RSA complex interfacial, RSA complex non interfacial
	df = pd.DataFrame(columns=['slope', 'intercept', 'slope_error',
							   'intercept_error', 'PCC',
							   'r_squared', 'chi2',
							   'ndof'],index= ["rsa_monomer_all","rsa_complex_all","rsa_complex_interfacial","rsa_complex_non_interfacial"])

	# 1) Monomer RSA all
	df_temp = read_dnds_values(inDir+"/rsa_monomer/rsa_monomer_all", 'rsa_monomer_all')
	binCenters, yVals, yErr = df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_monomer_all'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	
	# 2) Complex RSA all
	df_temp = read_dnds_values(inDir+"/rsa_complex/rsa_complex_all", 'rsa_complex_all')
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_complex_all'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	
	# 2) Complex RSA interfacial
	df_temp = read_dnds_values(inDir+"/rsa_complex/rsa_complex_interfacial", 'rsa_complex_interfacial')
	# Drop the last row error bars are too large
	df_temp.drop(df_temp.tail(1).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_complex_interfacial'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	
	# 2) Complex RSA non interfacial
	df_temp = read_dnds_values(inDir+"/rsa_complex/rsa_complex_non_interfacial", 'rsa_complex_non_interfacial')
	# Drop the last row error bars are too large
	df_temp.drop(df_temp.tail(1).index,inplace=True)
	binCenters, yVals, yErr =  df_temp['bin_center'].tolist(), df_temp['dnds'].tolist(), df_temp['dnds_err'].tolist()
	fit = linear_fit(binCenters, yVals, yErr)
	df.loc['rsa_complex_non_interfacial'] = [fit["slope"],fit["intercept"],fit["slope_error"],fit["intercept_error"],fit["PCC"],fit["r_squared"],fit["chi2"],fit["ndof"]]
	pd.set_option('display.max_columns', 500)
	
	path = os.path.join(plotDir,'linearFits_RSAs.feather')
	feather.write_dataframe(df, path)
	return df

def load_data_multiple_linear_regression():
	"""Load data for multiple liear regression """
	# 1) 100 random bins
	random_bins_100_dnds = os.path.join(inDir, '100_random_bins')
	random_bins_100_struct_prop = os.path.join(FASTA_DIR, '100_random_bins','README.txt')
	df1 = read_dnds_values_random_bins(random_bins_100_dnds, random_bins_100_struct_prop) 
	path = os.path.join(plotDir,'random_bins_100.feather')
	feather.write_dataframe(df1, path)

	return df1

def main():

	if not os.path.exists(plotDir):
		os.makedirs(plotDir)

	inPath = modeled_proteins
	# # # 1) Load data and write feather file for substitution and consurf score
	df = load_data(inPath)
	
	# 2) load data and write feather files to study the relationship between RSAs (monomer and complex) and dnds for INTERFACIAL residues
	df1,df2 = load_RSA_dnds_interfacial()
	
	# # # 3) Load data and write feather files to study the relationship between RSAs (monomer and complex) and dnds for ALL residues
	# # Fig 4
	df3,df4 = load_RSA_dnds()
	
	# # 3) load data and write feather files to study the relationship between interface involvment and dnds for interfacial residues
	df5,df6,df7,df8, = load_interf_inv_dnds()
	
	# # # 4) Get the python linear fits as a feather file to compare to work in R for INTERFACES
	df_fits_interf = get_python_linear_fits_interface()

	# # 5) Get the python linear fits as a feather file to compare to work in R RSA monomer all vs RSA complex all, and RSA complex interfacial vs RSA complex non-interfacial (Fig 4)
	df_fits = get_python_linear_fits()
	
	# # 6) Load data for multiple linear regression with random bins
	df_mult_reg_100 = load_data_multiple_linear_regression()
	return


if __name__ == '__main__':
	main()




