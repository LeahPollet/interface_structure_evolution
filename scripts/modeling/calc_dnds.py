"""
Script to run run_codeml and calculate dN/dS scores for the codon bins listed in '../../data/evolutionary_rates/binned_codons/dir_list.txt
Relationship and tree file for the species of interest to the dN/dS caluclations are required

File outputs:
	- many fasta files for each binning of the codons
	- many files containing average RSA values in different bins
	- csv files of protein and interface properties

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os, shutil, sys, subprocess
import glob
import time
from tqdm import tqdm
import tarfile

relationship = "(((((Scer,Spar),(Smik,Sbay)),Ncas),Cgla),((Agos,Klac),Calb));" # Relationship between the 3 yeast species
tree = '../../data/external/yeasts.tree'


def download_codeml():
	
	oPath = "../data"
	url = "http://abacus.gene.ucl.ac.uk/software/paml4.9j.tgz"
	# useful if not compiling properly: http://abacus.gene.ucl.ac.uk/software/paml.html

	zipped = oPath+"/"+os.path.basename(url)
	path = oPath+"/"+os.path.basename(url)[:8]
	executable = oPath+"/"+os.path.basename(url)[:8]+"/bin/codeml"
	if not os.path.exists(executable):
		print("Downloading codeml")
		# Download and unzip
		os.system('wget -P '+ oPath +' '+url)
		tar = tarfile.open(zipped)
		tar.extractall(oPath)
		tar.close()
		os.remove(zipped)
		# Compile and move to bin
		args = ('rm', '-rf', path+'/bin/*.exe')
		subprocess.call('%s %s %s' % args, shell=True)
		source_files = path+"/src"
		subprocess.call('make', shell=True,cwd=source_files)
		subprocess.call('mv baseml basemlg codeml pamp evolver yn00 chi2 ../bin', shell=True,cwd=source_files)

	absolute_path = os.path.abspath(executable)
	
	return absolute_path

def read_dir_list(fPath, fastaDir):
	"""Read text file containing list of dirs to run over.

	Args:
		fPath (str): path to file listing categories to run over.
		fastaDir (str): directory containing sequence files.

	Returns:
		list(str): directories containing the categories to run over.

	"""
	f = open(fPath, 'r')
	dirs = []
	for l in f.read().splitlines():
		dirs.append(l)
	f.close()
	return dirs

def ignore_files(dir, files):
	""" Ignore function, to copy a directory structure without the files in it 
	https://stackoverflow.com/questions/15663695/shutil-copytree-without-files """
	return [f for f in files if os.path.isfile(os.path.join(dir, f))]

def main():

	"""Submit dN/dS calculation jobs to computing cluster.

	Directories to submit jobs for are listed in a separate text file.
	A single job is submitted for each fasta file.

	"""
	# 1) Download codeml
	path_to_codeml = download_codeml()
	print(path_to_codeml)

	# 2) Set up
	path_closely_related = '../data/processed/modeled_ppis.pkl'
	# path_additional_closely_related = '../data/processed/modeled_ppis_additional_closely_related.pkl'
	path_additional_distantly_related = '../data/processed/modeled_ppis_additional_distantly_related.pkl'

	name_closely_related= ""
	# name_additional_closely_related ="_AdditionalCloseSpecies"
	name_additional_distantly_related ="_AdditionalDistantSpecies"

	# 3) Select the model we are running calculations for
	modelPath = path_additional_distantly_related
	name = name_additional_distantly_related
	
	# 4) Set up
	tree = "../data/external/yeasts"+name+".tree"
	OUT_DIR = '../data/evolutionary_rates'+name
	FASTA_DIR = os.path.join(OUT_DIR, 'binned_codons')
	DNDS_DIR = os.path.join(OUT_DIR, 'dnds')
	if not os.path.exists(DNDS_DIR): # Copy full directory structure of the fasta dir to the dnds directory, but without the files!
		shutil.copytree(FASTA_DIR,DNDS_DIR, ignore=ignore_files)

	# Select the folder we are running calculations for
	# When running for additional models:
	# only run name_additional_distantly_related, only run with low boot number
	nBootstraps = '5'
	dirs = read_dir_list(os.path.join(FASTA_DIR,'dir_list8.txt'), FASTA_DIR)
	nJobs = 0

	for dir in dirs:
		files = os.listdir(dir)
		subDir = dir.replace(FASTA_DIR + '/', '')
		print ("---------------------------")
		print (">>>", subDir)
		print ("---------------------------")
		for f in files:
			print ("\n")
			print ("	>>>", f)
			print ("	---------------------------")
			

			# Path to the fasta file
			fastaPath = os.path.join(dir, f)
			# Path to the corresponding dnds file:
			dndsPath = (fastaPath.replace('binned_codons' ,'dnds')).replace('.fasta', '.dnds')
			
			# check input not empty
			if os.stat(fastaPath).st_size == 0:
				print ('	no codons in input file: ', fastaPath)
				continue

			if f == "README.txt":
				print ("	Readme file, Skipping")
				continue

			if f[-3:] == "png":
				print ("	Graph file, Skipping")
				continue

			if os.path.exists(dndsPath):
				print ("	file " + dndsPath + " exists. Skipping")
				continue
			else:
				nBs = nBootstraps
			
			name = '_'.join([f.replace('.fasta', ''), subDir.replace('/', '_')])

			cmd = 'python run_codeml.py'
			# cmd = cmd + ' -c ' + PATH_TO_CODEML
			cmd = cmd + ' -a ' + fastaPath
			cmd = cmd + ' -t ' + tree
			cmd = cmd + ' -o ' + dndsPath
			cmd = cmd + ' -b ' + nBs
			cmd = cmd + ' -w ' + (fastaPath.replace('binned_codons' ,'dnds'))
			cmd = cmd + ' -p ' + path_to_codeml
			os.system(cmd)
			nJobs += 1
			
			

	print ('submited ', nJobs, ' jobs')

	return
	
if __name__ == '__main__':
	main()
