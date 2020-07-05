"""
	Step 1 in the pipeline
	
	Sets up the initial data folder structure needed for the analysis (external and processed data folders) 
	Automated download of the external data from Biogrid, PDB and Ensembl. Need to: Update the ftp dowload links when new data is made available 

	Requirements: 
				tqdm==4.42.0
				bio==0.1.0

Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os
import stat
import sys
import re
import gzip
import zipfile
import subprocess
import csv
import pickle
from collections import defaultdict, OrderedDict,Counter
from tqdm import tqdm
from Bio import SeqIO
from Bio.PDB import PDBParser, PDBIO, Select,PDBList
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from statistics import mean 
import tarfile


def download_initial_data(dir):
	"""
	Get initial input files from SGD, PDB, BioGrid, Ensembl.
	Args:
		dir (str): directory to download them into (/data/external directory)

	"""
	print ("----------- Download initial data -----------")
	remoteFiles = [
	('http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.5.180/BIOGRID-ALL-3.5.180.tab2.zip',"bioGrid_all.txt"),				# Using December 2019 data 
	('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt','pdb_entry_type.txt'),
	('ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_seqres.txt.gz','pdb_seqres.txt'),
	('ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/resolu.idx','resolu.idx'),
	('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz','pdb_chain_uniprot.csv'),								# Download manually if issues with the REST API (https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html)
	('ftp://ftp.ensembl.org/pub/release-95/fasta/saccharomyces_cerevisiae/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz','Scer.nt'),		# S. cere
	('ftp://ftp.ensembl.org/pub/release-95/fasta/saccharomyces_cerevisiae/pep/Saccharomyces_cerevisiae.R64-1-1.pep.all.fa.gz','Scer.aa'),
	('https://downloads.yeastgenome.org/sequence/fungi/S_uvarum/archive/MIT/orf_dna/orf_genomic.fasta.gz','Sbay.nt'),							# S. bay
	('https://downloads.yeastgenome.org/sequence/fungi/S_uvarum/archive/MIT/orf_protein/orf_trans.fasta.gz','Sbay.aa'),
	('https://downloads.yeastgenome.org/sequence/fungi/S_mikatae/archive/MIT/orf_dna/orf_genomic.fasta.gz','Smik.nt'),							# S. mik
	('https://downloads.yeastgenome.org/sequence/fungi/S_mikatae/archive/MIT/orf_protein/archive/orf_trans.20041119.fasta.gz','Smik.aa'),
	('https://downloads.yeastgenome.org/sequence/fungi/S_paradoxus/archive/MIT/orf_dna/orf_genomic.fasta.gz','Spar.nt'),						# S. par
	('https://downloads.yeastgenome.org/sequence/fungi/S_paradoxus/archive/MIT/orf_protein/orf_trans.fasta.gz','Spar.aa'),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_castellii_cbs_4309_gca_000237345/cds/Naumovozyma_castellii_cbs_4309_gca_000237345.ASM23734v1.cds.all.fa.gz',"Ncas.nt"), 	# N. Cas
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_castellii_cbs_4309_gca_000237345/pep/Naumovozyma_castellii_cbs_4309_gca_000237345.ASM23734v1.pep.all.fa.gz',"Ncas.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/_candida_glabrata_gca_001466525/cds/_candida_glabrata_gca_001466525.ASM146652v1.cds.all.fa.gz',"Cgla.nt"), 	#C. Glab
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/_candida_glabrata_gca_001466525/pep/_candida_glabrata_gca_001466525.ASM146652v1.pep.all.fa.gz',"Cgla.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/ashbya_gossypii/cds/Ashbya_gossypii.ASM9102v1.cds.all.fa.gz',"Agos.nt"), 	#A. Gos
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/ashbya_gossypii/pep/Ashbya_gossypii.ASM9102v1.pep.all.fa.gz',"Agos.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/cds/Kluyveromyces_lactis_gca_000002515.ASM251v1.cds.all.fa.gz',"Klac.nt"), #K. Lact
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/kluyveromyces_lactis_gca_000002515/pep/Kluyveromyces_lactis_gca_000002515.ASM251v1.pep.all.fa.gz',"Klac.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota2_collection/candida_albicans_sc5314_gca_000784635/cds/Candida_albicans_sc5314_gca_000784635.Cand_albi_SC5314_V4.cds.all.fa.gz',"Calb.nt"), # C. Albi
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota2_collection/candida_albicans_sc5314_gca_000784635/pep/Candida_albicans_sc5314_gca_000784635.Cand_albi_SC5314_V4.pep.all.fa.gz',"Calb.aa"),
	]

	for url,shortName in remoteFiles:
		if not os.path.exists(os.path.join(dir, shortName)):
			os.system('wget -P '+dir+' '+url)
			# Rename file to their corresponding shortName, unzip if needed
			outfile = os.path.join(dir,url.split('/')[-1])
			if (url.split('.')[-1] == "zip"):
				with zipfile.ZipFile(outfile, 'r') as zip_ref:
					zip_ref.extractall(dir)
				outfile = outfile[:-4]+".txt"
			if (url.split('.')[-1] == "gz"):
				os.system('gunzip '+outfile)
				outfile = outfile[:-3]
			os.rename(outfile,os.path.join(dir,shortName))

def download_clustal():
	"""
	Downloads and compiles clustal omega for macosx

	Returns: 
		path_to_executable (str): path to the installed ClustalW executable

	Useful if clustal is not compiling properly, or to install on other OS:
	http://www.clustal.org/clustal2/

	"""
	
	oPath = "../../data"
	url = "http://www.clustal.org/omega/clustal-omega-1.2.3-macosx"
	zipped = oPath+"/clustal/"+os.path.basename(url)
	path_to_clustal = oPath+"/clustal"
	path_to_executable = oPath+"/clustal/clustalo"
	if not os.path.exists(path_to_clustal):
		os.makedirs(path_to_clustal)
		print("Downloading clustalW")
		os.system('wget -P '+ path_to_clustal +' '+url)
		os.rename(zipped,path_to_executable)
		st = os.stat(path_to_executable)
		os.chmod(path_to_executable, st.st_mode | stat.S_IEXEC)
	
	return path_to_executable

def format_biogrid(ipath, opath1, opath2):
	"""
	Extract physical interactions between s. cer. proteins from BioGrid.
	Args: 
		ipath (str): path to the Biogrid file as downloaded by download_initial_data
		opath (str): path to the desired location for the filtered output file

	Returns:
		ORFs_in_multiple_reports (set([str,str])): Set of lists of PPIs reported more than once in yeast from the biogrid database

	"""
	print ("----------- Format Biogrid -----------")
	if os.path.exists(opath1):
		print ('Filtered BioGrid file already exists. Skipping.')

		ORF_list = set()
		with open( opath2) as f:
			for items in csv.reader( f, dialect="excel-tab" ):
				ORF_list.update([items[0],items[1]])
				
		return ORF_list

	if not os.path.exists(ipath):
		raise UserWarning('Input file does not exist: '+ipath)
	
	ints = {}
	orfpattern = "Y[A-Z]{2}[0-9]{3}[A-Z]"

	with open( ipath ) as f:
		for items in tqdm(csv.reader( f, dialect="excel-tab" )):
			if len( items ) == 24:
				orf1, orf2 = items[5], items[6]
				method, pubmedid = items[12], items[14]
				tax1, tax2 = items[15], items[16]
				if re.search( orfpattern, orf1 ) and re.search( orfpattern, orf2 ):
					if tax1 == "559292" and tax2 == "559292":
						if method == "physical":
							pair = ( orf1, orf2 ) if orf1 < orf2 else ( orf2, orf1 )
							ints.setdefault( pair, {} )[pubmedid] = 1
	with open(opath1, 'w') as fh1:
		with open(opath2, 'w') as fh2:
			count_multiple_reports = 0
			ORFs_in_all_PPIs, ORFs_in_multiple_reports= set(), set()
			
			for orf1, orf2 in sorted( ints.keys() ):
				ORFs_in_all_PPIs.update([orf1,orf2])
				fh1.write('\t'.join([orf1, orf2, str(len(ints[(orf1, orf2)]))])+'\n')
				if len(ints[(orf1, orf2)]) > 1:
					count_multiple_reports += 1
					ORFs_in_multiple_reports.update([orf1,orf2])
					fh2.write('\t'.join([orf1, orf2, str(len(ints[(orf1, orf2)]))])+'\n')
	# > 111802  PPIs, amongst  5762 proteins downloaded from Biogrid
	print("	>", len(ints.keys())," PPIs, amongst ",len(ORFs_in_all_PPIs),"proteins downloaded from Biogrid")
	# >  17007  PPIs, amongst  4021 proteins are reported more than once	 
	print("	> ", count_multiple_reports," PPIs, amongst ",len(ORFs_in_multiple_reports),"proteins are reported more than once")	
	return ORFs_in_multiple_reports
			
def format_pdb_seqs(ipath, opath):
	"""
	Filter PDB sequence file for proteins.

	Args:
		ipath (str): path to PDB file containing all amino acid sequences
		opath (srt): path to write out filtered file

	"""
	print ("----------- Format PDB seqs -----------")
	if os.path.exists(opath):
		print ('Filtered PDB aa sequence file already exists. Skipping.')
		return
	if not os.path.exists(ipath):
		raise UserWarning('Input file does not exist: '+ipath)
	fh = open(opath, "w")
	for record in tqdm(SeqIO.parse(ipath, "fasta")) :
		if "mol:protein" in record.description:
			record.description = record.name
			SeqIO.write(record, fh, "fasta")
	fh.close()

def get_aa_sequences(ORFList, scerSeqs,opath):
	"""
	Subset the full s. cer. protein aa sequences to only the ORFs of interest to the analyis 

	Args:
		ORFList (set): set of ORF names of interest.
		scerSeqs (str): path to file containing s. cer. protein aa sequences. 
		opath (str): path to output file containing result fasta file.
		
	"""
	print ("----------- Format aa seqs -----------")
	if os.path.exists(opath):
		print ('Filtered aa sequence file already exists. Skipping.')
		return
		
	fh = open(opath, "w")
	for record in SeqIO.parse(scerSeqs, "fasta"):
		if (record.name.split('_')[0]) in ORFList:
			SeqIO.write(record, fh, "fasta")	
	return

def run_blast(databaseFile, queryFile, opath, eValue):
	"""
	Run a blast protein alignment 

	Args:
		databaseFile (str): path to file containing the protein aa sequences used as the target database
		queryFile (str): path to file containing the protein aa sequences used as the query 
		opath (str): path to output file containing results of blast.
		eValue (str): e-value cutoff for results.

	"""
	print ("----------- Run Blast -----------")
	if not os.path.exists(databaseFile+'.pin'):
		os.system('makeblastdb -in ' + databaseFile + ' -dbtype prot')
	if os.path.exists(opath):
		print ('Orthologue matching statistics already exist. Skipping.')
		return
	fmt = '"6 sseqid qseqid slen qlen length sstart send qstart qend nident positive gaps evalue"'
	os.system('blastp -query ' + queryFile +
			  ' -db ' + databaseFile +
			  ' -evalue ' + eValue +
			  ' -outfmt ' + fmt +
			  ' -num_threads 4 > ' + opath +
			  ' 2>/dev/null') # Redirect stderr to /dev/null to supress warnings
	return


def coverage_cutoff_on_BLAST(ipath, opath, coverage,oFormat):
	"""
	Read blast results file and apply coverage cutoff.

	Args:
		ipath (str): path to the file containing the results of blast
		opath (str): path to output file with blast results filtered using the coverage cutoff
		coverage (num): coverage cutoff
		oFormat (str): format wanted for the outfile with target = target name, query = query name, evalue = evalue in the blast results
	"""
	print ("----------- Run coverage cutoff -----------")
	if os.path.exists(opath):
		print ('Coverage cutoff file already exists. Skipping.')
		return
	if not os.path.exists(ipath):
		raise UserWarning('Input file does not exist: '+ipath)
	with open( ipath ) as fh1:
		with open( opath, "w" ) as fh2:
			for items in csv.reader( fh1, dialect="excel-tab" ):

				target = items[0]
				query = items[1]
				evalue = float( items[12] )

				positives = int( items[10] )
				gaps = int( items[11] )
				qlen = float( items[2] )
				slen = float( items[3] )
				qstart = float( items[5] )
				qend = float( items[6] )
				sstart = float( items[7] )
				send = float( items[8] )
				if ( (((qend - qstart) + 1) - gaps) / qlen > coverage ) \
				   and ( (((send - sstart) + 1) - gaps) / slen > coverage ):
					fh2.write(oFormat(target,query,evalue))

		
def get_pdb_structures(rawDir, outDir, homologsPath, pdbMethodsPath):
	"""
	Download, quality control and filter pdb files for all the s. cer. AA sequences

	Args:
		rawDir (str): directory to download pdb files into
		outDir (str): directory to save curated pdb files into
		homologsPath (str): file with PDB structures homology mapped to yeast proteins
		pdbMethodsPath (str): file with experimental method of each PDB structure

	"""
	print ("----------- Get PDB structs -----------")
	if not os.path.exists(rawDir):
		os.makedirs(rawDir)
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	ftpDir = 'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/'
	format = lambda pdbid: pdbid.lower() + ".pdb1"
	# determine prot-only xray structures
	xray = {}
	with open(pdbMethodsPath) as fh: # File the the entry type of each PDB
		for items in csv.reader( fh, dialect='excel-tab' ):
			if 'prot' in items[1] and items[2] == 'diffraction':
				xray[items[0]] = 1
	# determine structures to download
	homologs = {}
	with open(homologsPath) as fh:
		for items in csv.reader( fh, dialect="excel-tab" ):
			homologs[items[1]] = 1
	
	# determine valid structures
	yeastlike = set( homologs.keys() ).__and__( set( xray.keys() ) )
	qcPath = os.path.join(outDir, '../pdb_quality_control_results.csv') # processed/assemblies/../pdb_quality_control_results
	
	if os.path.exists(qcPath):
		print("Download testing and curation already done or partially done.\nSkipping PDB files already downloaded")
		with open(qcPath) as f:
			alreadyDone =[items[0] for items in csv.reader(f)]
	else:
		print ('\nDownloading, checking and filtering PDB files')
		with open(qcPath, 'w') as qcFile:
			qcFile.write('PDB_ID,Download_status\n')
		alreadyDone = []
		
	with open(qcPath, 'a') as qcFile:
		try:
			for pdbID in tqdm(yeastlike):
				if (pdbID not in alreadyDone):
					print ("> ", pdbID)
					ftpPath = os.path.join(ftpDir, format(pdbID)+".gz")
					rawPath = os.path.join(rawDir, format(pdbID))       # external/assemblies/pdbID
					unzippedPath = os.path.join(rawDir,format(pdbID)) 
					outPath = os.path.join(outDir, format(pdbID))       # processed/assemblies/pdbID
					# PDB file not downloaded
					if not os.path.exists(outPath): 
						if not os.path.exists(rawPath):  		# Download and unzip the PDB file              
							os.system('wget -q -P %s %s' % (rawDir, ftpPath))
							os.system('gunzip '+rawPath)
						if not os.path.exists(unzippedPath): 	# Dowload failed error   
							print ('Download failed \n')        
							qcFile.write(pdbID+',failed_download\n')    
						else:     								# Dowload sucessful or curation error                                      
							passed, code = curate_pdb_file(rawPath, outPath)
							if not passed:
								print ('Does not pass curation')
								qcFile.write(pdbID+','+code+'\n')
							else:
								qcFile.write(pdbID+','+"OK"+'\n')
					# PDB file already downloaded (e.g. in previous run of the analysis)
					else:
						qcFile.write(pdbID+','+"OK"+'\n')

		finally:
			with open(qcPath) as f:
				all_PDBs_tested = [items[0] for items in csv.reader( f)]
			with open(qcPath) as f:
				all_PDBs_success = [items[0] for items in csv.reader( f) if items[1] =='OK']
			# print("Done:",len(set(all_PDBs_success)), "/",len(set(all_PDBs_tested)),",PDB files passed the download, testing and curation steps.")
	return

def curate_pdb_file(ipath, opath):
	"""
	Quality check and filter PDB file.

	Args:
		ipath: location of input PDB file.
		opath: location of output filtered PDB file.

	Returns:
		passed (bool): True if structure passed quality checks.
		code (str):  reason for failure of quality check.

	"""
	parser = PDBParser(PERMISSIVE=1, QUIET=1)
	pdbId = os.path.split(ipath)[1].split('.')[0]
	structures = parser.get_structure('name', ipath)
	try:
		passed, code = pdb_quality_tests(structures[0])
	except:
		print (pdbId + ' broke quality control tests')
		raise
	if passed:
		io = PDBIO()
		io.set_structure(structures)
		io.save(opath, Decider())
		# io.save((opath, 'wb'), Decider())
	return passed, code


def pdb_quality_tests(structure):
	"""
	Check if inter CA distances and heavy atom counts are reasonable for a given PDB structure
	
	Args:
		structure: PDB structure
	"""

	# inter CA distance cutoff
	maxAlphaSep = 4.0
	# Cutoff on fraction of bad residues
	maxBadness = 0.05
	# Heavy atoms per residue
	heavyCounts = {'ALA':  5,
				   'ARG': 11,
				   'ASN':  8,
				   'ASP':  8,
				   'CYS':  6,
				   'GLN':  9,
				   'GLU':  9,
				   'GLY':  4,
				   'HIS': 10,
				   'ILE':  8,
				   'LEU':  8,
				   'LYS':  9,
				   'MET':  8,
				   'PHE': 11,
				   'PRO':  7,
				   'SER':  6,
				   'THR':  7,
				   'TRP': 14,
				   'TYR': 12,
				   'VAL':  7}
	distances = []
	missingAtoms = []
	for chain in structure:
		residues = [r for r in chain]
		goodResidues = [r for r in residues if ('CA' in r)
						and (r.get_resname() in heavyCounts)]
		for i in range(len(goodResidues) - 1):
			r1 = goodResidues[i]
			r2 = goodResidues[i+1]
			distances.append(r1['CA'] - r2['CA'])
		for r in goodResidues:
			name = r.get_resname()
			expected = heavyCounts[name]
			observed = len([a for a in r])
			missingAtoms.append(expected - observed)
	if len(distances) == 0:
		return False, 'no_valid_residues'
	badDistances = list(filter((lambda x: x > maxAlphaSep), distances)) # Change for list comprehension?
	if len(badDistances) / float(len(distances)) > maxBadness:
		return False, 'distance'
	badMissingAtoms = list(filter((lambda x: x > 0), missingAtoms))
	if len(badMissingAtoms) / float(len(missingAtoms)) > maxBadness:
		return False, 'counts'
	return True, ''

class Decider (Select):
	"""
	Bio.PDB method for subsetting structures.

	"""

	def accept_model(self, model):
		if model.get_id() == 0:
			return 1
		else:
			return 0

	def accept_residue(self, residue):
		threeLetterAA = set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU',
							 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
							 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'])
		if (residue.get_resname() in threeLetterAA) and ('CA' in residue):
			return 1
		else:
			return 0


def update_homology_map(pdbDir, homologsPath, curatedChains, curatedHomologs):
	"""
	Update the list of PDB proteins homology matched to yeast to include only those that passed the curation tests.
	Args:
			pdbDir (str): 			directory with the downloaded PDB structures
			homologsPath (str): 	file with PDB structures homology mapped to yeast proteins
			curatedChains (str):	file with the list of non-zero length chains from each PDB structure
			curatedHomologs (str)	file with PDB structures that passed curation tests homology mapped to yeast proteins
	"""
	print ("----------- Update homology map -----------")
	if os.path.exists(curatedHomologs):
		print ('Curated Yeast ORF - PDB homologs file already exists. Skipping.')
		return
	okChains = homolog_chains(pdbDir, curatedChains)
	# Intersect curated chains with yeastlike
	with open(homologsPath) as fhIn:
		with open(curatedHomologs, 'w') as fhOut:
			for line in fhIn:
				line = line.strip() # Remove \n
				orf, pdbid, chain = line.split()
				if pdbid in okChains and chain in okChains[pdbid]:
					fhOut.write(line + '\n')
	return

def homolog_chains(pdbDir, outPath):
	"""
	Remove length zero chains
	Args:
		pdbDir (str):  directory with the downloaded PDB structures
		outPath (str): file with the list of non-zero length chains from each PDB structure
	Returns:
		okChains (dict): Chains that passed curation for each PDB ID
	""" 
	print ("----------- Homolog chains -----------")
	okChains = {}
	if os.path.exists(outPath):
		print ('Output file: ' + outPath + ' already exists. Skipping.')
		with open(outPath) as fh:
			for line in fh:
				pdbid, chains = line.strip().split( "\t" )
				okChains[pdbid] = chains
		if len(okChains) == 0:
			raise UserWarning('Empty file: '+outPath)
		return okChains

	pdbFiles = os.listdir(pdbDir)
	with open(outPath, 'w') as fh:
		for i, pdbFile in enumerate(tqdm(pdbFiles)):
			pdbid = pdbFile.split( "." )[0]
			# (PDB file to stdout | Extract only the column with chain information | sort removing duplicates | grep only chains with correct name format | insert a blank line after new line)
			cmd = "cat %s | cut -c22 | sort -u | grep '[A-Za-z0-9]' | perl -pe 's/\n//g'" % (os.path.join(pdbDir, pdbFile))
			p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
			for line in p.stdout:
				fh.write('\t'.join([pdbid, line.decode("utf-8")]) + '\n') #add the decode to cmd?
				okChains[pdbid] = line.decode("utf-8")
	return okChains

def read_yeast_orthologs(externalData):
	"""
	Maps each S.cer ORF to ortholog ORFs in 3 other yeast species

	Returns:
		{str: {str: str}}: s. cer ORF --> species --> ortholog ORF

	"""
	# BLAST S.cer against other yeast species (S.bay, S.mik, S.par)
	print ("----------- Read orthologues in other species -----------")
	eValue = '1e-5'
	for species in ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
		opath= os.path.join(externalData, 'Scer-' + species + '-blast_stats.best.txt')
		if not os.path.exists(opath):
			species1_path = os.path.join(externalData,'Scer.aa')
			print ("> Running BLAST for S.cer against "+species)
			species2_path = os.path.join(externalData,species+'.aa')
			run_blast(species1_path, species2_path, opath, eValue)
		else:
			print ("BLAST of S.cer against ",species ," already done. Skipping.")

	# Apply coverage cutoff to the results
	minimumCoverage = 0.5
	for species in ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']:
		ipath = os.path.join(externalData, 'Scer-'+species+'-blast_stats.best.txt')
		opath = os.path.join(externalData, 'Scer-'+species+'-blast_stats_coverageCutoff.best.txt')
		if not os.path.exists(opath):
			# Format wanted for the output: a the target name and b the query name in the blast results file
			oFormat = lambda a, b, e:('\t'.join([a.split('_')[0]]+[b]+[str(e)])+'\n')		
			coverage_cutoff_on_BLAST(ipath,opath, minimumCoverage,oFormat)
		else:
			print ("Curation of BLAST results already done. Skipping.")
			
	# Read the names of the yeast orthologs.
	# Returns: {str: {str: str}}: s. cer ORF --> species --> orthologS ORF with the best evalue in the BLAST alignment
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

def build_yeast_aligment(externalData,orthologs,opath,path_to_clustalo):
	"""
	Build alignment between S.cer nt, S.cer nt, S.bay nt, S.mik nt and S.par nt
	Args:
		externalData: 					path to the folder containing the raw data
		orthologs ({str: {str: str}}): 	s. cer ORF --> species --> ortholog ORF
		opath (str): 					path to save the codon alignment file

	"""
	print ("----------- Build yeast alignments -----------")
	# Read all of the sequence data
	Scer_ORFs = list(orthologs.keys())
	
	sequences_dict = read_all_sequence_files(externalData)

	# Write a FASTA file per yeast ORF with the sequences (nt) from each strain
	fasta_path = write_fasta_MSA_files(externalData,sequences_dict,Scer_ORFs,orthologs)
	
	# Run a MSA for each ORF (using the fasta file)
	msa_run_clustalw(fasta_path,path_to_clustalo)
	# Combine everything
	combined_codon_alignment_path = os.path.join(externalData,"codon_alignment.txt") # this is the same as opath
	build_combined_file(externalData,opath,fasta_path)
	
	return

def read_all_sequence_files(externalData):
	"""
	Read NT sequence files from all 6 yeast strains into a dictionary
	ORF names formatted to match the format of the orthologues dictionary
	Args: 
		externalData: path to the folder containing the raw data
	Returns:
		sequences_dict = {strain.xx -> {ORF_name -> [seq]}}
	"""
	sequences_dict ={}
	# For each strain, for each ORF update sequences_dict (Sequence dict have all of the ORF sequences for all strains)
	for strain in ['Scer.nt', 'Sbay.nt','Smik.nt', 'Spar.nt','Ncas.nt','Cgla.nt','Agos.nt','Klac.nt','Calb.nt']:
		sequences_dict.update({strain:{}})
		file = SeqIO.parse(open(os.path.join(externalData,strain)),'fasta')
		for fasta in file:
			name, sequence = fasta.id, str(fasta.seq)
			# Modify the ORF name to fit the format in the orthologue dictionary
			if strain[:4] =='Scer':  # strip the trailing "_mRNA"
				name = name[:-5]
			elif (strain[:4]) in ["Sbay","Smik","Spar"]: # For Sbay, Smik, Spar, switch from ORFN:xxxx to ORFP:xxxx
				s= list(name)
				s[3]='P'
				name= ''.join(s)
			else: #Other yeast species as is
				name = name

			seq = [sequence[i:i+3] for i in range(0, len(sequence), 3)][:-1]
			
			sequences_dict[strain].update({name:seq})
	
	return sequences_dict

def write_fasta_MSA_files(externalData,sequences_dict,Scer_ORFs,orthologs):
	"""
	For each S.cer ORF with orthologues in all 3 other yeast species (S.bay, S.mik,S.par),
	create the directory structure and write a fasta files containing the nt sequences of the ORF in all 4 species

	Args:
		externalData: 		path to the folder containing the raw data
		sequences_dict: 	Dicionary with aa and nt sequence info for all strains
		Scer_ORFs:			List of S.cer ORFs
		orthologs: 			Dictionary with mapping between yeast ORFs and orthologue ORFs in the 3 other yeast species
	Returns:
		path to a folder containing a FASTA file with aa sequences in the 4 yeast strains, and a 
		FASTA file with nt sequences in the 4 yeast strains, for each ORF with oththologue in all strains
	"""

	count = 0
	test_ORF = Scer_ORFs[0]
	if not os.path.exists(os.path.join(externalData,"MSA",test_ORF,test_ORF+".nt.fa")):
		
		# For each S.cer ORF
		for ORF in tqdm(Scer_ORFs): 

			# Check if it has an orthologue in the 6 other strains
			if all (k in orthologs[ORF] for k in ('sbay','smik','spar','ncas','cgla','agos','klac','calb')):

				# Create the directory to store the 2 fasta files for the ORF
				file_dir = os.path.join(externalData,"MSA",ORF)
				if not os.path.exists(file_dir):
					os.makedirs(file_dir)

				sbay_ORF = orthologs[ORF]['sbay']
				smik_ORF = orthologs[ORF]['smik']
				spar_ORF = orthologs[ORF]['spar']
				ncas_ORF = orthologs[ORF]['ncas']
				cgla_ORF = orthologs[ORF]['cgla']
				agos_ORF = orthologs[ORF]['agos']
				klac_ORF = orthologs[ORF]['klac']
				calb_ORF = orthologs[ORF]['calb']

				# Get all of the sequences needed for the alignment file
				scer_nt_seq =  sequences_dict['Scer.nt'][ORF]
				sbay_nt_seq =  sequences_dict['Sbay.nt'][sbay_ORF]
				smik_nt_seq =  sequences_dict['Smik.nt'][smik_ORF]
				spar_nt_seq =  sequences_dict['Spar.nt'][spar_ORF]
				ncas_nt_seq =  sequences_dict['Ncas.nt'][ncas_ORF]
				cgla_nt_seq =  sequences_dict['Cgla.nt'][cgla_ORF]
				agos_nt_seq =  sequences_dict['Agos.nt'][agos_ORF]
				klac_nt_seq =  sequences_dict['Klac.nt'][klac_ORF]
				calb_nt_seq =  sequences_dict['Calb.nt'][calb_ORF]

				# A fasta file for each ORF with sequences in all 4 strains
				o_file_path_NT = os.path.join(file_dir,ORF+".nt.fa") 
				
				if not os.path.exists(o_file_path_NT):
					with open(o_file_path_NT, "w" ) as o_file_NT:
						for i,(seq,name) in enumerate([(scer_nt_seq,'Scer.nt'),(sbay_nt_seq,'Sbay.nt'),(smik_nt_seq,'Smik.nt'),(spar_nt_seq,'Spar.nt'),(ncas_nt_seq,'Ncas.nt'),(cgla_nt_seq,'Cgla.nt'),(agos_nt_seq,'Agos.nt'),(klac_nt_seq,'Klac.nt'),(calb_nt_seq,'Calb.nt')]):
							o_file_NT.write(">"+'\t'.join([ORF,name])+'\n')
							o_file_NT.write(''.join(seq)+'\n')
					count = count +1

		print ("Fasta files created for ", count," ORFs")

	else:
		print ("Fasta files already created. Skipping.")
	

	return os.path.join(externalData,"MSA")

def msa_run_clustalw(fasta_path,path_to_clustalo):
	""" 
	Run CLUSTALW to generate an MSA from a fasta file
	
	Args: 
		fasta_path (str): Path to a directory containing the fasta files
		path_to_clustalo (str): Path to the clustal W executable
	"""
	test_ORF = os.listdir(fasta_path)[0]
	if not os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".nt.aln")) :
		for ORF in tqdm(os.listdir(fasta_path)):
			nt_fasta = os.path.join(fasta_path,ORF,ORF+".nt.fa")

			if not os.path.exists(nt_fasta[:-3]+'.aln'):
				os.system(path_to_clustalo+" -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')
	else:
		print ("CLUSTALW MSA already performed. Skipping.")
		
	return 

def build_combined_file(externalData,alignment_path,fasta_path):
	"""
	Write a combined alignment file for all ORFs and all species
	
	"""
	aligned_seq_dict_pickle = os.path.join(externalData,"combined_aligned_seq_dict.pickle")

	if not os.path.exists(alignment_path):
		with open(alignment_path,'w') as f:
			print ("Writing combined alignment data file")
			combined_dict = {}
			# For each ORF
			for ORF in tqdm(os.listdir(fasta_path)):
				if ORF == ".DS_Store":
					continue
				nt_aln = os.path.join(fasta_path,ORF,ORF+".nt.aln")
				combined_dict.update({ORF:{}}) # Add ORF as key in combined dict
				parsed_aln = SeqIO.parse(nt_aln, "fasta")
				
				temp = {} # Temporary dict to store the parsed info for this ORF
				for rec in parsed_aln:
					ORF = rec.id
					description = rec.description.split(' ')[1]
					strain = description.split('.')[0]
					sequence = "".join(list(rec.seq))
					temp.update({(strain,'nt'):sequence})

				# Find indices with gaps in the scer sequence
				gap_indexes = set([i for i, ltr in enumerate(temp['Scer','nt']) if ltr == "-"])
				# Remove position with gaps in Scer from all nt sequences (use map instead?)
				temp['Scer','nt'] = "".join([char for idx, char in enumerate(temp['Scer','nt']) if idx not in gap_indexes])
				temp['Spar','nt'] = "".join([char for idx, char in enumerate(temp['Spar','nt']) if idx not in gap_indexes])
				temp['Smik','nt'] = "".join([char for idx, char in enumerate(temp['Smik','nt']) if idx not in gap_indexes])
				temp['Sbay','nt'] = "".join([char for idx, char in enumerate(temp['Sbay','nt']) if idx not in gap_indexes])
				temp['Ncas','nt'] = "".join([char for idx, char in enumerate(temp['Ncas','nt']) if idx not in gap_indexes])
				temp['Cgla','nt'] = "".join([char for idx, char in enumerate(temp['Cgla','nt']) if idx not in gap_indexes])
				temp['Agos','nt'] = "".join([char for idx, char in enumerate(temp['Agos','nt']) if idx not in gap_indexes])
				temp['Klac','nt'] = "".join([char for idx, char in enumerate(temp['Klac','nt']) if idx not in gap_indexes])
				temp['Calb','nt'] = "".join([char for idx, char in enumerate(temp['Calb','nt']) if idx not in gap_indexes])

				# Translate scer_nt to get the aa alignment
				scer_aa = Seq(temp['Scer','nt'] ,generic_dna).translate()
				temp.update({('Scer','aa'):scer_aa})

				if not (len(temp['Scer','aa'])*3 == len(temp['Scer','nt'])==len(temp['Spar','nt'])==len(temp['Smik','nt'])==len(temp['Sbay','nt'])==len(temp['Ncas','nt'])==len(temp['Cgla','nt'])==len(temp['Agos','nt'])==len(temp['Klac','nt'])==len(temp['Calb','nt'])):
					raise UserWarning('Non-matching lengths when building codon alignment file')

				# Write to the combined_dict and to the summary file
				for key in [('Scer','aa'),('Scer','nt'),('Spar','nt'),('Smik','nt'),('Sbay','nt'),('Ncas','nt'),('Cgla','nt'),('Agos','nt'),('Klac','nt'),('Calb','nt')]: # Forcing the order we want
					if key[1] == 'aa':
						sequence = list(temp[key])
					else:
						initial = ["".join(temp[key][i:i + 3]) for i in range(0, len(temp[key]), 3)]
						sequence = [codon if not "-" in codon else "---" for codon in initial] # if gaps in codon change the whole codon to a gap

					combined_dict[ORF].update({(key[0],key[1]): sequence})
					f.write(ORF+"\t"+key[0]+"."+key[1]+"\t")
					f.write("\t".join(sequence))
					f.write("\n")
					
			# For testing: if need dictionnary 
			print ("Pickling combined alignment data dictionnary for speed")
			pickle_out = open(os.path.join(aligned_seq_dict_pickle),"wb")
			pickle.dump(combined_dict, pickle_out)
			pickle_out.close()
	
	else:
		print ("Combined alignment data file already exists. Skipping.")
		pickle_in = open(aligned_seq_dict_pickle,"rb")
		combined_dict = pickle.load(pickle_in)

	return 



def main():
	eValue = '1e-5'
	minimumCoverage = 0.5
	rawDir = '../../data/external'
	outDir = '../../data/processed'
	blastResultsPath = os.path.join(outDir, 'blast_scer_vs_pdb.txt')
	homologsPath = os.path.join(outDir, 'scer_pdb_homologs.txt')
	curatedPdbDir = os.path.join(outDir, 'assemblies')
	if not os.path.exists(outDir):
		os.makedirs(outDir)

	download_initial_data(rawDir)

	ORFs_with_PPIs = format_biogrid(os.path.join(rawDir, 'bioGrid_all.txt'),
				os.path.join(outDir, 'scer_protein_interactions.txt'),
				os.path.join(outDir, 'scer_protein_interactions_multiple_reports.txt'))

	format_pdb_seqs(os.path.join(rawDir, 'pdb_seqres.txt'),
				os.path.join(outDir, 'pdb_reduced.faa'))

	get_aa_sequences(ORFs_with_PPIs,os.path.join(rawDir, 'Scer.aa'),os.path.join(outDir, 'Scer_reduced.aa'))

	run_blast(os.path.join(outDir, 'pdb_reduced.faa'),
				os.path.join(outDir, 'Scer_reduced.aa'),
				blastResultsPath,
				eValue)

	# Format wanted for the output: a the target name, b the query name, and e the evalue in the blast results file
	oFormat = lambda a, b, e: ('\t'.join([b]+a.split('_'))+'\n')	
	coverage_cutoff_on_BLAST(blastResultsPath,
				homologsPath,
				minimumCoverage, oFormat)
	
	get_pdb_structures(os.path.join(rawDir, 'assemblies'),
				curatedPdbDir,
				homologsPath,
				os.path.join(rawDir, 'pdb_entry_type.txt'))
	
	update_homology_map(curatedPdbDir,
					 homologsPath,
					 os.path.join(outDir, 'curated_chains.txt'),
					 os.path.join(outDir, 'curated_scer_pdb_homologs.txt'))
	
	path_to_clustalo = download_clustal()

	orthologs = read_yeast_orthologs(rawDir)

	# # # Testing orthologues, get the number of Scer proteins with 1,2,3,4,5,6 othologs
	# # print(Counter([len(orthologs[key]) for key in orthologs]))
	
	build_yeast_aligment(rawDir,orthologs,os.path.join(rawDir,'codon_alignment.txt'),path_to_clustalo)
	

if __name__ == '__main__':
	main()
