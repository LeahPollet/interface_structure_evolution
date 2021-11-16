"""
	Step 1 in the pipeline
	
	Sets up the initial data folder structure needed for the analysis (external and processed data folders) 
	Automated download of the external data from Biogrid, PDB and Ensembl. Need to: Update the ftp dowload links when new data is made available 
	Requirements: 
				tqdm==4.42.0
				bio==0.1.0
				brew install p7zip
Author Leah Pollet
leah.pollet@mail.mcgill.ca
"""

import os
import zipfile
import subprocess
from tqdm import tqdm
import csv
import re
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.PDB import PDBParser, PDBIO, Select
from collections import defaultdict, Counter


def download_initial_data(dir):
	"""
	Get initial input files from SGD, PDB, BioGrid, Ensembl.
	Args:
		dir (str): directory to download them into (/data/external directory)
	"""
	print ("----------- Download initial data -----------")
	remoteFiles = [
	# ('http://thebiogrid.org/downloads/archives/Release%20Archive/BIOGRID-3.5.180/BIOGRID-ALL-3.5.180.tab2.zip',"bioGrid_all.txt"),				# Using December 2019 data 
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
	
	# Additional close species (Budding yeasts/Saccharomycetes)
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/zygosaccharomyces_rouxii_gca_000026365/cds/Zygosaccharomyces_rouxii_gca_000026365.ASM2636v1.cds.all.fa.gz',"Zrou.nt"), #Z. Roux
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/zygosaccharomyces_rouxii_gca_000026365/pep/Zygosaccharomyces_rouxii_gca_000026365.ASM2636v1.pep.all.fa.gz',"Zrou.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_dairenensis_cbs_421_gca_000227115/cds/Naumovozyma_dairenensis_cbs_421_gca_000227115.ASM22711v2.cds.all.fa.gz',"Ndai.nt"), #N. Dair
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/naumovozyma_dairenensis_cbs_421_gca_000227115/pep/Naumovozyma_dairenensis_cbs_421_gca_000227115.ASM22711v2.pep.all.fa.gz',"Ndai.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/cyberlindnera_fabianii_gca_001983305/cds/Cyberlindnera_fabianii_gca_001983305.ASM198330v1.cds.all.fa.gz',"Cfab.nt"), # C. Fabi
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/cyberlindnera_fabianii_gca_001983305/pep/Cyberlindnera_fabianii_gca_001983305.ASM198330v1.pep.all.fa.gz',"Cfab.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/tetrapisispora_phaffii_cbs_4417_gca_000236905/cds/Tetrapisispora_phaffii_cbs_4417_gca_000236905.ASM23690v1.cds.all.fa.gz',"Tpha.nt"), #T. Phaf
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/tetrapisispora_phaffii_cbs_4417_gca_000236905/pep/Tetrapisispora_phaffii_cbs_4417_gca_000236905.ASM23690v1.pep.all.fa.gz',"Tpha.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/candida_tropicalis_mya_3404_gca_000006335/cds/Candida_tropicalis_mya_3404_gca_000006335.ASM633v3.cds.all.fa.gz',"Ctro.nt"), #C. Trop
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/candida_tropicalis_mya_3404_gca_000006335/pep/Candida_tropicalis_mya_3404_gca_000006335.ASM633v3.pep.all.fa.gz',"Ctro.aa"),
	
	# Additional distant species (Pezizomycotina/Taphrimomycotina)
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fusarium_graminearum/cds/Fusarium_graminearum.RR1.cds.all.fa.gz',"Fgra.nt"), #F. Gram
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fusarium_graminearum/pep/Fusarium_graminearum.RR1.pep.all.fa.gz',"Fgra.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/neurospora_crassa/cds/Neurospora_crassa.NC12.cds.all.fa.gz',"Ncra.nt"), #N. Cras
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/neurospora_crassa/pep/Neurospora_crassa.NC12.pep.all.fa.gz',"Ncra.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/aspergillus_nidulans/cds/Aspergillus_nidulans.ASM1142v1.cds.all.fa.gz',"Anid.nt"), #A. Nidu
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/aspergillus_nidulans/pep/Aspergillus_nidulans.ASM1142v1.pep.all.fa.gz',"Anid.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/schizosaccharomyces_pombe/cds/Schizosaccharomyces_pombe.ASM294v2.cds.all.fa.gz',"Spom.nt"), #S. Pomb
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/schizosaccharomyces_pombe/pep/Schizosaccharomyces_pombe.ASM294v2.pep.all.fa.gz',"Spom.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/neolecta_irregularis_dah_3_gca_001929475/cds/Neolecta_irregularis_dah_3_gca_001929475.NeoIirr1.0.cds.all.fa.gz',"Nirr.nt"), #N. Irre
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/neolecta_irregularis_dah_3_gca_001929475/pep/Neolecta_irregularis_dah_3_gca_001929475.NeoIirr1.0.pep.all.fa.gz',"Nirr.aa"), 
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota4_collection/protomyces_lactucaedebilis_gca_002105105/cds/Protomyces_lactucaedebilis_gca_002105105.Prola1.cds.all.fa.gz',"Plac.nt"), #P. Lact
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota4_collection/protomyces_lactucaedebilis_gca_002105105/pep/Protomyces_lactucaedebilis_gca_002105105.Prola1.pep.all.fa.gz',"Plac.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/pneumocystis_jirovecii_ru7_gca_001477535/cds/Pneumocystis_jirovecii_ru7_gca_001477535.Pneu_jiro_RU7_V2.cds.all.fa.gz',"Pjir.nt"), #P. Jiro
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota3_collection/pneumocystis_jirovecii_ru7_gca_001477535/pep/Pneumocystis_jirovecii_ru7_gca_001477535.Pneu_jiro_RU7_V2.pep.all.fa.gz',"Pjir.aa"),
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/pneumocystis_murina_b123_gca_000349005/cds/Pneumocystis_murina_b123_gca_000349005.Pneumo_murina_B123_V4.cds.all.fa.gz',"Pmur.nt"), #P. Muri
	('ftp://ftp.ensemblgenomes.org/pub/release-47/fungi/fasta/fungi_ascomycota1_collection/pneumocystis_murina_b123_gca_000349005/pep/Pneumocystis_murina_b123_gca_000349005.Pneumo_murina_B123_V4.pep.all.fa.gz',"Pmur.aa")
	]
	# IntAct data, manual search for now:
	# Advanced search with 'Advanced field option' = Organism, "tax ID"= 4932 (S.cere)
	# to data/external: species_49.txt (S.cere)

	additionalClose_dir = os.path.join(dir,"AdditionalCloseSpecies")
	additionalDist_dir = os.path.join(dir,"AdditionalDistantSpecies")
	if not os.path.exists(additionalClose_dir):
		os.makedirs(additionalClose_dir)
	if not os.path.exists(additionalDist_dir):
		os.makedirs(additionalDist_dir)
	
	for url,shortName in remoteFiles:
		if shortName[0:4] in ["Zrou","Ndai","Cfab","Tpha","Ptan","Ctro"]:
			dir = additionalClose_dir
		if shortName[0:4] in ["Sscl","Fgra","Ncra","Anid","Spom","Nirr","Scom","Plac","Pjir","Pmur"]:
			dir = additionalDist_dir
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
	
def format_PPI_data(ipath1,ipath2, orfpattern,taxID,intAct_orfpattern,summary_files_path,text,multiple_reports = False):
	"""
	Extract physical interactions between proteins in species of interest from BioGrid and IntAct files
	Args: 
		ipath1 (str): path to the Biogrid file as downloaded by download_initial_data
		ipath2 (str): path to the IntAct file as downloaded for the species (downloaded manually from the IntAct website)
		orfpattern (str): regex pattern(s) for the ORF name in the species
		taxID (str): taxonomic identifier(s) for the species 
		text(srt): Name of the species formatted for printing
		multiple_reports (bool): whether or not we want to print reports only for PPIs reported more than once (different pumed IDs)
	"""
	ints = {} 	#Format: {('ORF1','ORF2'):set(pubmedIDs)} 			with ('ORF1','ORF2') sorted alphabetically
	PPIs = {}	#Format: {('ORF1','ORF2'):int(number of reports)} 	with ('ORF1','ORF2') sorted alphabetically
	
	# if not os.path.exists(ipath1):
	# 	raise UserWarning('Input file does not exist: '+ipath1)
	# if not os.path.exists(ipath2):
	# 	raise UserWarning('Input file does not exist: '+ipath2)

	if not os.path.exists(summary_files_path):

		# 1) Get all PPIs in the BioGrid file which: (1) are reported as physical interactions (2) match the specie's taxID (3) match the orf pattern in the species (removes RNA, mitochondrial DNA...) 
		print ("----------- Format Biogrid -----------")
		# Get line number for progress bar
		output = subprocess.check_output(['wc', '-l', ipath1])
		line_num = int(output.split()[0])

		# Parse through the Biogrid file
		with open( ipath1 ) as f:
			for items in tqdm(csv.reader( f, dialect="excel-tab" ),total=line_num):
				if len( items ) == 24:
					orf1, orf2 = items[5], items[6]
					method, pubmedid = items[12], items[14]
					tax1, tax2 = items[15], items[16]
					if method == "physical":
						if (re.search(taxID,tax1) and re.search(taxID,tax2)):
							if re.search( orfpattern, orf1 ) and re.search( orfpattern, orf2 ):
								pair = ( orf1, orf2 ) if orf1 < orf2 else ( orf2, orf1 )
								ints.setdefault(pair, set()).add(pubmedid)
		# print(len(ints.keys()),"PPIs from Biogrid", len([key for key in ints.keys() if len(ints[key])>1]),"are reported more than once")

		# 2) Get all PPIs in the IntAct file which: (1) are reported as physical interactions (+remove colocalization-MI:0403, proximity-MI:2364, rna cleavage-MI:0902, cleavage reaction-MI:0194), (2) match the specie's taxID (3) match the orf pattern in the species (removes RNA, mitochondrial DNA...) 
		print ("----------- Format IntAct -----------")
		# Get line number for progress bar
		output = subprocess.check_output(['wc', '-l', ipath2])
		line_num = int(output.split()[0])
		
		# Parse through the IntAct file
		with open( ipath2 ) as f:
			for items in tqdm(csv.reader( f, dialect="excel-tab" ),total=line_num):
				interaction_type = "MI:0403|MI:2364|MI:0902|MI:0194" #Interaction types to remove: colocalization,proximity, cleavage
				if(not re.search(interaction_type,items[11])):
					if (re.search(taxID,items[9]) and re.search(taxID,items[10])): # In our species of interest
						orf1_check, orf2_check = re.search( intAct_orfpattern, items[4]),re.search( intAct_orfpattern, items[5])
						if(orf1_check and orf2_check):# Regex match to the ORF pattern in the species
							orf1,orf2 = orf1_check.group(0).split("(")[0],orf2_check.group(0).split("(")[0]
							pubmedid = (re.search("pubmed:[0-9]{1,}", items[8])).group(0)[7:]
							pair = ( orf1, orf2 ) if orf1 < orf2 else ( orf2, orf1 )
							ints.setdefault(pair, set()).add(pubmedid)
				
		# 3) Write summary files
		PPIs = {key:len(value) for (key,value) in ints.items()}
		write_PPI_data_summary_files(PPIs, summary_files_path)

	else:
		with open( summary_files_path ) as f:
			for items in csv.reader(f, dialect="excel-tab" ):
				orf1,orf2 = items[0],items[1]
				num_of_reports = int(items[2])
				PPIs.update({(orf1,orf2):num_of_reports})

	uniqueORFs = set([orf for ppi in PPIs.keys() for orf in ppi])
	uniqueORFs_multiple = set([orf for ppi in PPIs.keys() for orf in ppi if (PPIs[ppi])>1])
	
	# Do we want all PPIs or the ones reported more than once?
	if multiple_reports:
		print(len(PPIs.keys()),text,"PPIs between",len(uniqueORFs) ,"unique",text," proteins downloaded from Biogrid+Intact")
		print(len([key for key in PPIs.keys() if (PPIs[key])>1]),"PPIs between",len(uniqueORFs_multiple),"unique",text," proteins are reported more than once")
	else:
		print(len(PPIs.keys()),text,"PPIs between",len(uniqueORFs) ,"unique",text," proteins downloaded from Biogrid+Intact")
	
	return uniqueORFs_multiple

def write_PPI_data_summary_files(PPIs, opath):
	"""
	Write two summary files for a species: (1) the list of all PPIs gathered (2) the list of all PPIs reported more than once
	Args:
		PPIs (dict): dictionary of PPIs in the species. Keys= sorted tuples of interacting ORF names, Values= number of unique pubmed ids reporting the interaction
		opath (str): path to the desired location for the filtered output file (all PPIs)
	"""
	opath1 = opath 								# Path to the summary file (all PPIs)
	opath2 = opath[:-4]+'_multiple_reports.txt' # Path to the summary file (PPIs reported more than once)
	with open(opath1, 'w') as fh1:
		with open(opath2, 'w') as fh2:
			for (orf1,orf2) in PPIs.keys():
				num_of_reports = PPIs[(orf1, orf2)]
				fh1.write('\t'.join([orf1, orf2, str(num_of_reports)])+'\n')
				if(num_of_reports) > 1:
					fh2.write('\t'.join([orf1, orf2, str(num_of_reports)])+'\n')
	return		

def get_scer_pdb_structs(PDBs_Scer):
	"""
	Get PDB IDs for structures in s. cerevisiae into a list. Will compare it with PDBs downloaded later
	Args:
		PDBs_Scer (str): Manually curated list of PDB IDs for structures with Source Organism Type = "natural" and Scientific Name of the Source Organism = "Saccharomyces cerevisiae S288C" OR "Saccharomyces cerevisiae"
	"""
	scer_PDBs=[]
	with open( PDBs_Scer ) as f:
		for items in csv.reader(f):
			for id in items:
				scer_PDBs.append(id)
	scer_PDBs=list(set(scer_PDBs))	
	return scer_PDBs

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
	# Get line number for progress bar
	output = subprocess.check_output(['wc', '-l', ipath])
	line_num = int(output.split()[0])/2
	
	fh = open(opath, "w")
	for record in tqdm(SeqIO.parse(ipath, "fasta"),total=line_num) :
		if "mol:protein" in record.description:
			record.description = record.name
			SeqIO.write(record, fh, "fasta")
	fh.close()
	return

def get_aa_sequences(ORFList, fullSeqFile,opath,chars_to_strip,text):
	"""
	Subset the full s. cer. protein aa sequences to only the ORFs of interest to the analyis 
	Args:
		ORFList (set): set of ORF names of interest.
		fullSeqFile (str): path to file containing all aa sequences in the species. 
		opath (str): path to output file containing result fasta file.
		chars_to_strip (int): number of characters to strip from orfname in the full sequence file to match the ORFList format
		text (str): species name formatted for printing
	"""
	print ("----------- Format aa seqs -----------")
	recordsToKeep = []
	if os.path.exists(opath):
		recordsToKeep = [record.id for record in SeqIO.parse(opath, "fasta")]
		print("Got aa sequences for",len(recordsToKeep),"/",len(ORFList),"orfs in",text,"(removed",len(ORFList)-len(recordsToKeep),"SPNCRNAs, pseudogenes, transposable elements, deleted ORFs...)")
	else:
		fh = open(opath, "w")
		for record in SeqIO.parse(fullSeqFile, "fasta"):
			orf = record.id[:-chars_to_strip]
			sequence = record.seq
			if(orf in ORFList):
				newRecord = SeqRecord(sequence, orf,'','')
				recordsToKeep.append(newRecord.id)
				SeqIO.write(newRecord, fh, "fasta")	
		print("Got aa sequences for",len(recordsToKeep),"/",len(ORFList),"orfs in",text,"(removed",len(ORFList)-len(recordsToKeep),"SPNCRNAs, pseudogenes, transposable elements, deleted ORFs...)")
		# print([orf for orf in ORFList if orf not in recordsToKeep])
	return set(recordsToKeep)

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


def coverage_cutoff_on_blast(ipath, opath, coverage,oFormat):
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
		PDB_ids = list(set([row[1] for row in csv.reader( open( opath ), dialect="excel-tab" )]))
		return PDB_ids
	if not os.path.exists(ipath):
		raise UserWarning('Input file does not exist: '+ipath)
	with open( ipath ) as fh1:
		with open(opath, "w" ) as fh2:
			for items in csv.reader( fh1, dialect="excel-tab" ):
				target = items[0]
				query = items[1]
				evalue = float( items[12] )
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
	PDB_ids = list(set([row[1] for row in csv.reader( open( opath ), dialect="excel-tab" )]))
	return PDB_ids

def get_pdb_structures(rawDir, outDir, homologsPath, pdbMethodsPath):
	"""
	Download, quality control and filter pdb files for all S. cerevisiae orfs
	Once we have the processed/assemblies folder, can run fully without the external/assemblies folder => can be deleted for space
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

	# FTP directory with biological unit coordinate files for all proteins in PDB
	ftpDir = 'ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/'
	format = lambda pdbid: pdbid.lower() + ".pdb1"

	# determine prot-only structures (x-ray,NMR,EM.. other all included)
	xray = {}
	with open(pdbMethodsPath) as fh: # File with the entry type of each PDB
		for items in csv.reader( fh, dialect='excel-tab' ):
			if 'prot' in items[1] :
			# if 'prot' in items[1] and items[2] == 'diffraction':
				xray[items[0]] = 1
	# determine structures to download
	homologs = {}
	with open(homologsPath) as fh:
		for items in csv.reader(fh, dialect="excel-tab" ):
			homologs[items[1]] = 1

	allPDBIDs = [row[1] for row in csv.reader(open(homologsPath), dialect="excel-tab" )]	

	# determine valid structures
	yeastlike = set( homologs.keys() ).__and__( set( xray.keys() ) )
	qcPath = os.path.join(outDir, '../pdb_quality_control_results.csv') # outDir/pdb_quality_control_results.csv
	
	print("Total number of PDB files to download: ", len(yeastlike))
	
	# Read the QC file and add to already done
	if os.path.exists(qcPath):
		print("Download testing and curation already done or partially done.\nSkipping PDB files already downloaded")
		with open(qcPath) as f:
			alreadyDone =[items[0] for items in csv.reader(f)]
	else:
		print ('\nDownloading, checking and filtering PDB files')
		with open(qcPath, 'w') as qcFile:
			qcFile.write('PDB_ID,Download_status\n')
		alreadyDone = []
	
	# Read the already curated files
	list_of_curated = os.listdir(outDir)
	curated = [filename[0:4] for filename in list_of_curated]
	
	# print(len(curated))
	# print(len(alreadyDone))
	# print("intersection: ",len(set(curated).intersection(set(alreadyDone))))
	# print("In curated, not in QC: ",len(set(curated) - set(alreadyDone)))
	# print("In QC, not in curated: ",len(set(alreadyDone) - set(curated)))
	# print(list(set(alreadyDone) - set(curated))[0:3])
	# exit()
	# exit()

	# If the PDB was previously downloaded but not needed anymore, remove it to save space:
	to_remove = [pdb for pdb in alreadyDone if pdb not in yeastlike]
	for pdbID in to_remove:
		rawPath = os.path.join(rawDir, format(pdbID))
		outPath = os.path.join(outDir, format(pdbID))
		if os.path.exists(rawPath):
			os.remove(rawPath)
		if os.path.exists(outPath):
			os.remove(outPath)

	count_alreadyDone,count_tryDownload = 0,0
	with open(qcPath, 'a') as qcFile:
		try:
			for pdbID in tqdm(yeastlike):
				if (pdbID not in alreadyDone):
					print ("> ", pdbID)
					ftpPath = os.path.join(ftpDir, format(pdbID)+".gz")
					rawPath = os.path.join(rawDir, format(pdbID))       # external/assemblies/pdbID
					unzippedPath = os.path.join(rawDir,format(pdbID)) 
					outPath = os.path.join(outDir, format(pdbID))       # processed/assemblies/pdbID
					# PDB file not in the final curated folder of files
					if not os.path.exists(outPath): 
						if not os.path.exists(rawPath):  		# Download and unzip the PDB file              
							os.system('wget -q -P %s %s' % (rawDir, ftpPath))
							os.system('gunzip '+rawPath)
						if not os.path.exists(unzippedPath): 	# Download failed error   
							print ('Download failed \n')        
							qcFile.write(pdbID+',failed_download\n')    
						else:     								# Download sucessful or curation error                                      
							passed, code = curate_pdb_file(rawPath, outPath)
							if not passed:
								print ('Does not pass curation')
								qcFile.write(pdbID+','+code+'\n')
							else:
								qcFile.write(pdbID+','+"OK"+'\n')
					# PDB file already downloaded (e.g. in previous run of the analysis), should be 0 because of the already done check
					else:
						qcFile.write(pdbID+','+"OK"+'\n')
			
				
		finally:
			with open(qcPath) as f:
				all_PDBs_tested = [items[0] for items in csv.reader( f)]
			with open(qcPath) as f:
				all_PDBs_success = [items[0] for items in csv.reader( f) if items[1] =='OK']
			print("Done:",len(set(all_PDBs_success)), "/",len(set(all_PDBs_tested)),",PDB files passed the download, testing and curation steps.")
			
	return
def curate_pdb_file(ipath, opath):
	"""
	Quality check and filter PDB file.
	Args:
		ipath: location of input PDB file.
		opath: location of output filtered PDB file.
	Returns:
		bool: True if structure passed quality checks.
		str:  reason for failure of quality check.
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

def print_summary(summaryfile):
	
	with open(summaryfile) as fh: # QC file from get_pdb_structures
		allPDBs = [row[0] for row in csv.reader(open(summaryfile)) if row[0] != "PDB_ID"]
		# Curated as ok
		curatedPDBs = [row[0] for row in csv.reader(open(summaryfile)) if row[1] == "OK"]
		# Failed download
		noBiolAssemblyPDBs =[row[0] for row in csv.reader(open(summaryfile)) if row[1] == "failed_download"]
		# Heavy atom counts are unreasonable
		countsPDBs =[row[0] for row in csv.reader(open(summaryfile)) if row[1] == "counts"]
		# InterCA distances are unreasonable
		distancePDBs =[row[0] for row in csv.reader(open(summaryfile)) if row[1] == "distance"]
		# No valid residues
		no_valid_residuesPDBs =[row[0] for row in csv.reader(open(summaryfile)) if row[1] == "no_valid_residues"]
		
		print("Total number of PDB biological assembly files we tried to curate: ",len(allPDBs),len(set(allPDBs)))
		print(len(curatedPDBs),"/",len(allPDBs), "files passed curation.")
		print(len(noBiolAssemblyPDBs),"/",len(allPDBs), "do not have a biological assembly file.")
		print(len(countsPDBs)+len(distancePDBs)+len(no_valid_residuesPDBs),"/",len(allPDBs), "failed curation (heavy atoms,CA distances...)")
		print("total: ", len(curatedPDBs)+len(noBiolAssemblyPDBs)+len(countsPDBs)+len(distancePDBs), "=", len(allPDBs))
	return

def update_homology_map(pdbDir, homologsPath, curatedChains, curatedHomologs):
	"""
	Update the list of PDB proteins homology matched to yeast to include only those that passed the curation tests.
	Args:
			pdbDir (str): 			directory with the downloaded PDB structures
			homologsPath (lst): 	file with PDB structures homology mapped to yeast proteins 
			curatedChains (str):	file with the list of non-zero length chains from each PDB structure
			curatedHomologs (lst)	file with PDB structures that passed curation tests homology mapped to yeast proteins
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
				orf, pdbid, chain, evalue = line.split()
				if pdbid in okChains and chain in okChains[pdbid]:
					fhOut.write(line + '\n')
	return 

def homolog_chains(pdbDir, outPath):
	"""
	Remove length zero chains
	Args:
			pdbDir (str):  directory with the downloaded PDB structures
			outPath (str): file with the list of non-zero length chains from each PDB structure
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

def download_clustal():
	# useful if clustal is not compling properly:
	# http://www.clustal.org/clustal2/
	oPath = "../data"
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

def blast_closely_related_species(ortho_species_dir,scer_DB_dir,outDir,summaryFile,scer_closely_related):
	"""
	Maps each S. cervisiae ORF to ortholog ORFs in 8 closely related yeast species
	Args:
		ortho_species_dir (str): path to the directory with amino acid sequences for closely related species to S. cerevisiae (downloaded by process_data.py)
		scer_DB_dir (str): path to file with amino acid sequences for S. cerevisiae orfs in PPIs
		outDir (str): path to directory with results of the alignment 
		summaryFile (str): name of summary file for the alignement (will be in outDir ))
		scer_closely_related (lst): list of species included
	""" 
	
	print ("----------- Find orthologues in closely related species -----------")
	
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	
	for species2 in scer_closely_related:
		run_ortholog_mapping("Scer", species2,scer_DB_dir,os.path.join(ortho_species_dir,species2+'.aa'),outDir)
	
	
	scer_orthologs = read_scer_orthologs([outDir],[scer_closely_related],os.path.join(outDir,summaryFile))

	return scer_orthologs

def read_scer_orthologs(dirs_list,species_list,summaryFile):
	"""
	Read othos in a  dictionnary and print summary file
	dirs_list (lst(str)): list of directories where the alignments to closely related species are located
	species_list (lst(lst)): list of species names for each directory in dirs, in the same order
	summaryFile (str): path to the file where to write a summary of the ortho alignment
	"""

	scer_orthos = defaultdict(dict)
	scer_ortho_species =[]
	for dir,species in zip(dirs_list,species_list):
		for spec in species:
			scer_ortho_species.append(spec)
			blastFile = os.path.join(dir, 'Scer-' + spec + '-blast_stats_coverageCutoff.best.txt')
			with open(blastFile,"r") as f:
				for line in f:
					scer_orthos[line.split()[1]].setdefault(spec,[]).append((line.split()[0],float(line.split()[2])))

	for scer_orf in scer_orthos:
		for orth_species in scer_orthos[scer_orf]:
			(scer_orthos[scer_orf][orth_species]).sort(key=lambda x: float(x[1]))
			scer_orthos[scer_orf][orth_species]= scer_orthos[scer_orf][orth_species][0][0]


	print_ortho_summary_file(summaryFile,scer_ortho_species,scer_orthos,print_sum=False)

	return scer_orthos

def print_ortho_summary_file(summaryFile,scer_ortho_species,scer_orthos,print_sum=False):

	# Print summary file
	if not os.path.exists(summaryFile):
		print("Writing summary file for orthologs of S. cerevisiae ORFs in "+', '.join(scer_ortho_species))
		with open(summaryFile, 'w') as f:
			f.write("Scer\t"+'\t'.join(scer_ortho_species)+"\n")
			for orf in scer_orthos.keys():
				row = [orf]
				for orthoSpecies in scer_ortho_species:
					if orthoSpecies in scer_orthos[orf]:
						row.append(scer_orthos[orf][orthoSpecies])
					else:
						row.append("NA")
				f.write('\t'.join(row)+"\n")
		if (print_sum):
			print("Number of S. cerevisiae ORFs with x /",str(len(scer_ortho_species)),"orthologs:")
			print(Counter([len(scer_orthos[key]) for key in scer_orthos]))
	
	else:
		if (print_sum):
			print("Summary file for orthologs of S. cerevisiae ORFs in "+', '.join(scer_ortho_species),"already exists. Skipping.")
			print("Number of S. cerevisiae ORFs with x /",str(len(scer_ortho_species)),"orthologs:")
			print(Counter([len(scer_orthos[key]) for key in scer_orthos]))
	
	return 
	

def run_ortholog_mapping(species1, species2, species1_path,species2_path,out_folder):
	"""
	Runs BLAST for species 1 against species 2, applies coverage cutoff
	Args:
			species1 (str): short name for species 1 (example: Scer, Spom)
			species2 (str): short name for species 2 (example: Smik, Spar, Sbay, Sjap, Soct, Scry)
			species1_path (str): path to the sequences for species 1
			species2_path (str): path to the sequences for species 2
			out_folder(str): path to the folder where to store results

	""" 
	eValue = '1e-5'
	minimumCoverage = 0.5

	# Run BLAST
	blast_results_path= os.path.join(out_folder,species1+ '-' + species2 + '-blast_stats.best.txt')
	blast_coverage_cutoff_path= os.path.join(out_folder,species1+ '-' + species2 + '-blast_stats_coverageCutoff.best.txt')
	if not (os.path.exists(species1_path) and os.path.exists(species2_path)):
		raise UserWarning('Input aa sequence file does not exist for '+species1+" or "+species2)
	if not os.path.exists(os.path.join(out_folder,species1+'-'+species2+'-blast_stats.best.txt')):
		print ("> Running BLAST for "+species1+" against "+species2)
		run_blast(species1_path, species2_path, blast_results_path, eValue)
	else:
		print ("BLAST of "+species1+" against "+species2+" already done. Skipping.")

	# Apply coverage cutoff
	if not os.path.exists(blast_coverage_cutoff_path):
		oFormat = lambda a, b, e: ('\t'.join([b]+a.split('_')+[str(e)])+'\n')	# Format wanted for the output: a the target name, b the query name, and e the evalue in the blast results file
		coverage_cutoff_on_blast(blast_results_path,blast_coverage_cutoff_path, minimumCoverage,oFormat)
	else: 
		print ("Curation of BLAST results of "+species1+" against "+species2+" already done. Skipping.")

	return

def print_summary_orthos(ortho_mapping,txt):
	print("----------- Print Summary -----------")
	print(txt,":")
	species = set([species for d in ortho_mapping.values() for species in list(d.keys())])
	print(", ".join(species))
	print("Total number of ORFs = ", len(ortho_mapping.keys()))
	count = Counter([len(ortho_mapping[key]) for key in ortho_mapping])
	print("Number of ORFs with mapping in all",len(species),"related species = ", count[len(species)])
	print("Number of no-mapping for each species:")
	for spec in species:
		test = [orf for orf,mapping in ortho_mapping.items() if spec not in mapping.keys()]
		print(spec,":",len(test),"/",len(ortho_mapping.keys()))
	print("Number of time an ORf is dropped just because of this species:")
	for spec in species:
		test = [orf for orf,mapping in ortho_mapping.items() if spec not in mapping.keys() if len(mapping.keys())== len(species)-1]
		print(spec,":",len(test),"/",len(ortho_mapping.keys()))

	
	return

def build_yeast_aligment(sequences_dirs,MSA_dir,orthologs,opath,path_to_clustalo):
	"""
	Build alignment between S cer and related species
	Args:
		sequences_dirs (lst): 			list of path(s) to the folder(s) containing the sequences of related species
		orthologs ({str: {str: str}}): 	s. cer ORF --> species --> ortholog ORF
		opath (str): 					path to save the codon alignment file

	"""
	print ("----------- Build yeast alignments -----------")
	# Read all of the sequence data
	Scer_ORFs = list(orthologs.keys())

	species = set([species for d in orthologs.values() for species in list(d.keys())])
	print("Building MSA for Scer ORFs using related species: \n"+ ", ".join(species))
	print(len({ key:value for (key,value) in orthologs.items() if len(value) == len(species)}),"/",len(orthologs),"ORFs have mapping in all",len(species),"related species.")
	
	sequences = ["Scer.aa","Scer.nt"]+[spec+".nt" for spec in species]+[spec+".aa" for spec in species]

	test_ORF = Scer_ORFs[0]
	fasta_example = os.path.join(MSA_dir,test_ORF,test_ORF+".aa.fa") # Only run read all sequence if needed (long)
	sequences_dict = read_all_sequence_files(sequences_dirs,sequences,fasta_example)

	# Write a FASTA file per yeast ORF with the sequences (nt and aa) from each strain
	fasta_path = write_fasta_MSA_files(MSA_dir,sequences_dict,Scer_ORFs,orthologs)

	# Run a MSA for each ORF (using the fasta file)
	msa_run_clustalw(fasta_path,path_to_clustalo)

	# Combine everything
	build_combined_file(sequences_dirs,opath,fasta_path,["Scer"]+list(species))
	
	return

def read_all_sequence_files(sequences_dirs,sequences_name,fasta_example):
	"""
	Read NT and AA sequence files from all yeast strains into a dictionary
	ORF names formatted to match the format of the orthologues dictionary
	Args: 
		sequences_dirs (lst): 	list of path(s) to the folder(s) containing the sequences of related species
		fasta_example (str):	path to example fasta file, if exists, already done and don't need sequences
	Returns:
		sequences_dict = {strain.xx -> {ORF_name -> [seq]}}
	"""
	sequences_dict ={}

	if not os.path.exists(fasta_example):
		# For each strain, for each ORF update sequences_dict (Sequence dict have all of the ORF sequences for all strains)
		for strain in sequences_name:
			sequences_dict.update({strain:{}})
			for dir in sequences_dirs:
				if os.path.isfile(os.path.join(dir,strain)):
					file = SeqIO.parse(open(os.path.join(dir,strain)),'fasta')
			for fasta in file:
				name, sequence = fasta.id, str(fasta.seq)
				# Modify the ORF name to fit the format in the orthologue dictionary
				if strain[0:4] =='Scer':  # strip the trailing "_mRNA"
					name = name[:-5]

				seq = [sequence[i:i+3] for i in range(0, len(sequence), 3)][:-1]
				
				sequences_dict[strain].update({name:seq})
	
	return sequences_dict

def write_fasta_MSA_files(MSA_dir,sequences_dict,Scer_ORFs,orthologs):
	"""
	For each S.pom ORF with orthologues in all 8 other yeast species ('Sjap', 'Soct', 'Scry','Nirr','Pjir','Pmur','Plac','Scom'),
	create the directory structure and write a fasta files containing the nt sequences and aa sequences of the ORF in all 9 species

	Args:
		MSA_dir: 			path to the folder where we will run the MSA
		sequences_dict: 	Dicionary with aa and nt sequence info for all strains
		Scer_ORFs:			List of S.cer ORFs
		orthologs: 			Dictionary with mapping between yeast ORFs and orthologue ORFs in the other yeast species
	Returns:
		path to a folder containing a FASTA file with aa sequences in the 9 yeast strains, and a 
		FASTA file with nt sequences in the 9 yeast strains, for each ORF with orthologue in all strains
	"""

	count = 0
	test_ORF = Scer_ORFs[0]
	species = list(set([species for d in orthologs.values() for species in list(d.keys())]))

	if not os.path.exists(os.path.join(MSA_dir,test_ORF,test_ORF+".aa.fa")):
		
		# For each S.pom ORF
		for ORF in tqdm(Scer_ORFs): 

			# Check if it has an orthologue in all other species
			if all (k in orthologs[ORF] for k in species):

				# Create the directory to store the 2 fasta files for the ORF
				file_dir = os.path.join(MSA_dir,ORF)
				if not os.path.exists(file_dir):
					os.makedirs(file_dir)

				# 2 fasta files for each ORF with sequences in all species
				o_file_path_NT = os.path.join(file_dir,ORF+".nt.fa") 
				o_file_path_AA = os.path.join(file_dir,ORF+".aa.fa") 

				species_for_fasta = ["Scer"]+species

				# Write the NT fasta file
				if not os.path.exists(o_file_path_NT):
					with open(o_file_path_NT, "w" ) as o_file_NT:
						for spec in species_for_fasta:
							if spec == "Scer":
								ortho_ORF = ORF
							else:
								ortho_ORF = orthologs[ORF][spec]
								# If using the ORFP:/ORFN notation, switch to match the seq dictionnary
								if ortho_ORF[0:4] == "ORFP":
									ortho_ORF = "ORFN"+ortho_ORF[4:]
								# If using the :pep notation, remove it to match the seq dictionnary
								if ortho_ORF[-4:] == ":pep":
									ortho_ORF = ortho_ORF[:-4]
							seq = sequences_dict[spec+".nt"][ortho_ORF]
							
							o_file_NT.write(">"+'\t'.join([ortho_ORF,spec+".nt"])+'\n')
							o_file_NT.write(''.join(seq)+'\n')

				# Write the AA fasta file
				if not os.path.exists(o_file_path_AA):
					with open(o_file_path_AA, "w" ) as o_file_AA:
						for spec in species_for_fasta:
							if spec == "Scer":
								ortho_ORF = ORF
							else:
								ortho_ORF = orthologs[ORF][spec]
							seq = sequences_dict[spec+".aa"][ortho_ORF]

							o_file_AA.write(">"+'\t'.join([ortho_ORF,spec+".aa"])+'\n')
							o_file_AA.write(''.join(seq)+'\n')

				count = count +1


		print ("Fasta files created for ", count," ORFs")

	else:
		print ("Fasta files already created. Skipping.")
	

	return MSA_dir

def msa_run_clustalw(fasta_path,path_to_clustalo):
	""" 
	Run CLUSTALW to generate an MSA from a fasta file
	Uses local version of CLUSTALW for mac.
	Downloaded in download_clustal 
	Need to add skipping DS_Store files
	"""
	print ("----------- Run MSA -----------")
	test_ORF = os.listdir(fasta_path)[0]
	if not os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".nt.aln")) :
		for ORF in tqdm(os.listdir(fasta_path)):
			if ORF == ".DS_Store":
					continue
			nt_fasta = os.path.join(fasta_path,ORF,ORF+".nt.fa")

			if not os.path.exists(nt_fasta[:-3]+'.aln'):
				os.system(path_to_clustalo+" -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')
				# os.system("../../data/clustal/clustalo -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')

	if not os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".aa.aln")) :
		for ORF in tqdm(os.listdir(fasta_path)):
			if ORF == ".DS_Store":
					continue
			aa_fasta = os.path.join(fasta_path,ORF,ORF+".aa.fa")

			if not os.path.exists(aa_fasta[:-3]+'.aln'):
				os.system(path_to_clustalo+" -i "+ aa_fasta + ' -o '+aa_fasta[:-3]+'.aln')
				# os.system("../../data/clustal/clustalo -i "+ nt_fasta + ' -o '+nt_fasta[:-3]+'.aln')

	if os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".aa.aln")) and os.path.exists(os.path.join(fasta_path,test_ORF,test_ORF+".nt.aln")) :
		print ("CLUSTALW MSA already performed. Skipping.")
		
	return 

def build_combined_file(sequences_dirs,opath,fasta_path,species):
	"""
	Args:
		sequences_dirs (lst): 			list of path(s) to the folder(s) containing the sequences of related species
		opath (str): 					path to save the codon alignment file
		fasta_path (str):				path to the folder where with the MSA
		species (lst)					list of all species for combined file
	"""

	aligned_seq_dict_pickle = os.path.join(os.path.split(opath)[0],os.path.basename(opath).split(".")[0]+"_dict.pickle")

	if not os.path.exists(opath):
		with open(opath,'w') as f:
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
					orf = rec.id
					description = rec.description.split(' ')[1]
					strain = description.split('.')[0]
					sequence = "".join(list(rec.seq))
					temp.update({(strain,'nt'):sequence})
				
				# Find indices with gaps in the scer sequence
				gap_indexes = set([i for i, ltr in enumerate(temp['Scer','nt']) if ltr == "-"])

				# Remove position with gaps in Scer from all nt sequences
				for spec in species:
					temp[spec,"nt"] = "".join([char for idx, char in enumerate(temp[spec,"nt"]) if idx not in gap_indexes])
				
				# Translate scer_nt to get the aa alignment
				scer_aa = Seq(temp['Scer','nt'] ,generic_dna).translate()
				temp.update({('Scer','aa'):scer_aa})
				
				if not all(length == len(temp['Scer','aa'])*3 for length in [len(temp[spec,'nt']) for spec in species]):
					raise UserWarning('Non-matching lengths when building codon alignment file')

				# Write to the combined_dict and to the summary file
				species_lst = species[1:] # Remove scer for species list
				for key in [('Scer','aa'),('Scer','nt')]+[(spec,"nt") for spec in species_lst]:
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
		# pickle_in = open(aligned_seq_dict_pickle,"rb")
		# combined_dict = pickle.load(pickle_in)

	return 



def main():
	eValue = '1e-5'
	minimumCoverage = 0.5
	rawDir = '../data/external'
	outDir = '../data/processed'
	PDB_blast_dir = os.path.join(rawDir, 'blastToPDB')
	scer_intAct_file = os.path.join(rawDir, 'species_49.txt')
	biogrid_file = os.path.join(rawDir, 'bioGrid_all.txt')
	scer_PPIs_summary_file = os.path.join(outDir, 'scer_protein_interactions.txt')
	blastResultsPath = os.path.join(PDB_blast_dir, 'blast_scer_vs_pdb.txt')
	homologsPath = os.path.join(outDir, 'scer_pdb_homologs.txt')
	curatedPdbDir = os.path.join(outDir, 'assemblies')
	if not os.path.exists(outDir):
		os.makedirs(outDir)
	if not os.path.exists(PDB_blast_dir):
		os.makedirs(PDB_blast_dir)

	multiple_reports = True

	# 1) Download initial PPI data 
	download_initial_data(rawDir)

	# 2) Extract interactions for s. cerevisiae from BioGrid and IntAct and write summary files
	orfpattern = "Y[A-Z]{2}[0-9]{3}[A-Z\-]{1,}"
	intAct_orfpattern = "Y[A-Z]{2}[0-9]{3}[A-Z\-]{1,}\(locus name\)|R[0-9]{4}[A-Z]{1}\(locus name\)"
	tax = "559292|4932|764097|643680|307796"
	S_cere_orfs = format_PPI_data(biogrid_file,scer_intAct_file,orfpattern,tax,intAct_orfpattern,scer_PPIs_summary_file,"S. cerevisiae",multiple_reports)

	# 3) Try and get PPI data from PDB structures: (TESTING) ***
	# All PDB structure for S.cerevisiae:
	# PDB Statistics: PDB Data Distribution by Source Organism (Natural Source) (https://www.rcsb.org/stats/distribution-source-organism-natural)
	# Saccharomyces cerevisiae S288C (555)
	# Saccharomyces cerevisiae (502)
	# All already included except 38 PDB_IDs
	# 5bq6 (1 prot + antibiotic) (1)
	# 3u2f,3ud0,4j2w,1aqr,1a4e,2h53,1fmy,4j31,3u32,1aqq,1aqs,3u2y,2h50,5bqj,5bqa,1aoo,5bps, (1 protein) (17)
	# 1lu3,1lmv,1tra,3izd,6tna,1fcw,1tn2,1yfg,1i9v,1tn1,486d,6lvr,1lpw,4tra,2tra,4tna,1ehz,3tra, (nucleic acid - tRNA) (18)
	# PPIs between different organisms: 3n2d,1vad (2)
	# + some other PDB ids from this list are removed because low coverage
	
	# Got all of the PDBIDs manually into ../data/external/PDB_scer.txt
	# Will check from this list how many of the PDB structs here are not included in the final PDBs downloaded!
	scer_PDBs = get_scer_pdb_structs(os.path.join(rawDir, 'PDB_scer.txt'))

	# 4) Filter all PDB sequences for only proteins sequences
	format_pdb_seqs(os.path.join(rawDir, 'pdb_seqres.txt'),os.path.join(PDB_blast_dir, 'pdb_reduced.faa'))

	# 5) Subset the full aa sequence file in the S. cerevisiae to only sequences for orfs in PPIs
	S_cere_orfs = get_aa_sequences(S_cere_orfs,os.path.join(rawDir, 'Scer.aa'),os.path.join(PDB_blast_dir, 'Scer_reduced.aa'),5,"S. cerevisiae")
	
	# 6) Blast sequences for orfs in PPIs vs. sequences on PDB 
	run_blast(os.path.join(PDB_blast_dir, 'pdb_reduced.faa'),os.path.join(PDB_blast_dir, 'Scer_reduced.aa'),blastResultsPath,eValue)

	# 7) Format BLAST results (No coverage cutoff)
	# ../data/processed/blast_scer_vs_pdb.txt contains all results from BLAST of scer sequences in PPIs vs PDB proteins
	oFormat = lambda a, b, e: ('\t'.join([b]+a.split('_')+[str(e)])+'\n')	# Format wanted for the output: a the target name, b the query name, and e the evalue in the blast results file
	PDBs = coverage_cutoff_on_blast(blastResultsPath,os.path.join(outDir, 'blast_scer_vs_pdb.txt'),0,oFormat)

	# 8) Coverage cutoff on BLAST, only keep blast results with coverage >= minimumCoverage 
	# (NOT NEEDED, we have a composite coverage cutoff in build_protein_models.py?)
	# ../data/processed/scer_pdb_homologs.txt contains all results from BLAST of scer sequences in PPIs vs PDB proteins with coverage >=50%
	oFormat = lambda a, b, e: ('\t'.join([b]+a.split('_')+[str(e)])+'\n')	# Format wanted for the output: a the target name, b the query name, and e the evalue in the blast results file
	PDBs_coverage = coverage_cutoff_on_blast(blastResultsPath,homologsPath,minimumCoverage, oFormat)

	# Test is some PDB structure from yeast *** should be added or does BLAST already find them all!
	# scer_PDBs = [x.lower() for x in scer_PDBs]
	# test1 = set([x for x in scer_PDBs if x not in PDBs])
	# test2 = set([x for x in scer_PDBs if x not in PDBs_coverage])

	# # 9) Download all of the PDB structures 
	# get_pdb_structures(os.path.join(rawDir, 'assemblies'),
	# 	curatedPdbDir,
	# 	homologsPath,
	# 	os.path.join(rawDir, 'pdb_entry_type.txt'))

	# 10) print summary of PDB quality control results
	print_summary(os.path.join(outDir, 'pdb_quality_control_results.csv'))

	# 11) Final homology mapping of S.cerevisiae proteins involved in PPIs to PDB structures (only biological assemblies that passed curation)
	update_homology_map(curatedPdbDir,homologsPath,
					 os.path.join(outDir, 'curated_chains.txt'),
					 os.path.join(outDir, 'curated_scer_pdb_homologs.txt'))

	# 12) Download clustal exectuable 
	path_to_clustalo = download_clustal()

	# 13) Find ortholgs for each S. cerevisiae protein in PPI in closely related species (Set up needed to calculate dN/dS) 
	# closely related species for S. cerevisiae: S. paradoxus, S. mikatae, S. bayanus, N. castellii, C. glabrata, E. gossypii, K. lactis, C. albicans
	# Format scer_closely_related: {scerORF:{"Smik":smikORF,"Spar":sparORF,"Sbay":sbayORF...}} with scerORF all ORFs in S. cerevisiae PPIs, smikORF,sparORF,sbayORF the orthologous ORFs with the best evalue in the BLAST alignment
	scer_closely_related = ['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb']
	closely_related_ortho_mapping = blast_closely_related_species(rawDir,os.path.join(PDB_blast_dir, 'Scer_reduced.aa'),os.path.join(rawDir,"blastToCloselyRelated"),"blastToCloselyRelated.summary.txt",scer_closely_related)
	
	# Adding more closely related species
	scer_additional_closely_related = ['Ctro', 'Tpha','Cfab','Ndai','Zrou']
	additionalClose_dir = os.path.join(rawDir,"AdditionalCloseSpecies")
	additional_closely_related_mapping = blast_closely_related_species(additionalClose_dir,os.path.join(PDB_blast_dir, 'Scer_reduced.aa'),os.path.join(rawDir,"blastToAdditionalCloselyRelated"),"blastToAdditionalCloselyRelated.summary.txt",scer_additional_closely_related)
	
	# Merge current and additional closely related dir, print summary file
	# TESTING (All keys in additional_closely_related_mapping are also in closely_related_ortho_mapping)
	combined_close =  {}
	for key,value in closely_related_ortho_mapping.items():
		if key in additional_closely_related_mapping:
			val = {k:v for (k,v) in value.items()}
			val.update(additional_closely_related_mapping[key])
			combined_close[key] = val
		else:
			val = {k:v for (k,v) in value.items()}
			combined_close[key] = val

	file = os.path.join(rawDir,"blastToAdditionalCloselyRelated/blastToAdditionalCloselyRelated_original.summary.txt")
	print_ortho_summary_file(file,['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Ctro', 'Tpha','Cfab','Ndai','Zrou'],combined_close,print_sum=False)
			
	# Adding more distantly related species
	scer_additional_distantly_related = ['Fgra','Ncra','Anid','Spom','Nirr','Plac','Pjir','Pmur']
	additionalDistant_dir = os.path.join(rawDir,"AdditionalDistantSpecies")
	additional_distantly_related_mapping = blast_closely_related_species(additionalDistant_dir,os.path.join(PDB_blast_dir, 'Scer_reduced.aa'),os.path.join(rawDir,"blastToAdditionalDistantlyRelated"),"blastToAdditionalDistantlyRelated.summary.txt",scer_additional_distantly_related)
	
	
	# Merge current and additional distantly related, print summary file
	# TESTING (All keys in additional_distantly_related_mapping are also in closely_related_ortho_mapping)
	combined_distant =  {}
	for key,value in closely_related_ortho_mapping.items():
		if key in additional_distantly_related_mapping:
			val = {k:v for (k,v) in value.items()}
			val.update(additional_distantly_related_mapping[key])
			combined_distant[key] = val
		else:
			val = {k:v for (k,v) in value.items()}
			combined_distant[key] = val

	file = os.path.join(rawDir,"blastToAdditionalDistantlyRelated/blastToAdditionalDistantlyRelated_original.summary.txt")
	print_ortho_summary_file(file,['Sbay', 'Smik', 'Spar','Ncas','Cgla','Agos','Klac','Calb','Fgra','Ncra','Anid','Spom','Nirr','Plac','Pjir','Pmur'],combined_distant,print_sum=False)
	
	# 14) Checking which species produces the most NAs in a summary file
	# print_summary_orthos(closely_related_ortho_mapping,"Orginal closely related species")
	# print_summary_orthos(additional_closely_related_mapping,"Additional closely related species")
	# print_summary_orthos(combined_close,"Orginal closely related species + additional closely related species")
	# print_summary_orthos(additional_distantly_related_mapping,"Additional distantly related species")
	# print_summary_orthos(combined_distant,"Orginal closely related species + additional distantly related species")
	
	# Running 3 analysis: original, original + addition close species, distant species
	# NT run multiples for the MSA,
	# Afterwards just add a different name for each to models?
	# Will NT run evolutionary rate calculations multiple times

	# 15) Build combined file with the codon alignment for all proteins/all species
	# Base set of species
	# build_yeast_aligment([rawDir],os.path.join(rawDir,'MSA'),closely_related_ortho_mapping,os.path.join(rawDir,'codon_alignment.txt'),path_to_clustalo)
	
	# Base set + additional closely related
	# build_yeast_aligment([rawDir,os.path.join(rawDir,"AdditionalCloseSpecies")],os.path.join(rawDir,'MSA_AdditionalCloseSpecies'),combined_close,os.path.join(rawDir,'codon_alignment_AdditionalCloseSpecies.txt'),path_to_clustalo)
	
	# Base set + Distantly related species 
	build_yeast_aligment([rawDir, os.path.join(rawDir,"AdditionalDistantSpecies")],os.path.join(rawDir,'MSA_AdditionalDistantSpecies'),combined_distant,os.path.join(rawDir,'codon_alignment_AdditionalDistantSpecies.txt'),path_to_clustalo)
	
	

if __name__ == '__main__':
	main()



	
