The scripts in this folder are the only ones needed to run the pipeline from scratch:

---------------------------------
To run the pipeline from scratch:
---------------------------------
Download the 'scripts' folder into a local directory of your choosing (here interface_structure_evolution)

#Initial clean up
	cd interface_structure_evolution/scripts/modeling
	find . -name '.DS_Store' -type f -delete
#Modeling
	python process_data.py
	python pdb_coverage.py
	python build_protein_models.py
#Evolutionary rate calculations
	python consurfDB.py
	python rate4site.py
	python bin_codons.py 
	python calc_dnds.py
#Set up for plotting	
	python write_for_R.py

-----------
Requirements:
-----------
Additional scripts:
----------
- amino_acid.py: Requirement for build_protein_models.py: Class definitions for amino acids
- modeled_proteins.py: Requirement for build_protein_models.py: Class definitions for the PPI models.
- pattern_align_local.py: Requirement for build_protein_models.py: Aligns two sequences and transfers annotations
- sasa_scan.py: Requirement for build_protein_models.py: Wrapper script for doing solvent accessibility calculations using DSSP.
- run_codeml.py: Requirement for calc_dnds.py: Wrapper for doing dN/dS calculations with codeml.

Software and packages:
----------
(Unless otherwise specified all software and packages used in the pipeline were downloaded and compiled as per website instructions, and the path to the final executable needed was added to the PATH environment variable)

Python 3.7.3
Protein-Protein BLAST 2.7.1+
Clustal Omega 1.2.3 - Downloaded and compiled locally in the process_data.py script
DSSP 2.0.4 - Downloaded and compiled locally in the build_protein_models.py script
paml4.9j - Downloaded and compiled locally in the calc_dnds.py script

Libraries:
----------
Standard Python 3.7.3 libraries 
feather 0.4.0
scipy 1.4.1
json 2.0.9
tqdm 4.42.0
Bio 1.76
numpy 1.18.1
pandas 1.0.1
matplotlib 3.1.2

Databases:
----------
http://thebiogrid.org
https://pdb101.rcsb.org
https://www.ebi.ac.uk
http://ensemblgenomes.org
https://yeastgenome.org
https://consurfdb.tau.ac.il

Running the pipeline will download and generate files into a data and a plot folder

Any additional information can be requested by email at leah.pollet@mail.mcgill.ca.
