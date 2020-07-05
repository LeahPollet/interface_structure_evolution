
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
	python consurf.py
	python bin_codons.py 
	python calc_dnds.py
#Set up for plotting	
	python write_for_R.py

-----------
Requirements:
-----------

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

-------------------------
Description of data files:
-------------------------
External:
--------
Raw data used as input in this analysis. The files in this folder are downloaded from the databases mentioned above by scripts in the pipeline and stored here with little to no processing. Those include: 
- Raw downloads from Biogrid (bioGrid_all.txt), PDB (assemblies, pdb_entry_type.txt, pdb_seqres.txt,resolu.idx) and EBI (pdb_chain_uniprot.csv) 
- Downloaded nucleotide and amino acid sequences for species of interest to this analysis (*.nt, *.aa)
- Results of BLAST alignments between those species (*blast-stats*.txt), and 
- Multiple sequence alignment results between those species for various yeast ORFs (MSA, codon_alignment.txt, combined_aligned_seq_dict.pickle).

Processed:
--------
Processed output files from various scripts in this analysis. The files in this folder summarize results obtained at intermediary steps in the pipeline and could be useful to other analysis and research work. Those include: 
- Curated lists of PPIs reported in S. cerevisiae by at least one (scer_protein_interactions.txt) or more than one (scer_protein_interactions_multiple_reports.txt) experiments on the BioGrid database 
- Results of homology mapping between S.cerevisiae proteins and PDB ids and chains (scer_pdb_homologs.txt,curated_scer_pdb_homologs.txt)
- Quality control reports on PDB structure files (pdb_quality_control_results.csv,curated_chains.txt,pdb_lengths.csv)
- The final structural models of PPIs generated in this analysis in a pickled python dictionnary format (modeled_proteins.pkl)

evolutionary_rates:
------------------
Intermediary steps and results of the dN/dS calculations on residues in our PPI models. 
- The binned_codons folder contains the results of binning residues in our models randomly and according to various structural properties. The files in this folder are used to compute dN/dS ratios for the bins.
- The RSA_distributions folder contains summary of the average values of RSA (monomer and complex RSA) for the residue bins.
- The dnds folder contains the results of running dN/dS calculations on each residue bin. the *.dnds files contain summary statistics such as the number of bootstrap used, number of nucleotides used, mean and standard error for the dN/dS calculations.

clustal/dssp-2.0.4/consurf/paml4.9j:
--------
Local version of software and tools used in this analysis. Downloaded as per website instructions/

Any additional information can be requested by email at leah.pollet@mail.mcgill.ca.
