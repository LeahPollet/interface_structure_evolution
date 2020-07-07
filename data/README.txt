-------------------------
Description of data files:
-------------------------
External:
--------
Raw data used as input in this analysis. The files in this folder are downloaded various databases by scripts in the pipeline and stored here with little to no processing. Those include: 
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

consurf :
------------------
Intermediary steps and results from ConSurf score calculations obtained from the ConSurf database for this analysis.
- The downloaded files are in the consurfDB_files
- Input file used is in the input_file folder
- Log files for both the download from the database and the match with structural models of PPIs are also available

clustal/dssp-2.0.4/paml4.9j:
--------
Local version of software and tools used in this analysis. Downloaded as per website instructions/

Any additional information can be requested by email at leah.pollet@mail.mcgill.ca.