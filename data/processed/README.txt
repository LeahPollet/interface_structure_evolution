-------------------------
Description of processed data files:
-------------------------
The files in this folder are output files from various scripts in this analysis. They summarize results obtained at intermediary steps in the pipeline and could be useful to other analysis and research work. 

PPI data :
----------
scer_protein_interactions.txt (Curated lists of PPIs reported in S. cerevisiae by at least one experiment on the BioGrid database)
scer_protein_interactions_multiple_reports.txt (Curated lists of PPIs reported in S. cerevisiae by more than one experiment on the BioGrid database)

Structure data :
----------
pdb_quality_control_results.csv (Quality control results on PDB files in this analysis)
curated_chains.txt (Curated list of chains for PDB files in this analysis)
pdb_lengths.csv (Structure and sequence lengths for PDB files in this analysis)

Alignment data : 
----------
blast_scer_vs_pdb.txt (Results of BLAST alignment between S.cere and PDB sequences)
scer_pdb_homologs.txt (Reformatted mapping between S.cere ORF names and PDBid/chain ids)
curated_scer_pdb_homologs.txt (Curated mapping with only chains of interest/high quality)
pdb_reduced.faa* (BLAST database files used in the alignment)

Structural models :
----------
modeled_proteins.pkl (Structural models of PPIs in this analysis in a pickled python dictionnary format)

Any additional information can be requested by email at leah.pollet@mail.mcgill.ca. 
