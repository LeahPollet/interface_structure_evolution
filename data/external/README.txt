-------------------------
Description of external data files:
-------------------------
The files in this folder are the raw data used as input in this analysis. They are files automatically downloaded from various databases by scripts in the pipeline and stored here with little to no processing.
 
PPI data :
----------
bioGrid_all.txt (raw download from the BioGrid, list of all PPIs on the database)

Sequence data :
----------
****.aa (amino acid sequence data for a yeast species of interest to this analysis)
****.nt (nucleotide sequence data for a yeast species of interest to this analysis)

Structure data :
----------
pdb_entry_type.txt (raw download from the PDB, entry type for each PDB id on the database)
pdb_seqres.txt (raw download from the PDB, refseq for each PDB id on the database)
resolu.idx (raw download from the PDB, resolution for each PDB id on the database)
pdb_chain_uniprot.csv (raw download from the EBI, chain information for PDB ids)

Alignment data :
----------
*-*-blast_stats.best.txt (results of BLAST alignments between two yeast species)
*-*-blast_stats_coverageCutoff.best.txt (results of BLAST with coverage cutoff applied)
Scer.aa.* (BLAST database files used in the alignment)
MSA/*.nt.fa (fasta file with sequence in 9 yeast species for various ORFs) 
MSA/*.nt.aln (ClustalW MSA file between 9 yeast species for various ORFs) 
codon_alignment.txt (Combined MSA results for all ORFs, reformatted as a codon alignment)
combined_aligned_seq_dict.pickle (Combined MSA results for all ORF, formatted as a pickled python dictionary)

