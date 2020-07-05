#! /usr/bin/env python

"""
Wrapper for doing dN/dS calculations
===============================================
Author: Eric Franzosa (eric.franzosa@gmail.com)
"""

import os, sys, re, glob, argparse
import shutil, random, json
from numpy import mean, std
from Bio.Phylo.PAML import codeml
from Bio import Seq, SeqIO


def read_fasta ( path ):
	fdict = {}
	with open( path ) as fh:
		for record in SeqIO.parse( fh, "fasta" ):
			fdict[record.name] = str( record.seq )
	return fdict

def make_codons ( fdict ):
	headers = list(fdict.keys())
	sequences = list(fdict.values())
	assert len( set( map( len, sequences ) ) ) == 1, \
		"Not all coding sequences were the same length"
	seqlen = len( sequences[0] )
	assert seqlen % 3 == 0
	for h, s in fdict.items():
		fdict[h] = [s[k:k+3] for k in range( 0, seqlen, 3 )]
	bad_columns = {}
	for codons in fdict.values():
		for i, codon in enumerate( codons ):
			if re.search( "[^ACGT]", codon ) or codon in "TAA|TGA|TAG":
			   bad_columns[i] = 1
	for h, codons in fdict.items():
		fdict[h] = [codon for i, codon in enumerate( codons ) if i not in bad_columns]
	return fdict

def resample_codons ( fdict ):
	n = len( list(fdict.values())[0] )
	choices = [random.choice( range( n ) ) for k in range( n )]
	fdict2 = {}
	for h, codons in fdict.items():
		fdict2[h] = [codons[i] for i in choices]
	return fdict2

def drop_fdict_alignment ( fdict, path ):
	fdict2 = {h:"".join( codons ) for h, codons in fdict.items()}
	n = len( fdict2 )
	seqlen = len( list(fdict2.values())[0] )
	with open( path, "w" ) as fh:
		fh.write(str(n)+"\t"+str(seqlen)+"\n")
		for h, seq in fdict2.items():
			fh.write(str(h)+"\n")
			fh.write(str(seq)+"\n")
		   

def run_codeml ( fdict, tree, c_working_dir ,pathToCodeML):
	if not os.path.exists(c_working_dir):
		os.makedirs(c_working_dir)
	c_temp_align = os.path.join( c_working_dir, "temp.nuc" )
	c_temp_out = os.path.join( c_working_dir, "temp.out" )
	c_temp_tree = os.path.join( c_working_dir, "tree" )

	
	# gotcha: paml wants everything in the working dir
	shutil.copy( tree, c_temp_tree )
	# make the temp alignment file
	drop_fdict_alignment( fdict, c_temp_align )
	
	# command
	cml = codeml.Codeml(
		alignment = c_temp_align,
		tree = c_temp_tree,
		out_file = c_temp_out,
		working_dir = c_working_dir,
		)
	# inherited from earlier projects
	cml.set_options( verbose=True )
	cml.set_options( noisy=0 )
	cml.set_options( runmode=0 )
	cml.set_options( seqtype=1 )
	cml.set_options( model=0 )
	cml.set_options( getSE=0 )
	cml.set_options( cleandata=0 )
	# execute
	results = cml.run(command = pathToCodeML)
	
	# results = cml.run()
	
	try:
		dn = results["NSsites"][0]["parameters"]["dN"]
		ds = results["NSsites"][0]["parameters"]["dS"]
		dnds = dn / float( ds )
	except:
		dn, ds, dnds = None, None, None
	return dn, ds, dnds

def main ( ):

	# argument parsing (python argparse)
	parser = argparse.ArgumentParser()
	parser.add_argument( "-a", '--alignment', help='aligned nucleotides in fasta format' )
	# parser.add_argument( "-c", '--codefile', help='absolute path to the codeml executable' )
	parser.add_argument( "-t", '--tree', help='tree file for seqs from alignment' )
	parser.add_argument( "-o", '--output', default=None, help='output file (default stdout)' )
	parser.add_argument( "-b", '--bootstraps', default=0, type=int, help='# of bootstraps (default 0)' )
	parser.add_argument( "-w", '--workingdir', default=None, help='working directory' )
	parser.add_argument( "-p", '--pathToCodeML', default=None, help='path to codeml executable' )
	args = parser.parse_args()


	# load the fastas
	fdict = make_codons( read_fasta( args.alignment ) )
	seqlen = len(list(fdict.values())[0] ) * 3 # this is a list of codons
	
	# get defaults
	dn, ds, dnds = run_codeml( fdict, args.tree, args.workingdir ,args.pathToCodeML)

	# do bootstraps
	boot_dn, boot_ds, boot_dnds = [], [], []
	for trial in range( args.bootstraps ):
		fdict2 = resample_codons( fdict )
		dn2, ds2, dnds2 = run_codeml( fdict2, args.tree, args.workingdir ,args.pathToCodeML)
		boot_dn.append( dn2 )
		boot_ds.append( ds2 )
		boot_dnds.append( dnds2 )
	
	# format the output
	def bootform ( boot_list ):
		if len( boot_list ) > 0:
			boot_list = [k for k in boot_list if k is not None]
			return [str( len( boot_list ) ), mean( boot_list ), std( boot_list )]
		else:
			return []
	def makestring ( items ):
		return "\t".join( [k if type( k ) is str else "%.4f" % ( k ) for k in items] )

	headers = ["param", "val"] + ( ["trials", "mean", "std"] if args.bootstraps > 0 else [] )
	with open( args.output, "w" ) if args.output is not None else sys.stdout as fh:
		fh.write("# of NTs: "+str(seqlen)+"\n")
		fh.write(makestring( headers )+"\n") 
		fh.write(makestring( ["dN", dn] + bootform( boot_dn ) )+"\n") 
		fh.write(makestring( ["dS", ds] + bootform( boot_ds ) )+"\n")
		fh.write(makestring( ["dN/dS", dnds] + bootform( boot_dnds ) )+"\n") 
		

if __name__ == "__main__":
	main( )
