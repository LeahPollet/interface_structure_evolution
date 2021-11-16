"""
Models of yeast PPIs that contain structural and
evolutionary residue-level annotations.

"""
import numpy as np
from Bio import Seq
from amino_acid import AminoAcid

class Residue:
	def __init__(self, id, indexInSeq, aaLetter,cdns,num_closely_related,rsaMonomer,rsaComplex,drsa):
		self.id = id                    # Identifier notation for the residue in the PDB file (to match with dssp output: format examples MET1:A (AaPosition:chain))
		self.indexInSeq = indexInSeq    # Index of the residue with respect to the sequence of the entire protein as saved in the model (0:prot length)
		self.aa = AminoAcid(aaLetter)
		if len(cdns) != num_closely_related+1:
			raise ValueError('Wrong number of codons. Expected '+str(num_closely_related+1))
		if not all([len(c) == 3 for c in cdns.values()]):
			raise ValueError('codons need to be 3 letters long')
		self.codons = cdns                       # Dict of yeast strain to 3-letter string
		self.rsaMonomer = rsaMonomer             # Relative solvent accessibility in uncomplexed state
		self.rsaComplex = rsaComplex             # Relative solvent accessibility in uncomplexed state
		if (rsaMonomer-rsaComplex) != drsa:
			raise ValueError('Error in solvent accessibility calculations')
		self.interfacial = False
		if (drsa != 0):
			self.interfacial = True
		self.dRsa = drsa                   # Delta RSA
		self.interRRC = np.nan             # Residue contacts with interaction partners.
		self.dCenter = np.nan              # Distance to the geometric center of the interface
		self.dEdges = np.nan               # Distance to the edges of the interface (most exposed surface upon complex formation)
		
		self.substitution = np.nan                      # 1 if there is a substitution in closely related species, 0 otherwise
		self.consurfDB_score = np.nan                   # Pre-computed Consurf score downloaded from ConSurf DB
		self.consurfDB_conf_int = np.nan                  
		self.consurf_score = np.nan                     # Manually computed ConSurf score (Unsing rate4site with MSA and tree as input)
		self.conf_int_consurf_score = (np.nan,np.nan)
		self.consurf_std = np.nan
		self.dN = np.nan
		self.dS = np.nan
		self.dNdS = np.nan                              # dN/dS score

	def __str__(self):
		if self.interfacial:
			txt= ""
		else:
			txt= "non-"
		return("Modeled "+txt+"interfacial residue \n" +
				"ID: "+self.id+" , IndexInSeq: "+str(self.indexInSeq)+" , AminoAcid: "+self.aa.letter+" , Codons:  "+str(self.codons)+"\n"+
				"MonomerRSA: "+str(self.rsaMonomer)+" , ComplexRSA: "+str(self.rsaComplex)+" , DeltaRSA: "+str(self.dRsa) +"\n"+
				"interRRC: "+str(self.interRRC)+" , dCenter: "+str(self.dCenter)+" , dEdges: "+str(self.dEdges)+ "\n"+
				"Substitution: "+str(self.substitution)+" , Consurf Score (DB): "+str(self.consurfDB_score)+" , Consurf Score (rate4Site): "+str(self.consurf_score)+" , dN/dS: "+str(self.dNdS))

class Interface:
	def __init__(self,orfName, interactionPartner,pdb,pdbPartner ):
		self.orf = orfName
		self.interactionPartner = interactionPartner
		self.pdb = pdb                      # PDB ID of aligned structure
		self.pdbPartner = pdbPartner        # PDB ID of aligned structure for partner
		self.residues = []

	def add_residue(self, res):
		self.residues.append(res)

	def get_size(self):
		return len(self.residues)
		
	def __str__(self):
		return('Modeled interface for ' + self.orf + "\n"+
				"The interface is modeled for the interaction with "+self.interactionPartner+" using the PDB file "+self.pdb[0]+" (chains "+self.pdb[1]+" and "+self.pdbPartner[1]+")\n"+
				"The interface contains "+ str(len(self.residues))+" modeled interfacial residues")
class Protein:
	def __init__(self,orfName,interactionPartner,species,pdb,evalue,coverage, pdbPartner,aaSeq,seqPDB,pdbWholeProteinSeqs,orthologs,resolution):
		self.orf = orfName
		self.interactionPartner = interactionPartner
		self.species = species
		self.pdb = pdb     # PDB ID of aligned structure
		self.evalue = evalue
		self.coverage = coverage
		self.pdbPartner = pdbPartner     # PDB ID of aligned structure for partner
		self.seq = aaSeq                # original aa sequence of s. cer. protein
		self.pdbStructureSeq = seqPDB       # aa seq of just the structurally imaged residues of aligned PDB chain
		self.pdbWholeProteinSeq = pdbWholeProteinSeqs    # aa seq of entire protein of PDB chain
		self.resolution = resolution
		self.orthologs = orthologs 		# Mapping to orthologs in related species
		self.orthologsSpecies = list(orthologs.keys())	
		self.residues = []
		self.interface = []  

	def add_interface(self, interf):
		self.interface.append(interf)

	def add_residue(self, res):
		# Adding an interfacial residue to a protein
		if res.interfacial:
			if len(self.interface)==0:
				raise ValueError('Add an interface before adding an interfacial residue')
			elif len(self.interface)>1:
				raise ValueError('Multiple interfaces in the protein?')
			else:
				self.interface[0].add_residue(res)
		# For all residues we add
		self.residues.append(res)
		
	
	def __str__(self):
		if len(self.interface)>0:
			txt= "("+str(self.interface[0].get_size())+" interfacial residues)"
		else:
			txt = ""

		if self.interactionPartner == None:
			return("Modeled S. cerevisiae protein: "+str(self.orf)+"\n"+
				"The protein is modeled as a single protein (interacting partner lost in the species) using the PDB file "+str(self.pdb[0])+" (chain "+str(self.pdb[1])+")\n"+
				"The protein has closely related orthologs: "+str(self.orthologs)+"\n"+
				"The protein has "+str(len(self.residues))+" modeled residues\n"+
				"The protein has "+str(len(self.interface))+" modeled interface "+txt)

		else:
			return("Modeled S. cerevisiae protein: "+str(self.orf)+"\n"+
				"The protein is modeled interacting with "+str(self.interactionPartner)+" using the PDB file "+str(self.pdb[0])+" (chains "+str(self.pdb[1])+" and "+str(self.pdbPartner[1])+")\n"+
				"The protein has closely related orthologs: "+str(self.orthologs)+"\n"+
				"The protein has "+str(len(self.residues))+" modeled residues\n"+
				"The protein has "+str(len(self.interface))+" modeled interface "+txt)

class PPI:
	
	def __init__(self,orf1, orf2, reports):
		self.orf1 = orf1
		self.orf2 = orf2
		self.proteins = []
		self.numberOfReports = int(reports)      	# Number of reports for the PPI 
		self.multipleReports = False            	# True of the PPI is reported more than once
		if(int(self.numberOfReports) > 1 ):
			self.multipleReports = True   
		self.noPartner = False
		if(self.orf1 == None or self.orf2 == None):
			self.noPartner = True
			
	def add_protein(self, prot):
		self.proteins.append(prot)
		
		
	def __str__(self):
		
		if (self.noPartner == True):
			return ('Modeled S. cerevisiae single protein between '+ str(self.orf1) +" and "+ str(self.orf2)+"\n"+
				"Reported "+ str(self.numberOfReports)+ " times by independent experiments \n")
				
		else:
			return ('Modeled S. cerevisiae PPI between '+ str(self.orf1) +" and "+ str(self.orf2)+"\n"+
				"It was reported "+ str(self.numberOfReports)+ " times by independent experiments \n")
				
				

