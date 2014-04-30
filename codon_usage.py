#! /usr/bin/python3

#Module that provides various codon usage stats. It is dependent on Bio.Seq
from Bio.Seq import Seq
import Bio
class CodonUsageTable:
	CODON_DICT={'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0, 
	'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0, 
	'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0, 
	'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0, 
	'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0, 
	'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0, 
	'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0, 
	'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0, 
	'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0, 
	'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0, 
	'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0, 
	'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0, 
	'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0} 

	TRANSLATIONS={'CYS': ['TGT', 'TGC'], 
	'ASP': ['GAT', 'GAC'], 
	'SER': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
	'GLN': ['CAA', 'CAG'], 
	'MET': ['ATG'], 
	'ASN': ['AAC', 'AAT'], 
	'PRO': ['CCT', 'CCG', 'CCA', 'CCC'], 
	'LYS': ['AAG', 'AAA'], 
	'STOP': ['TAG', 'TGA', 'TAA'], 
	'THR': ['ACC', 'ACA', 'ACG', 'ACT'], 
	'PHE': ['TTT', 'TTC'], 
	'ALA': ['GCA', 'GCC', 'GCG', 'GCT'], 
	'GLY': ['GGT', 'GGG', 'GGA', 'GGC'], 
	'ILE': ['ATC', 'ATA', 'ATT'], 
	'LEU': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
	'HIS': ['CAT', 'CAC'], 
	'ARG': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
	'TRP': ['TGG'], 
	'VAL': ['GTA', 'GTC', 'GTG', 'GTT'], 
	'GLU': ['GAG', 'GAA'], 
	'TYR': ['TAT', 'TAC'] }
	def __init__(self,seq=None,table=None, drop=False):
		'''
		CodonUsageTable(seq=<BioSeq object>,table=CodonUsageTable.CODON_DICT.copy(), drop=False)
		Creates codon table for a given Bio.seq object or takes precomputed codon dictionary.
		Either seq or table must be provided, and module will crash if dictionary has wrong structure.
		If table is provided, seq is not analysed.
                if drop is set to true, codons that contain symbols not from [ACTGactg] will cause ValueError exception
		Currently assumes seq to be DNA, RNA support to be added later'''
		#Notice that self.usage_table is added either here or in _count_int()
		if not table==None:
			if not isinstance(table,dict):
				raise TypeError ('Codon usage table must be a dictionary!')
			self.usage_table=table
		else:#No table provided, so let's compute it
			if seq==None:
				raise AttributeError ('Neither sequence nor table provided')#No seq either, nothing to do here
			if not isinstance(seq, Bio.Seq.Seq):
				raise TypeError ('Sequence must be Seq object from Bio.Seq.')
			self.dnaseq=str(seq)
			if self.dnaseq=='':
				raise ValueError ('Empty sequence object was used as input')
			self.dnaseq=self.dnaseq.upper()
			if len(self.dnaseq)%3>0:
				raise ValueError ('Sequence must contain whole number of codons, ie 3*n nucleotides.')
			self._count_int(drop)
		self.relative_usage=CodonUsageTable.CODON_DICT.copy()
		self._count_rel()

	def _count_int(self, drop):
		'''internal; counts codons in a string, not seqobject'''
		self.usage_table=CodonUsageTable.CODON_DICT.copy()
		for j in range(0, len(self.dnaseq),3):
			codon=self.dnaseq[j:j+3]
			if codon in self.usage_table:
				self.usage_table[codon]+=1
			elif drop:
				raise ValueError('Invalid codon: %s'%codon)

	def _count_rel(self):
		'''internal; computes relative frequencies assuming raw ones are ready'''
		for j in CodonUsageTable.TRANSLATIONS.keys():
			aa_total=0
			for k in CodonUsageTable.TRANSLATIONS[j]:
				aa_total+=self.usage_table[k]
			for k in CodonUsageTable.TRANSLATIONS[j]:
				if aa_total>0:self.relative_usage[k]=self.usage_table[k]/aa_total

	def __str__(self):
		'''Prints a 4-column TSV: AA - codon - raw count - relative frequency'''
		output=''
		for j in CodonUsageTable.TRANSLATIONS.keys():
			for k in CodonUsageTable.TRANSLATIONS[j]:
				output+='{0}\t{1}\t{2}\t{3}\n'.format(j,k,self.usage_table[k],self.relative_usage[k])
		return output

	def codon_count(self,codon):
		'''Tells how many times a given codon was seen by this object'''
		return self.usage_table[codon]


	def relative_count(self,codon):
		'''Tells how many times a given codon was seen by this object'''
		return self.relative_usage[codon]


	def __add__(self,right):
		'''Adds two codon usage tables: adds raw frequencies of codons
		and recounts relative ones'''
		#Add type checks later
		codons=CodonUsageTable.CODON_DICT.copy()
		for j in CodonUsageTable.CODON_DICT.keys():
			codons[j]=self.codon_count(j)+right.codon_count(j)
		c=CodonUsageTable(table=codons)
		return c
		

