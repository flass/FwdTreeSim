#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Models for simulation of phylogenetic tree(s) forward in time, specifically recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

from FwdTreeSim import IOsimul


genetypes = ['core', 'accessory-slow', 'orfan-slow', 'accessory-fast', 'orfan-fast']

replicontypes = ['chromosome', 'plasmid', 'chromid']

genecontentprofiles = {'chromosome':{'core':0.5, 'accessory-slow':0.15, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.15}, \
						  'plasmid':{'core':0.0, 'accessory-slow':0.35, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.35}, \
						  'chromid':{'core':0.3, 'accessory-slow':0.25, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.15} \
					  }

### not functional ; still under development

class Pangenome(object):
	"""Top simulation object, includes all the genomes and gene families"""
	def __init__(self, speciessim, lgenefams=[], t=0, **kwargs):
		verbose = kwargs.get('verbose', !kwargs.get('silent', False))
		# the species simulation at time t is taken to define the starting lineages in which to distribute gene families
		self.speciessim = speciessim
		self.t = t
		self.size = len(lgenefamsim) if lgenefamsim else kwargs.get('size', 10000)
		if lgenefams:
			self.genefams = lgenefams
		else:
			panprofile = kwargs.get('panprofile', IOsimul.MetaSimulProfile(profiles=[(p*self.size, IOsimul.DTLSimulProfile(type=t)) for t,p in genecontentprofiles['chromosome'].iteritems()]))
			self.genefams = []
			for n in range(self.size):
				# initiate a gene simulation
				bddtlmodel = models.BirthDeathDTLModel(self.size)
				bddtlprof = panprofile.sampleprofile(verbose=verbose)
				bddtlsim = simulators.DTLtreeSimulator(model=bddtlmodel, refsimul=moransim, profile=bddtlprof, noTrigger=True)
				# create a gene family
				genefam = GeneFam(name='g%06d'%n, type=bddtlprof.type, history=bddtlsim)
				self.genefams.append(genefam)
			

class Genome(object):
	"""the genome map of a single lineage"""
	def __init__(self, dreplicons):
		self.dreplicons = dreplicons
		
	def get_all_genes(self):
		lgenes = []
		for r in self.dreplicons.values():
			lgenes += r.lgenes
		return lgenes
		
class Replicon(object):
	""""""
	
	def __init__(self, lgenes, circular=True, **kwargs):
		self.lgenes = lgenes
		self.circular = circular
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "chromosome")
		
	def insert_gene(self, gene, locus=None):
		"""add gene to the ordered list of genes; by defautls at its end, or if specified at its index 'locus'"""
		if not locus is None: self.lgenes.insert(locus, gene)
		else: self.lgenes.append(gene)
		gene.attach_to_replicon(self)
		
	def remove_gene(self, gene, locus=None):
		"""remove gene from the ordered list of genes; by defautls use the object as its own identifier, but instead use its 'locus' index if specified"""
		if locus is None:
			g = gene
			self.lgenes.remove(g)
		else:
			g = self.lgenes[locus]
		g.attach_to_replicon(None)
		
	## all excise / integrate operations use the Python-style indexing with slice notation
	def integrate_segment(self, lsegment, k, orientation=1):
		"""insert a list of genes (genomic segment) into replicon's gene list at the given position i"""
		if orientation<0: lsegment.reverse()
		self.lgenes.insert(k, lsegment)
		
	def excise_segment(self, i, j):
		"""deletes a list of genes (genomic segment) from replicon's gene list between the given positions i and j, and returns the the segment
		
		j can be higher than i if the replicon is circular, which indicate that the excised segment covers the origin (indexes 0 and -1 [last]);
		if not, an error is thrown.
		"""
		if i==j: return []
		if i<=j:
			lsegment = self.lgenes[i:j]
			self.lgenes = self.lgenes[:i] + self.lgenes[j:]
		elif self.circular:
			# i > j
			lsegment = self.lgenes[i:] + self.lgenes[:j]
			self.lgenes = self.lgenes[j:i]
		else:
			raise IndexError, "begin point (i=%d) must be before the end point (j=%d) of excised segment on the linear replicon"%(i, j)
		return lsegment
		
	def translocate_segment(self, i, j, k, orientation=1):
		"""excise and re-integrate a genomic segment in the same replicon, i.e. excise from i-j position and integrate at position k (as defined before excision)."""
		#~ if orientation>0 and ((k==i) or (k==j)): pass # should be seemless
		lsegment = self.excise_segment(i, j)
		# find the relative coordinate of insetrtion after excision
		if i<=j:
			if k < i:     K = k
			elif k >= j:  K = k - j + i
			else: raise IndexError, "insertion point (k=%d) cannot be located between boundaries of the excised segment (i=%d, j=%d)"%(k, i, j)
		else:
			# i > j and self.circular==True (given excise_segment() has been executed before without an error)
			if k >= j: K = k - j
			else: raise IndexError, "insertion point (k=%d) cannot be located between boundaries of the excised segment (i=%d, j=%d trough the origin 0)"%(k, i, j)
		self.integrate_segment(lsegment, K, orientation=orientation)
			
	def transfer_integrate_segment(self, recrepli, i, j, k, orientation=1):
		"""excise a genomic segment from a replicon at position i-j and integrate it into another at position k."""
		lsegment = self.excise_segment(i, j)
		recrepli.integrate_segment(lsegment, K, orientation=orientation)
		
	def transfer_replace_segment(self, recrepli, i, j, k, l):
		"""excise a genomic segment from a replicon at position i-j and integrate it into another by replacing the local segment at position k-l."""
		lsegment = self.excise_segment(i, j)
		recrepli.excise_segment(k, l)
		if k <= l:
			recrepli.integrate_segment(lsegment, k)
		else:
			# the deletion spaned the origin so the scar is equated with the new origin
			recrepli.integrate_segment(lsegment, 0)
		
	
class Gene(object):
	""""""
	
	def __init__(self, treenode, **kwargs):
		self.treenode = treenode
		self.name = kwargs.get('name')
		self.replicon = None
		self.fam = kwargs.get('fam')
		
	def attach_to_replicon(self, repli):
		self.replicon = repli
		#replicon.insert_gene(self, locus=locus)

		
	
class GeneFam(object):
	""""""
	
	def __init__(self, treenode, **kwargs):
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "core")
		self.history = kwargs.get('genesim')
		
