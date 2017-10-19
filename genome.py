#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Models for simulation of phylogenetic tree(s) forward in time, specifically recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

from FwdTreeSim import IOsimul, models, simulators
import copy
import random
import numpy as np

genetypes = ['core', 'accessory-slow', 'orfan-slow', 'accessory-fast', 'orfan-fast']

replicontypes = ['chromosome', 'plasmid', 'chromid']

genecontentprofiles = {'chromosome':{'core':0.5, 'accessory-slow':0.15, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.15}, \
						  'plasmid':{'core':0.0, 'accessory-slow':0.35, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.35}, \
						  'chromid':{'core':0.3, 'accessory-slow':0.25, 'orfan-slow':0.1, 'accessory-fast':0.1, 'orfan-fast':0.15} \
					  }
					  
repliratecorrectors = {'chromosome':{'rdup':1,'rtrans':1,'rloss':1} 'plasmid':{'rdup':1,'rtrans':2,'rloss':2}, 'chromid':{'rdup':1,'rtrans':1.5,'rloss':1.5}}

explinitrepliprofiles = [('chromosome', 'c1', True, 8000, 1), ('plasmid', 'p1', True, 500, 0.1), ('plasmid', 'p2', True, 500, 0.1), ('plasmid', 'p3', True, 500, 0.1), ('plasmid', 'p4',True,  500, 0.1)]

### not functional ; still under development

class Pangenome(object):
	"""Top simulation object, includes all the genomes and gene families
	
	Initially, all genomes in the pangenome have the same syntenic structure, even though not the same gene complement
	(i.e. the gene family order is the same on the common genomic backbone, but only a fraction of gene families are present in each genomes).
	The same applies to replicons.
	
	"""
	def __init__(self, refsimul, lgenefams=[], t=0, **kwargs):
		verbose = kwargs.get('verbose', !kwargs.get('silent', False))
		# the species simulation at time t is taken to define the starting lineages in which to distribute gene families
		self.t = t
		self.refsimul = refsimul
		self.timeslices = refsimul.get_timeslices()
		self.lineages = refsimul.trees
		self.treeroot = refsimul.treeroot
		self.refconbran = refsimul.contempbranches
		self.popsize = refsimul.popsize
		
		connectlen = kwargs.get('connectlen', 0)
		
		# list of initial replicons' properties ; provide index for matrix below
		self.initrepliprofiles = kwargs.get('initrepliprofiles', explinitrepliprofiles)
		# presence/absence matrix for initial replicons (lineage roots, initial replicons)
		self.init_repli_phyloprofile = np.zeros((self.popsize, len(initrepliprofiles)), dtype=bool)
		# presence/absence matrix for gene families (lineage roots, gene families)
		self.init_gene_phyloprofile = np.zeros((self.popsize, sum([x[3] for x in initrepliprofiles)), dtype=bool)
		# genome-wide list of gene families; provide index for matrix above
		self.genefams = []
		# slice boundaries of the list of gene families above giving the replicon-specific gene family content
		self.repligenefams = []
		
		
		self.refsimul.labeltreenodes(onlyExtants=False)
		# attach empty genomes to the roots of lineages in reference species history
		for l, lineageroot in enumerate(self.lineages):
			lineageroot.genome = Genome([], name=lineageroot.label())
			
			
		# Initial set of replicons to be observed in this pangenome, defined as (type, circularity, name, number of distinct gene families, expected frequency at root)
		gg = -1 # genome-wide gene family index
		rlast = 0
		for r, initrepliprof in enumerate(self.initrepliprofiles):
			rtype, rname, rcirc, rsize, refreq = initrepliprof
			repligenefamslice = (rlast, rlast+rsize)
			self.repligenefams.append(repligenefamslice)
			rlast = rlast+rsize
			
			# determine the occurence profile of these replicons in species lineages
			#~ self.init_repli_phyloprofile[rname] = [n for n in range(refsimul.popsize) if random.random() < rfreq]
			for n in range(self.popsize)
				self.init_repli_phyloprofile[n, r] = random.random() < rfreq
			
			repliprofile = IOsimul.MetaSimulProfile(profiles=[(p, IOsimul.DTLSimulProfile(type=t)) for t,p in genecontentprofiles[rtype].iteritems()])
			for rg in range(rsize):  # replicon-specific gene family index
				gg += 1
				gfname = 'g%06d'%gg
				# initiate a gene simulation
				bddtlmodel = models.BirthDeathDTLModel(self.size)
				bddtlprof = repliprofile.sampleprofile(verbose=verbose)
				# do not pick gene lineages yet
				bddtlsim = simulators.DTLtreeSimulator(model=bddtlmodel, refsimul=moransim, profile=bddtlprof, noTrigger=True, genelineages=False)
				# generate gene presence distribution, but not the trees
				#~ genefam_phyloprofile = bddtlsim.pickgenelineages(return_index=True)
				self.init_gene_phyloprofile[bddtlsim.pickgenelineages(return_index=True),gg] = True
				# intersects the gene presence profile with the replicon presence profile
				#~ self.init_gene_phyloprofile[gfname] = list(set(genefam_phyloprofile) & set(self.init_repli_phyloprofile[rname]))
				#~ self.init_gene_phyloprofile[gfname].sort()
				self.init_gene_phyloprofile[,gg] = np.logical_and(self.init_gene_phyloprofile[,gg], self.init_repli_phyloprofile[,r])
				# now generates the lineage trees for the gene family
				bddtlsim.pickgenelineages(from_index=list(self.init_gene_phyloprofile[gfname].nonzero()[0]))
				# and connect those trees
				bddtlsim.connecttrees(l=connectlen)
				
				# create a gene family
				genefam = GeneFam(name=gfname, type=bddtlprof.type, history=bddtlsim)
				# add fam to general ledger of gene families
				self.genefams.append(genefam)
		
			
			# fill-in genomes attached to the roots of lineages in reference species history
			for l, lineageroot in enumerate(self.lineages):
				Gname = lineageroot.label()
				# create empty replicon
				repli = Replicon(lgene=[], type=rtype, name="%s_%s"%(Gname, rname), circular=rcirc)
				# collect gene families present in that genome
				gfis = self.init_gene_phyloprofile[l,].nonzero()[0]
				for gfi in gfis:
					if not (gfi >= repligenefamslice[0] and gfi < repligenefamslice[1]):
						# filter gene fams out of replicon
						continue
					# create gene
					gf = self.genefams[gfi]
					gene = Gene(treenode=gf.history.treeroot[Gname], name="%s_%s_%d"%(Gname, rname, g))
					# attach gene to replicon
					repli.insert_gene(gene)
				# attach replicon to genome
				lineageroot.genome.add_replicon(repli)
			
	def __iter__(self):
		"""iterator yielding replicons"""
		for genefam in self.genefams:
			yield genefam
			
	def branch_get_genome(self, cb, t):
		

class Genome(object):
	"""the genome map of a single lineage"""
	def __init__(self, lreplicons):
		self.lreplicons = lreplicons
		self.name = kwargs.get('name')
		
	def __iter__(self):
		"""iterator yielding replicons"""
		for repli in self.lreplicons:
			yield repli
			
	def __getitem__(self, key):
		for repli in self:
			if repli.name==key:
				return repli
		else:
			raise KeyError, "no replicon named '%s' in this genome instance"%key
		
	def get_all_genes(self):
		lgenes = []
		for repli in self:
			lgenes += repli.lgenes
		return lgenes
			
	def __len__(self):
		"""genome size in gene number"""
		return len(self.get_all_genes())
		
	def __copy__(self):
		
			
	def add_replicon(self, repli):
		"""add replicon due to large-scale transfer event"""
		self.lreplicons.append(repli)
		
	def cleanup_empty_replicon(self):
		if len(self)==0:
			raise ValueError, "Genome named '%s' reached a void state!"%self.name
		i = 0
		while i < len(self.lreplicons):
			if len(self.lreplicons[i])==0:
				del self.lreplicons[i]
			else:
				i += 1
		
				
		
class Replicon(object):
	"""Representation of the DNA macromolecules bearing genes in a bacterial genomes (replicon = replication unit)
	
	Replicons objects are part of a genome object, and are made as an orderred set (a list) of genes, 
	which can be rearranged and in/from which genes can be inserted/deleted.
	"""
	
	def __init__(self, lgenes, circular=True, **kwargs):
		self.lgenes = lgenes
		self.circular = circular
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "chromosome")
		self.ratecorrectors = kwargs.get('ratecorrectors', repliratecorrectors[self.type])
	
	def __iter__(self):
		"""iterator yielding genes"""
		for gene in self.lgene:
			yield gene
			
	def __len__(self):
		"""replicon size in gene number"""
		return len(self.lgenes)
		
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
	
	def __init__(self, genesim, **kwargs):
		self.name = kwargs.get('name')
		self.type = kwargs.get('type', "core")
		self.history = genesim
		
