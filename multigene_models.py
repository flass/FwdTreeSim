#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Models for simulation forward in time of evolution of a pangenome, specifically modifying gene trees relative to a species tree 
by recording co-events of duplication, transfer and loss of multiple genes."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "14 Feb 2017"

#~ import copy
import tree2
import random
from numpy.random import poisson, geometric, exponential, normal
from FwdTreeSim import models, IOsimul

def poissonmp1(lam, size=None):
	"""skewed Pisson distribution to avoid drawing 0"""
	return (poisson(lam=(lam-1), size=size) + 1)


class BlockBirthDeathDTLModel(models.BirthDeathDTLModel):
	"""Model considering the co-evolution of multiple gene families, notably through DTL events involving block of several neighbouring genes.
	
	Simulations of of the several (many) gene families making up the pangenome have to be made concurrently so co-events can be applied to all 
	concerned gene families at the same time. Incidentally, the model encompasses the evolution of all the gene families at once.
	
	Thus, events are sampled not byt test each branch at time t within one gene family like for BirthDeathDTLModel, but by sampling 
	among the genes (corresponding to a branch) as they exist in genomes at time t, with the same rates as in BirthDeathDTLModel, times the number of genes.
	This allows us to sample groups of neighbour genes to be deleted, duplicated or transferred together.
	
	Sampling of events is made as follows: 
	0) one current branch of the species phylogeny is picked at a rate e*N (e the event rate, N the number of gene in that branch's current genome)
	1) one gene is picked at a uniform rate within the genome (or within replicons). 
	2) then a gene track length is sampled from a given discreet ditribution (e.g. Poisson) to define the gene block involved in the co-event
	3) all these genes are marked for the event, that is then idependently applied to their respective history.
	
	the initial event sampling rate must thus:
	- be uniform for all genes (no distinction of core, accessory-slow, accessory-fast, etc. for the DTL rates)
	- take into account the number of genes in the genome
	- take into account the mean gene number expectation from the gene block length ditribution
	"""
	
	def __init__(self, **kwargs):
		print 'invoke multigene_models._BlockBirthDeathDTLModel.__init__()'
		super(BlockBirthDeathDTLModel, self).__init__(**kwargs)
		# self.__rseed
		# self.tunit
		# self.popsize
		# self.rdup
		# self.rtrans
		# self.rloss
		self.blocklenmean = kwargs.get('blocklenmean', 1.5)
		self.blocklendist = kwargs.get('blocklendist', poissonmp1)
		
	def stepforward(self, currbranches, pangenome, evtidgen=None, **kwargs):
		"""generate event record that refer to the reference tree"""
		# NB: assumes timeslice numbering starts with 1; t=0 would induce null-rate at first step
		t = pangenome.t
		currrefbranches = pangenome.refconbran[t]
		timeslice = pangenome.timeslices[t]
		levents = []
		devents = {}
		trec = {}
		# parallelization should happen here, iterating over branches - but risk of non idependence of branches with linking transfer events?
		# transfer event must have the donor as focus branch, thus limiting the tree manipulation on this branch to a single operation
		for cb in currbranches:
			# fetch the genome attached to branch and time point
			genome = pangenome.branch_get_genome(cb, t)
			
			#### adapt what's below:
			
			
			for gene in genome:
				# consider events exclusive, each proba is counted cummulatively so every draw in [0;1] can only point to one event type
				if random.random() <= self.rloss:
					evtype = 'loss'
					self.lossEvent(cb, timeslice)
				elif random.random() <= self.rloss+self.rdup:
					evtype = 'dupl'
					self.duplicationEvent(cb, timeslice)
				elif random.random() <= self.rloss+self.rdup+self.rtrans:
					evtype = 'trans'
					# pick a recipient branch from the current REFERENCE tree branches
					rec = random.choice(currrefbranches)
					self.transferEvent(cb, rec, timeslice)
				else:
					evtype = None
				if evtype: 
					e = DTLevent(evtype, cb, t, evtidgen, levents, devents, trec)
					# annoates the node's subtree labels by appending a string that signifies the event
					if e.recgenenode:
						if evtype=='trans':
							# add a pre-tag to salvage the info about the line of events leading to this node
							# which is not present on the label copied from the reference tree
							prelab = cb.label().split('-', 1)[-1]
						else:
							prelab = ""
						# annotate the transfer/duplication node
						cb.go_father().edit_label("%s%s%d"%(prelab, DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
						# and the recipient with all its descendants
						e.recgenenode.edit_all_labels("%s%s%d"%(prelab, DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
					else:
						# annotate the node with the event (only for losses)
						cb.edit_label("%s%d"%(DTLevent.etshorts[evtype], e.id()), mode='a', sep="-")
						simul.extincts.append(cb)
		return (levents, devents, trec)
