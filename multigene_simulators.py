#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Simulate forward in time the evolution of a pangenome as a collection of genome maps and phylogenetic trees, 
recording in gene trees (relative to a species tree) the (co-)events of duplication, transfer and loss of (multiple) genes.
"""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "15 Feb 2017"

from FwdTreeSim import simulators, multigene_models, IOsimul

class DTLgenomeSimulator(simulators.DTLtreeSimulator):
		""""""
	def __init__(self, model, pangenome, noTrigger=True, **kwargs):
		#~ super(DTLgenomeSimulator, self).__init__(model=model, noTrigger=noTrigger, genelineages=False, **kwargs) 
		
		## from BaseTreeSimulator:
		# self.model = model
		# self.t = 0
		# self.times = kwargs.get('times', [])
		# self.eventsrecord = kwargs.get('eventsrecord', [])
		# self.eventsmap = kwargs.get('eventsmap', {})
		# self.extincts = kwargs.get('extincts', [])
		# self.contempbranches =  kwargs.get('contempbranches', [])
		# self.nodeClass = kwargs.get('nodeClass', tree2.Node)
		# self.ngen = kwargs.get('ngen')
		# self.eventidgen
		# self.dumppickle_warning
		
		## from MultipleTreeSimulator:
		# self.popsize = self.model.popsize
		# self.profile = kwargs.get('profile', IOsimul.SimulProfile())
		
		## from DTLtreeSimulator:
		# self.refsimul = refsimul	# keep hard link to reference/species tree simulator object; better pickle them together!
		# self.refroot = getattr(refsimul, 'treeroot', None)
		# self.reftrees = refsimul.trees
		# self.refconbran = refsimul.contempbranches
		# self.reftimeslices = refsimul.get_timeslices()
		# self.refextincts = refsimul.extincts
		# self.ngen = refsimul.t - 1
		# self.times = refsimul.times
		# self.eventsmap = {}
		# self.transferrec = {}
		
		self.pangenome = pangenome
		
		#~ for genefam in self.pangenome:
