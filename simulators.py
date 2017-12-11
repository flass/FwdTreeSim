#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Simulate phylogenetic tree(s) forward in time, notably under models recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import os
import copy
import random
from FwdTreeSim import models, IOsimul
from FwdTreeSim import nodelabelprefix
import tree2
from numpy.random import gamma

### utilitary functions

def _connecttrees(trees, l=0, returnCopy=False):
	"""create a common root for all trees in the given list, with branches of lenght l."""
	treeroot = trees[0].newnode()
	for tree in trees:
		if returnCopy: t = copy.deepcopy(tree.go_root())
		else: t = tree.go_root()
		treeroot.link_child(t, newlen=tree.lg()+l)
	return treeroot
	
def _get_timeslices(t, times):
	ts = []
	for t in range(t):
		low = 0 if t==0 else ts[t-1][1] # t low bound is t-1 up bound
		up = low + times[t]
		ts.append((low, up))
	return ts
	

# for dumppickle() funciton; special case of simulators.BaseTreeSimulator class
dumpwarningmsg =  "cannot pickle generator objects, have to delete event id generator 'self.eventidgen' first\n"
dumpwarningmsg += "(will discontinue numeration of events for further simulation).\n"
dumpwarningmsg += "Delete generator before pickling? (y/n) "
checkdic = {'simulators.BaseTreeSimulator':{'attrname':'eventidgen', 'prompt':dumpwarningmsg}}

######################################
# Tree Simulators
# wrappers for the simulation models; perfom the simulation from t0 to the end
######################################


class BaseTreeSimulator(object):
	"""
	Base wrapper class for forward tree simulations.
	
	Main attributes:
	- self.model           : a model instance from *Model classes from 'models' module.
	- self.tree            : a (list of) tree2.Node (or derived class) phylogenetic tree object; 
	or self.trees             uses the tree provided through 'starttree' argument or (default) a single-node tree (set) as a seed for the simulation.
	- self.t               : the number of the last simulation iteration / time slice.
	
	Attributes than can be provided at instance intitiation through keyword arguments (e.g. to resume previous simulation)
	- self.times           : list of evolutionary time spent within each time slice.
	- self.extincts        : list of leaf nodes representing lineages considered extinct.
	- self.events		   : record of all events
	                          structured as a dictionary (which keys are event ids) of event objects (models.BaseEvent children classes).
	- self.eventsrecord    : record events by time slice where they occurred;
	                          structured as a list (which indexes correspond to the time slice number) of lists of event ids.
	- self.eventsmap       : record events by tree branch on which they occurred (in case of transfer events, donor node is used);
	                          structured as a dictionary (which keys are reference tree node labels) of lists of event ids.
	- self.contempbranches : record by time slice of the extant branches present at this time that were subject to the evolutionary processes;
	                          structured as a list of list of of node objects.
	
	Falcultative keyword arguments:
	- ngen                 : the end number of generation to simulate (count includes potential previous iterations).
	- noTrigger            : when initiating an instance, automatic checkdata() and evolve(ngen) triggers from parent classes
	                          are deactivated when noTrigger=True.
	"""
	
	def __init__(self, model, **kwargs):
		print 'invoke _BaseTreeSimulator__init__'
		#~ print 'kwargs:', kwargs
		self.model = model
		self.t = 0
		self.times = kwargs.get('times', [])
		self.events = kwargs.get('events', {})
		self.eventsrecord = kwargs.get('eventsrecord', [])
		self.eventsmap = kwargs.get('eventsmap', {})
		self.extincts = kwargs.get('extincts', [])
		self.contempbranches =  kwargs.get('contempbranches', [])
		self.nodeClass = kwargs.get('nodeClass', tree2.AnnotatedNode)
		self.nodeAttr = kwargs.get('nodeAttr', [])+['extinct']
		self.ngen = kwargs.get('ngen')
		self.eventidgen = models.eventIdGen()	# generator object providing serial unique identifiers for event objects
		self.dumppickle_warning = '\n'.join(["Cannot pickle generator objects; to pickle simulator object instance %s,"%repr(self), \
											"you have to first delete its event id generator attribute 'eventidgen'", \
											"(will discontinue numeration of events for further simulation)."])
		self.logger = IOsimul.SimulLogger(simultype=type(self))
		
	def __getattr__(self, name):
		if name=='__dict__': return self.__dict__
		try:
			# search normal attributes
			return self.__dict__[name]
		except KeyError:
			try:
				# search private cached attirbutes
				return self.__dict__['_'+name]
			except KeyError:
				try:
					# look up an attribute-generating method
					return type(self).__dict__['get_'+name](self)
				except KeyError:
					raise AttributeError, "'%s' is not an attribute of %s instance of %s"%(name, repr(self), type(self))
		
	def checkdata(self):
		"""placeholder; method should be called from descendant classes"""
		for attr in ['times', 'eventsrecord', 'extincts', 'contempbranches']:
			assert len(getattr(self, attr)) == self.t
		return 0
		
	def get_timeslices(self):
		return _get_timeslices(self.t, self.times)
		
	def update_records(self, levents, dnode2eventids, contempbranches, brlen=None):
		leventids = []
		for event in levents:
			eid = event.evtid()
			self.events[eid] = event
			leventids.append(eid)
		self.eventsrecord.append(leventids)
		self.eventsmap.update(dnode2eventids)
		self.contempbranches.append(contempbranches)
		if brlen: self.times.append(float(brlen))
		
		
	@staticmethod
	def recordCollapsedNodeLabs(leaf, dcollapsednodelabs):
		"""assume pruning always uses Node.pop.newpop as root node cannot go extinct"""
		flab = leaf.go_father().label()
		blab = leaf.go_brother().label()
		dcollapsednodelabs[flab] = (blab, 'c') # record it was flab's child
		dcollapsednodelabs[leaflab] = (blab, 'b') # record it was leaflab's brother
		
	@staticmethod
	def copy_prune_dead_lineages(tree, extincts, collapsenodes=False, trimroot=False):
		"""removed branches leading to extinct lineages; require unique naming of all nodes!!!"""
		livetree = tree.deepcopybelow(shallow_copy_attr=['ref', 'event'], keep_lg=True)
		lleaves = livetree.get_leaf_labels()
		for leaf in extincts:
			leaflab = leaf.label()
			if leaflab in lleaves:
				if len(lleaves)==1:
					# last leaf to stand on the tree, cannot pop it ; rather return null result
					return None
				livetree.pop(leaflab, noCollapse=(not collapsenodes))
				lleaves.remove(leaflab)
		if trimroot:
			# set root branch lenght to zero
			livetree.set_lg(0)
		return livetree
		
	def get_tree(self, connecttrees=0):
		"""return simulated tree (collection) as a single tree
		
		i.e. if there are multiple trees, a single connecting root to all trees is enforced
		"""
		if hasattr(self, 'trees'):
			self.connecttrees(l=connecttrees, returnCopy=False)
			return self.treeroot
		elif hasattr(self, 'tree'):
			return self.tree
		
	def get_nodes_with_descendants(self, sampled=[]):
		tree = self.get_tree()
		lwithdescent = []
		sext = set(self.extincts)
		for node in tree:
			if (set(node.get_leaf_labels()) - sext):
				# some species below this nodes are extant
				lwithdescent.append(node.label())
		lwithdescent =  [l for l in lwithdescent if (not l in sampled)] # preserves the order ; not sure why i do remove these nodesin the first place...
		return lwithdescent
	
	def prepare_write_endlog(self, connecttrees=10):
		"""prepare for data export of simulation (a posteriori)"""
		tree = self.get_tree(connecttrees=connecttrees)
		tree.check_unique_labelling()
		dcollapsednodelabs = {}
		extanttree = self.get_extanttree(lentoroot=connecttrees, dcollapsednodelabs=dcollapsednodelabs, removelosses=False)
		# filter out the list of events
		
		if self.nodeClass == tree2.AnnotatedNode:
			# number nodes by increasing branch length distance from root, as ALE labels a species tree
			extanttree.complete_node_ids(order=[sorted(extanttree.get_all_children(), key=lambda x: x.distance_root())])
	
	def write_endlog(self, dirout, simultype):
		"""data export of simulation (a posteriori)"""
		logger = IOsimul.SimulLogger(simultype=simultype, bnfout="%/log_"%dirout)
		logger
			

################################################
# Classes describing how wether a single or many
# trees are simulated within the population
################################################
		
class SingleTreeSimulator(BaseTreeSimulator):
	def __init__(self, model, starttree=None, **kwargs):
		super(SingleTreeSimulator, self).__init__(model=model, **kwargs)
		self.tree = starttree or self.nodeClass(l=float(0), lab="Root")
		
		if not kwargs.get('noTrigger'):
			self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(ngen=self.ngen)
		
	def checkdata(self):
		"""assert (self-)consistency of data, i.e. support of tree object class and that nodes in eventsrecord are included in the tree"""
		assert isinstance(self.tree, tree2.Node) # any descendant class of tree2.Node
		lerrnodes = []
		for ex in self.extincts:
			if ex not in self.tree:
				#~ raise TreeReferenceError(ex, self.tree)
				lerrnodes.append(ex)
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.tree)
		for t in self.eventsrecord:
			evt = self.eventsrecord[t]
			for u in evt:
				evtu = evt[u]
				for v in evtu:
					if evtu[v] not in self.tree:
						#~ raise TreeReferenceError(evtu[v], self.tree)
						lerrnodes.append(evtu[v])
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.tree)
		return 0

	def evolve(self, ngen=None, verbose=False, nodeathspan=[], stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		if ngen is None:
			if self.ngen >= 0: ngener = self.ngen
			else: raise ValueError, "must provide an integer non-negative value for 'self.ngen' attribute (got %s)"%repr(self.ngen)
		else:
			ngener = ngen
		while self.t < ngener:
			self.t += 1
			levents, dnode2eventids, contempbranches, brlen = self.model.stepforward(self, allowdeath=(self.t not in nodeathspan), evtidgen=evtidgen)
			# record what happened
			self.update_records(levents, dnode2eventids, contempbranches, brlen)
			if verbose: self.verbevolve(dnode2eventids)
			# check if tree got full extinct
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
			elif set(self.tree.get_leaves())==set(self.extincts):
				print "Extinction of all lineages at time", self.t
				rval = 1
				break
		else:
			# set root branch lenght to zero
			self.tree.set_lg(0)
			# destroys cache outdated attributes
			if hasattr(self, '_extanttree'): self._extanttree = None
			rval = 0
		return rval
	
	def verbevolve(self, devents):
		print "\ntime:", self.t,
		nextanttm1 = len(self.contempbranches[-1])
		#~ if self.t > 1:
			#~ nextanttm1 = len(self.contempbranches[self.t-2])	# aims at timeslice t-1, but shifted due to Python 0-based numbering ;
			#~ nextanttm1 = len(self.contempbranches[-1])	# aims at timeslice t-1, but shifted due to Python 0-based numbering ;
			#~ # equivalent to self.contempbranches[-1], but only if we keep calling this function above the update of contempbranches in evolve()
		#~ else:
			#~ nextanttm1 = 1
		print "\tstarted with %d live branches"%(nextanttm1)
		if (isinstance(self.model, models.UniformDiscreetBirthDeathModel) or isinstance(self.model, models.GenericDiscreetBirthDeathModel)):
			for evtype in devents:
				netgrowth = sum(node.nb_children() for node in devents[evtype]) - len(devents[evtype])
				print "  %s:\t%d\t(%.2g per lineage)"%(evtype, netgrowth, float(netgrowth)/nextanttm1)
		elif isinstance(self.model, PartialMoranProcess):
			for evtype in devents:
				print "\t%s:\t%d"%(evtype, sum(int(not (type(e) is int)) for e in devents[evtype])),
	
	################################	
	## post-simulation tree-cleaning
	
	def labeltreenodes(self, dictprefix=nodelabelprefix, onlyExtants=True):
		"""puts distinctive labes at internal nodes and extinct and extant leaves (default prefixes are N, E, S,respectively). Labelling follows increasing order from root"""
		self.tree.complete_internal_labels(prefix=dictprefix['livetip'], onlyLeaves=True, exclude=self.extincts)
		if not onlyExtants:
			self.tree.complete_internal_labels(prefix=dictprefix['deadtip'], onlyLeaves=True)
			self.tree.complete_internal_labels(prefix=dictprefix['node'], excludeLeaves=True)
		
	def scaletree(self, relToExtantMrca=True):
		if relToExtantMrca: t = self.get_extanttree(compute=True)
		else: t = self.tree
		root2tipdist = self.get_extants()[0].distance_root()
		self.tree /= root2tipdist
		self.extanttree /= root2tipdist
		self.timeslices = [float(ts)/root2tipdist for ts in self.timeslices]
		
		
	#############################
	## result description methods
			
	def get_extanttree(self, compute=True, dcollapsednodelabs=None, **kw):
		"""returns phylogenetic tree 'cleaned' of its dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttree
		# only works with all nodes being labelled, to use labels as references rather than the node objects (which refer to the original tree)
		self.labeltreenodes()
		livetree = self.copy_prune_dead_lineages(self.tree, self.extincts, dcollapsednodelabs=dcollapsednodelabs)
		self._extanttree = livetree	# create or update cache attribute
		return livetree
				
	def nb_extant(self):
		return self.tree.nb_leaves() - len(self.extincts)
		
	def get_extants(self, depthsorted=False):
		#~ e = list(set(self.tree.get_leaves()) - set(self.extincts))
		e = [leaf for leaf in self.tree.get_leaves() if not (leaf.extinct is True)]
		if depthsorted: e.sort(key=lambda x: x.depth())
		return e
	
class MultipleTreeSimulator(BaseTreeSimulator):
	def __init__(self, model=None, **kwargs):
		print 'invoke _MultipleTreeSimulator__init__'
		#~ print 'kwargs:', kwargs
		super(MultipleTreeSimulator, self).__init__(model=model, **kwargs)
		self.profile = kwargs.get('profile', IOsimul.SimulProfile())
		if not self.model:
			# by default create a private model instance; these are better not shared between gene families 
			# as their rates might have to be updated at different time points
			self.model = eval('models.'+self.profile.modeltype)(**kwargs)
		self.popsize = self.model.popsize
		self.trees = [self.nodeClass(l=float(0), lab="Root_%d"%i, addAttr=self.nodeAttr) for i in range(self.popsize)]
		
		if not kwargs.get('noTrigger'):
			self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(ngen=self.ngen)
		
	def checkdata(self):
		"""assert (self-)consistency of data, i.e. support of tree object class and that nodes in eventsrecord are included in the trees"""
		for i in range(len(self.trees)):
			if not isinstance(self.trees[i], tree2.Node): raise TypeError, "element %d of 'trees' attribute is not a tree2.Node instance"%(i)
		lerrnodes = []
		allnodes = sum((t.get_all_children() for t in self.trees), [])
		for ex in self.extincts:
			if ex not in allnodes:
				lerrnodes.append(ex)
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.trees)
		for t in self.eventsrecord:
			evt = self.eventsrecord[t]
			for u in evt:
				evtu = evt[u]
				for v in evtu:
					if evtu[v] not in allnodes:
						lerrnodes.append(evtu[v])
		if lerrnodes: raise tree2.AggregateTreeReferenceError(lerrnodes, self.trees)
		return 0
		
		
	def get_extants(self, depthsorted=False):
		"""generate a single list of all extant leaves across all trees in  the self.trees set"""
		#~ sleave = set(sum((tree.get_leaves() for tree in self.trees), []))
		#~ extants = list(sleave - set(self.extincts))
		extants = sum(([leaf for leaf in self.tree.get_leaves() if not (leaf.extinct is True)] for tree in self.trees), [])
		if depthsorted: extants.sort(key=lambda x: x.depth())
		if isinstance(self.model, models.MoranProcess): assert len(extants) == self.popsize
		return extants
		
	def update_rates(self):
		"""update model according to the profile's rate schedule and the current time step"""
		if self.t in self.profile.rateschedule:
			# updates the rate attributes, e.g. 'rdup', 'rtrans', 'rloss'
			self.model.__dict__.update(self.profile.rateschedule[self.t])
		
	def evolve(self, ngen=None, verbose=False, nodeathspan=[], stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		if ngen is None:
			if self.ngen >= 0: ngener = self.ngen
			else: raise ValueError, "must provide an integer non-negative value for 'self.ngen' attribute (got %s)"%repr(self.ngen)
		else:
			ngener = ngen
		while self.t < ngener:
			self.t += 1
			# check if need to update model at this time step
			self.update_rates()
			levents, dnode2eventids, contempbranches, brlen = self.model.stepforward(self, allowdeath=(self.t not in nodeathspan), evtidgen=evtidgen)
			# record what happened
			self.update_records(levents, dnode2eventids, contempbranches, brlen)
			if verbose: self.verbevolve(devents)
			# check if simulation should stop, e.g. because tree got full extinct
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
		else:
			# destroys outdated cache attributes
			if hasattr(self, '_extanttrees'): self._extanttrees = None
			rval = 0
		# give labels to extant tips
		for i, extant in enumerate(self.get_extants()): extant.edit_label("%s%d"%(nodelabelprefix['livetip'], i))
		return rval
		
	#~ def get_most_extant_tree(self):
		#~ """choose the tree with most extant leaves for the focal output of simulation"""
		#~ pass
	
	def verbevolve(self, devents):
		print "\ntime:", self.t
		
	def connecttrees(self, l=0, returnCopy=False, treeattr='trees'):
		"""create a common root for all trees in the population, with branches of lenght l
		
		store the new single tree object in new attirbute 'treeroot'.
		"""
		l = float(l)
		assert l>=0
		trees = getattr(self, treeattr)
		if returnCopy:
			treeroot = _connecttrees(trees, l=l, returnCopy=returnCopy)
		elif not hasattr(self, 'treeroot'): 
			treeroot = _connecttrees(trees, l=l, returnCopy=returnCopy)
			self.treeroot = treeroot
			if treeattr=='genetrees' and hasattr(self, 'refroot'):
				# for the moment means isinstance(self, DTLTreeSimulator) is True,
				# but could be extanded to other classes with gene trees...
				self.treeroot.ref = self.refroot
			if l>0:
				# update self.times atribute
				for i, ti in enumerate(self.times):
					self.times[i] = ti + l
		else: 
			treeroot = self.treeroot
		return treeroot
	
	################################	
	## post-simulation tree grooming
	
	def labeltreenodes(self, dictprefix=nodelabelprefix, treesattrname='trees', onlyExtants=True, silent=True):
		"""puts distinctive labels at internal nodes and extinct and extant leaves (default prefixes are N, E, S,respectively). Labelling follows increasing order from root"""
		for t in getattr(self, treesattrname):
			t.complete_internal_labels(prefix=dictprefix['livetip'], onlyLeaves=True, exclude=self.extincts, silent=silent)
			if not onlyExtants:
				t.complete_internal_labels(prefix=dictprefix['deadtip'], onlyLeaves=True, silent=silent)
				t.complete_internal_labels(prefix=dictprefix['node'], excludeLeaves=True, silent=silent)
				
	def shadetreenodes(self, shadeeventstoextincts=True, silent=True, dimfactor=0.66):
		"""puts distinctive colors saturation on branches leading to extinct and extant leaves have different saturation"""
		def dimcolor(node):
			if hasattr(node, 'color'):
				if node.color():
					# color already set by an event
					if shadeeventstoextincts:
						node.edit_color([c*dimfactor for c in node.color()])
				else:
					 node.edit_color([255*dimfactor, 255*dimfactor, 255*dimfactor])
		
		#~ sextincts = set(self.extincts)
		for node in self.extincts:
			dimcolor(node)
			# below unecessary unless single-child nodes (rosary-like branches) are present
			#~ fat = node.go_father()
			#~ while fat:
				#~ ch = fat.get_leaves()
				#~ if (set(ch) - sextincts):
					#~ # father has non-extinct children
					#~ break
				#~ dimcolor(fat)
				#~ node = fat
				#~ fat = node.go_father()
			
	def get_extanttrees(self, compute=True, addextincts=[], collapsenodes=False, cache=True):
		"""returns a COPY of the list of all trees 'cleaned' of their dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttrees
		# only works with all nodes being labelled, to use labels as references rather than the node objects (which refer to the original tree) in Node.pop() 
		self.labeltreenodes(silent=False) # should not be necessary though
		# exclude null trees (None objects) from the returned list
		extanttrees = [self.copy_prune_dead_lineages(t, self.extincts+addextincts, trimroot=False, collapsenodes=collapsenodes) for t in self.trees if t]
		if cache: self._extanttrees = extanttrees
		return extanttrees

	def get_extanttree(self, compute=True, addextincts=[], lentoroot=0, cache=True, collapsenodes=False, **kw):
		"""wrapper for get_extanttrees() methods; returns a COPY of the single root connecting all trees 'cleaned' of their dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttree
		else:
			extanttrees = [t for t in self.get_extanttrees(compute=True, addextincts=addextincts, collapsenodes=collapsenodes, **kw) if t] # possibly use a child class' method
			if extanttrees: # not an empty list
				extanttree = _connecttrees(extanttrees, l=lentoroot, returnCopy=False)
				if cache: self._extanttrees = extanttrees
			else: extanttree = None
			return extanttree
		
	
	
class DTLtreeSimulator(MultipleTreeSimulator):
	"""Simulates gene trees from (species) reference trees under a DTL model (see Szollosi et al. 2013, Lateral Gene Transfer from the Dead, Systematic Biology 62(3):386–397)
	
	Attribute 'genetrees' (these gene trees could as well be a locus/plasmid tree as well as a gene tree) is generated orignially by deep-copying the reference tree set, 
	sampling a reference tree at proba 'rootfreq' (to model the heterogeneous ocurrence of accessory genes in genomes of a prokaryotic population).
	A reference to the the original node is kept for each node the gene tree copies under the additional attribute 'ref'.
	During the simulation proceess, modifications are made on gene trees, which can be: copying a subtree from itself and grafting it next to the original (duplication), 
	removing a subtree (loss), or removing a subtree and replacing it with a copy of a subtree __from the reference tree set__ (transfer).
	
	As for parents classes, events are recorded in the 'eventsrecord' and 'eventsmap' attributes, storing events object indexed by time slice and reference tree node origin, respectively. 
	Additional atrributes are:
	- 'transferrec' a dictionary of transfer event ids indexed by reference tree node destination.
	- 'extantevents' a list of DTL events that are received (occur) in species lineages with (sampled) descendants; only those will be listed in the pseudo-ALE output.
	- 'refnodeswithdescent' the list of species nodes with (sampled) descendants; this is computed at the start of the simulation 
	   or can be provided through argument 'refnodeswithdescent', in which case one can restrict it to a sublist of sampled species.
	
	Profile of a gene family can be passed on with keywor argument 'profile', determining the root frequency of the gene in the population 
	and the rates of evolution (edits the evolution model object), in a time-heterogeneous manner.
	"""
	def __init__(self, model=None, noTrigger=False, **kwargs):
		print 'invoke _DTLtreeSimulator__init__'
		#~ print 'kwargs:', kwargs
		super(DTLtreeSimulator, self).__init__(model=model, noTrigger=True, allow_multiple=kwargs.get('allow_multiple',False), **kwargs) ##W : dynamic allow_multiple rather than static False
		# automatic checkdata() and evolve(ngen) triggers from parent classses are deactivated by noTrigger=True
		#~ self.rootfreq = kwargs['rootfreq'] # frequency a which the gene family is found at the root of each tree in the multiple reference tree set (species lineages from a Moran process)
		refsimul = kwargs.get('refsimul')
		if refsimul:
			self.refsimul = refsimul	# keep hard link to reference/species tree simulator object; better pickle them together!
			self.refroot = refsimul.connecttrees(returnCopy=False)	# enforce existence of a root in reference simulation
			self.reftrees = refsimul.trees
			self.refconbran = refsimul.contempbranches
			self.refextincts = refsimul.extincts
			self.ngen = refsimul.t - 1
			self.times = refsimul.times
		else:
			self.reftrees = kwargs['reftrees']
			self.refconbran = kwargs['refcontempbranches']
			self.refextincts = kwargs['refextincts']
			self.times = kwargs['times']
			self.ngen = kwargs.get('ngen')
			
		self.reftimeslices = _get_timeslices(len(self.times), self.times)
		
		# generate a list of nodes of the reference tree 
		self.refnodeswithdescent = kwargs.get('refnodeswithdescent', refsimul.get_nodes_with_descendants())
		self.transferrec = {}
		self.extantevents = []
		
		
		self.profile = kwargs.get('profile', IOsimul.DTLSimulProfile(type='core'))
		self.genetrees = []
		self.pickgenelineages( allow_multiple=kwargs.get('allow_multiple',False) ) ## W : go search in kwargs rather than direct call to unknown
		
		if not noTrigger:
			# self.checkdata()
			# if ngen specified, launch simulation for ngen iterations
			if self.ngen: self.evolve(self.ngen, connecttrees=kwargs.get('connecttrees', False))
			
	def pickgenelineages(self, allow_multiple=False, return_index=False, from_index=None):
			"""with a proba rootfreq, copy a gene tree from each reference tree
			
			# gene trees are like the reference tree set, with an addition of a link of each gene tree node to its reference tree node	
			
			# when rootfreq > 1, do several rounds of sampling the reference trees, 1 per integer slice of the expected frequency,
			# at a probability that's first 1 and in the last round is the decimal remainder
			
			# multiple stable copies should be modelled as different gene families (duplicate lineage imply roblems in node labelling for instance);
			# the only benefit of modelling them as an originally multi-copy gene familly is that they can recombine...
			# not sure this is really much different from evolving separately though, and would complicate locus attribution at t=0 in the model with linkage
			# by default, should keep starting with <=1 copies/genome
			"""
			genetrees = []
			if self.profile.rootfreq > 1:
				if allow_multiple: print "Warning: 'allow_multiple' option is on, and expected root frequency of gene is > 1  => !!experimental!!"
				else:  raise ValueError, "Expected root frequency of gene should be <= 1 (for rootfreq > 1, must turn on 'allow_multiple' option => !!experimental!!)"
			p = float(self.profile.rootfreq)
			l = float(self.profile.rootlen)
			while p >= 0:
				if return_index: 
					genetrees += [n for n in range(len(self.reftrees)) if random.random() < p]
				else:
					if from_index:
						for n in from_index:
							genetrees.append(self.reftrees[n].deepcopybelow(keep_lg=True, add_ref_attr=True))
					else:
						genetrees += [rt.deepcopybelow(keep_lg=True, add_ref_attr=True) for rt in self.reftrees if random.random() < p]
				if genetrees:
					# if not tree was sampled, redo a sampling round: there is no point simulating an entirely absent gene family
					p -= 1.0
				# not ideal as that does not introduce much variance in copy number, as only the nth copy is not certain to be sampled.
				# Ideal for modelling genes with proba 1.x, i.e. core genes with possibility of transient duplicates arising
			if return_index: return genetrees
			else: self.genetrees = genetrees
		
	def get_current_branches(self, l):
		"""retrieves branches reaching the root-to-tip length of l"""
		currbranches = []
		for tree in self.genetrees:
			# beware! after evolution, the node 'tree' stored in the 'self.genetrees' list may not be the root of the tree, due to addition of nodes for transfer or duplication event
			currbranches += tree.go_root().get_comtemporary_branches(l)
		return currbranches
		
	def evolve(self, ngen=None, verbose=False, stopcondition=(lambda x: (None, None)), **kwargs):
		"""simulation engine, iterates step of the simulation model"""
		evtidgen = kwargs.get('evtidgen', self.eventidgen)
		if ngen is None:
			if self.ngen >= 0: ngener = self.ngen
			else: raise ValueError, "must provide an integer non-negative value for 'self.ngen' attribute (got %s)"%repr(self.ngen)
		else:
			ngener = ngen
		while self.t < ngener:
			# record time step advance
			self.t += 1
			# check if need to update model at this time step
			self.update_rates()
			print "t=%d"%self.t, 
			timeslice = self.reftimeslices[self.t]
			# get current bene tree branches
			currbranches = self.get_current_branches((timeslice[1]+timeslice[0])/2) # input
			#~ print ": currbranches", [cb.label() for cb in currbranches]
			levents, dnode2eventids, trec = self.model.stepforward(currbranches, self, evtidgen=evtidgen)
			for evt in levents:
				# record the events sent/received by surviving lineages
				for refnode in [evt.recrefnode, evt.donrefnode]:
					if refnode in self.refnodeswithdescent:
						self.extantevents.append(evt.evtid())
						break # for refnode loop
				if self.logger:
					# record on log what happened
					# write out all events
					self.logger.DTLsingleEventLog(evt) ## W : replaced :  DTLEventLog -> DTLsingleEventLog
			self.update_records(levents, dnode2eventids, currbranches)
			self.transferrec.update(trec)
			# check if simulation should stop
			s, v = stopcondition(self)
			if s:
				print "Stop at time %d because: %s"%(self.t, v)
				rval = 0
				break
		else:
			rval = 0
		return rval
		
	def finish(self, **kwargs):
		if kwargs.get('connecttrees', False):
			treeroot = self.connecttrees(l=self.profile.rootlen, returnCopy=False)
		multreelen = kwargs.get('multreelen', getattr(self.profile, 'multreelen', 0))
		if multreelen:
			if type(multreelen) is tuple:
				# implemented in the profile, format being: ('functionname', param1, pram2, ...)
				m = eval(multreelen[0])(*multreelen[1:])
			else:
				m = float(multreelen)
			if multreelen<=0: raise ValueError, "unvalid tree length multiplier value: %f; should have a value > 0"%(m)
			if m!=1:
				treeroot *= m
				# update self.times atribute
				for i, ti in enumerate(self.times):
					self.times[i] = ti * m
		return treeroot
		
	def connecttrees(self, l=0, returnCopy=False):
		"""create a common root for all gene trees in the population, with branches of lenght l
		
		store the new single tree object in new attirbute 'treeroot'.
		gene tree-specific method ads the ref attibute to specie tree root.
		"""
		return super(DTLtreeSimulator, self).connecttrees(l=l, returnCopy=returnCopy, treeattr='genetrees')
		
	def labeltreenodes(self, dictprefix=nodelabelprefix, onlyExtants=True, silent=True):
		super(DTLtreeSimulator, self).labeltreenodes(dictprefix=dictprefix, treesattrname='genetrees', onlyExtants=onlyExtants, silent=silent)

	def copy_prune_dead_lineages(self, tree, collapsenodes=False, trimroot=False):
		"""remove lineage ends that belong to extinct reference species (not loss events!); require unique labeling of all nodes!!!"""
		if not tree: return None
		livetree = tree.deepcopybelow(shallow_copy_attr=['ref', 'event'], keep_lg=True)
		lleaves = livetree.get_leaf_labels()
		drefleaves = {}
		for leaf in livetree.get_leaves():
			#~ drefleaves.setdefault(leaf.ref.label(), []).append(leaf.label())
			# take only the first component of leaf label to collect the leaves in extinct species
			# that descend from a transfered/duplicated lineage AND APPARENTLY HAVE LOST REFERENCE TO OWNER SPECIES
			### DIRTY WORKAROUND, HAVE TO FIX THE LOSS OF REFERENCE !!!!!
			drefleaves.setdefault(leaf.label().split('-')[0], []).append(leaf.label())
		for extref in self.refextincts:
			if extref.label() in drefleaves:
				for leaflab in drefleaves[extref.label()]:
					if len(lleaves)==1:
						# last leaf to stand on the tree, cannot pop it ; rather return null result
						return None
					livetree.pop(leaflab, noCollapse=(not collapsenodes))
					lleaves.remove(leaflab)
		if trimroot:
			# set root branch lenght to zero
			livetree.set_lg(0)
		return livetree
		
	def get_extanttrees(self, compute=True, addextincts=[], collapsenodes=False, removelosses=True):
		"""returns a COPY of the list of all trees 'cleaned' of their dead branches ; NB: original tree will have all its nodes labelled afterward"""
		if not compute: return self._extanttrees
		# only works with all nodes being labelled, to use labels as references rather than the node objects (which refer to the original tree) in Node.pop() 
		self.labeltreenodes(silent=False) # should not be necessary though
		# NB: will exclude null trees (None objects) from the returned list
		# first remove lineages associated to extinct species
		lineages = [self.copy_prune_dead_lineages(t, trimroot=False, collapsenodes=collapsenodes) for t in self.genetrees if t]
		# and second remove lineages associated to gene loss
		if removelosses:
			if collapsenodes:
				extanttrees = [BaseTreeSimulator.copy_prune_dead_lineages(t, self.extincts, trimroot=False, collapsenodes=collapsenodes) for t in lineages if t]
			else:
				# keep a "rosary" structure of branches where single-child nodes are annotated as speciation-loss events
				extanttrees = [IOsimul.annotateSpeciationLossEvents(extanttree=t, lossnodes=self.extincts, trimLosses=True) for t in lineages if t]
		else:
			extanttrees = lineages
		return extanttrees
		
		
