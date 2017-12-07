#!/usr/bin/python
# -*- coding: utf-8 -*-

"""Simulate phylogenetic tree(s) forward in time, notably under models recording events of duplication, transfer and loss for gene trees relative to a species tree."""

__author__ = "Florent Lassalle <florent.lassalle@imperial.ac.uk>"
__date__ = "27 July 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import os
import cPickle as pickle
#~ import pickle
import random
import json
from FwdTreeSim import checkDeleteGenneratorAttr, _byteify
#~ import FwdTreeSim.simulators as simulators
import tree2

## pickle serialized object dump/load
def loadpickle(fileorpath, autoclosefile=True):
	if type(fileorpath)==str and os.path.exists(os.path.dirname(fileorpath)):
		fpickle = open(fileorpath, 'r')
	elif isinstance(fileorpath, file):
		if not os.access(os.R_OK):
			raise ValueError, "argument file object is not readable"
		else:
			fpickle = fileorpath
	else:
		raise ValueError, "argument should be a file path or file-like object; '%s' is not"%(repr(fileorpath))
	obj = pickle.load(fpickle)
	if autoclosefile: fpickle.close()
	return obj

def dumppickle(obj, fileorpath, autoclosefile=True, prompt=False, silent=True):
	if type(fileorpath)==str and os.path.exists(os.path.dirname(fileorpath)):
		fpickle = open(fileorpath, 'w')
		fpathstr = fileorpath
	elif isinstance(fileorpath, file):
		if not os.access(os.W_OK):
			raise ValueError, "argument file object is not writeable"
		else:
			fpickle = fileorpath
		fpathstr = repr(fileorpath)
	else:
		raise ValueError, "argument should be a file path or file-like object; '%s' is not"%(repr(fileorpath))
	# cannot pickle generator objects, have to delete them first from the object attributes
	if not silent:
		dumpwarningmsg = getattr(obj, 'dumppickle_warning', None)
		if dumpwarningmsg: print dumpwarningmsg
	if checkDeleteGenneratorAttr(obj, prompt=prompt) < 1:
		# no generator to delete or deletion was performed
		pickle.dump(obj, file=fpickle, protocol=2)
		if not silent: print "saved %s in binary format in file '%s'"%(repr(obj), fpathstr)
		if autoclosefile: fpickle.close()
	else:
		print "could not save object %s due to failure to delete its attribute generator object(s)"%repr(obj)

def make_rateschedule(ltimes, lrates, drateschedule={}):
	assert len(times)==len(lrates)
	for t, i in enumerate(ltimes):
		if t in drateschedule:
			# verify that there is not redundant times listed for scheduled rate change
			raise ValueError, "input profile include several instance of the same time slice (t=%d) marked for changing evolutionary rates, or is redundant with a previously specified one"%t
		else:
			drateschedule[t] = lrates[i]				

class SimulProfile(object):
	"""define generic class of profiles for simulator instances
	
	by default the attached simulator model is 'MoranProcess'
	"""
	def __init__(self, **kwargs):
		print 'invoke simulator.SimulProfile.__init__()'
		self.modeltype = kwargs.get('modeltype', 'MoranProcess')
		self.rateschedule = kwargs.get('rateschedule', {}) 	# expects a dictionary with:
														#  keys = FIRST time slice of a range (the rates are checked for update at every time slice)
														#  values = a dictionary with: keys = rate names as in target model; values = the rate value for the range.
														# e.g. {0:{'rdup':.1, 'rtrans':.01, 'rloss':.02}, 100:{'rdup':.01, 'rtrans':.01, 'rloss':.02}}
														# starts at t=0 with relatively high duplication rate, leading to family expansion,
														# then from t=100 settles for evolution with an (expected) constant genome size
		self.read_compound_rateschedule_fmt(kwargs)
	
	def read_compound_rateschedule_fmt(self, dprof):
		ltimes = dprof.get('times', [])
		lrates = dprof.get('rates', [])
		make_rateschedule(ltimes, lrates, self.rateschedule)

class DTLSimulProfile(SimulProfile):
	"""define class of profiles for DTL simulator instances; specifically provides the expected root frequency of the gene/element
	
	by default the attached simulator model is 'BirthDeathDTLModel'
	Use shorthands such as 'core', 'accessory-slow', or 'orfan-fast' to access pre-defined profiles. 
	Profiles can be created and stored in the class attribute
	"""
	dtypes = { \
	           'core' : { 'rootfreq':1, 'rootlen':1, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'accessory-slow' : { 'rootfreq':0.5, 'rootlen':10, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'accessory-fast' : { 'rootfreq':0.5, 'rootlen':30, 'rateschedule': {0:{'rdup':.001, 'rtrans':.001, 'rloss':.002}} }, \
	           'orfan-slow' : { 'rootfreq':0.01, 'rootlen':10, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'orfan-fast' : { 'rootfreq':0.01, 'rootlen':30, 'rateschedule': {0:{'rdup':.001, 'rtrans':.001, 'rloss':.002}} } \
	          }
	
	def __init__(self, **kwargs):
		print 'invoke simulator.DTLSimulProfile.__init__()'
		print 'kwargs', kwargs
		super(DTLSimulProfile, self).__init__(**kwargs) 
		self.modeltype = kwargs.get('modeltype', 'BirthDeathDTLModel')
		self.type = kwargs.get('type', 'core')
		self.rootfreq = kwargs.get('rootfreq')
		self.rootlen = kwargs.get('rootlen')
		for attr in ['rateschedule', 'rootfreq', 'rootlen']:
			if not getattr(self, attr): setattr(self, attr, DTLSimulProfile.dtypes[self.type][attr])
		# particular metaparameter for tree length multiplier (emulates sequence diversity in gene family),
		# which value can be set, but otherwise is not assumed to be related to the gene type
		if not getattr(self, 'multreelen'): setattr(self, 'multreelen', ('gamma', 2, 0.5))

	
class MetaSimulProfile(object):
	"""Profiles for specifying evolutionary gene classes of genes and their respective frequencies in the genome. 
	
	This class is only used for I/O of profiles.
	"""
	
	def __init__(self, profiles=None, **kwargs):
		"""Profiles for specifying evolutionary gene classes of genes and their respective frequencies in the genome. Can be provied though:
		
		* 'json' keyword argument, containing a JSON-formated file specifying the profiles (see example);
		* 'profiles' argument, which must be a list of tuples, such [(A1,B1), (A2,B2), ...], with:
		  - Ax being a float/integer (expected frequency [i.e. number/fraction of gene families] 
		     at which a simulation profile is used (all frequencies are normalized to sum as 1, i.e. turned into probabilities);
		  - Bx being a SimulProfile object (specification of the gene profile).
		"""
		self.dprof = {}
		self.lfreq = []
		self.ngenes = 0
		
		if 'json' in kwargs:
			# load lfreq and dprof from JSON format file
			lngprof = self.loadJSON(kwargs['json'])
		elif profiles:
			lngprof = profiles
		sumfreq = sum([float(ng) for ng, prof in lngprof])
		cfreq = float(0)
		for i in range(len(lngprof)):
			ng, prof = lngprof[i]
			if ng<=0: raise ValueError, "profile frequencies (_expected_ number of gene families simulated under this model profile) must be strictly positive; got:\n%s"%repr(ng)
			if not isinstance(prof, SimulProfile): raise TypeError, "expected a IOSimul.SimulProfile class (or derivative) instance; got:\n%s"%repr(prof)
			# record cummulative frequencies over [0,1] so to define probability thresholds matching keys in the 'dprof' dict 
			if i<(len(lngprof)-1):
				f = cfreq + float(ng)/sumfreq
			else:
				# make sure upper bound is 1 in spite of sum of approximated fractions
				f = 1
			self.dprof[f] = prof
			self.lfreq.append(f)
			self.ngenes += ng
			cfreq = f 
		# enforce an integer value for number of gene families to simulate
		self.ngenes = max(int(self.ngenes), 1)
			
	def loadJSON(self, jsonfile):
		"""reads profile from JSON file or string; see example for format
		
		only the final tree length multiplier ('multreelen') is defined by default 
		(will draw from a Gamma ditribution with parameters (k=2, thetha=0.5)) 
		and can be omitted from the JSON file.
		"""
		with open(jsonfile, 'r') as fprofiles:
			d = json.load(fprofiles, object_hook=_byteify)
		clsname = d['simprofclass']
		lprofiles = d['profiles']
		lngprof = []
		#~ if clsname.startswith('IOsimul.'): simprofcls = eval(clsname)
		#~ else: simprofcls = eval('IOsimul.'+clsname)
		simprofcls = eval(clsname)
		for profile in lprofiles:
			ngenes = profile['ngenefams']
			dprof = profile['profile']
			prof = {'rootfreq':dprof['rootfreq'], 'rootlen':dprof['rootlen'], 'multreelen':dprof.get('multreelen', ('gamma', 2, 0.5)), 'rateschedule':dprof.get('rateschedule', {})}
			if not prof['rateschedule']:
				make_rateschedule(dprof.get('times', []), dprof.get('rates', []), prof['rateschedule'])
			lngprof.append((float(ngenes), simprofcls(**prof)))
		return lngprof
			
	def sampleprofile(self, verbose=False):
		p = random.random()
		for f in self.lfreq:
			if p<f:
				if verbose:
					print "use simulation profile: %s"%(repr(self.dprof[f]))
				return self.dprof[f]
			
class SimulLogger(object):
	"""Generate dump files for building a database of  simulation models and scenarios and resulting trees and events."""
	
	def __init__(self, table_fields=None, simultype=None, bnfout="log_", mode='w'):
		"""(create and) open connections to the output database dump file.
		
		'table2fields' argument is a dict with the desired database's table names as keys 
		and a list of table fields as values.
		"""
		table2fields = {}
		if simultype:
			if str(simultype) == "<class 'FwdTreeSim.simulators.DTLtreeSimulator'>":
				table2fields = {'species_tree_record':['species_branch_name', 'Duplications', 'Transfers', 'Losses', 'Originations', 'copies'], \
								'undated_transfer_record':['from', 'to', 'freq'], \
								'event_record':['evt_type', 'from_species_rank', 't_out', 'from_species_branch_name', 'from_gene_node_id', 'to_species_rank', 't_back', 'to_species_branch_name', 'to_gene_node_id'] \
								}
		else:
			table2fields = table_fields
			
		# collection of file-like objects to write the different log streams
		self.foutdict = {}
		for tablename in table2fields:
			tabledump = open(bnfout+"%s.tsv"%tablename, mode)
			# writes table header
			tabledump.write('\t'.join(table2fields[tablename])+'\n')
			# store file handle
			self.foutdict[tablename] = tabledump
		
	def close(self):
		"""close all file connections"""
		for tabledump in self.foutdict.values():
			tabledump.close()

	def DTLsingleEventLog(self, evt):
		self.foutdict['event_record'].write('\t'.join([str( evt.etshorts[evt.eventtype] ) , str(evt.donrefnode.nodeid()), str( evt.t ) , ])+'\n')
		
	def DTLsummaryEventLog(self, evt):
		self.foutdict['undated_transfer_record']
	
		
def annotateSpeciationLossEvents(trimLosses=False, **kw):
	"""takes an input tree and returns it annotated tree with SL events at speciation nodes which one child lineage ends with a loss event
	
	The input tree can be passed trhough kw argument 'extanttree' or this may be extracted from a simulation passed as kw argument 'simul'.
	optionally prunes off the loss-ending lineage (preserving the node when not stading alone).
	"""	
	def getheadnodes(node, removednodelabels=None):
		"""return the head node and the one just under it on a "rosary"-like chain of single-child nodes"""
		rosnode = node
		frosnode = rosnode.go_father()
		while frosnode and frosnode.nb_children()==1:
			if not (removednodelabels is None): removednodelabels.append(rosnode.label())
			rosnode = frosnode
			frosnode = rosnode.go_father()
		if frosnode:
			headnode = frosnode
			subheadnode = rosnode
		else:
			# rosnode has no father <=> is root
			headnode = rosnode
			subheadnode = rosnode.get_children()[0]
		return (headnode, subheadnode)
		
	def trimlosses(headnode, subheadnode):
		headnode.unlink_child(subheadnode)
		while headnode.nb_children()==0:
			# reach the next headnode with non-null descendance
			headnode, subheadnode = getheadnodes(headnode, removednodelabels=removednodelabels)
			headnode.unlink_child(subheadnode)
		return (headnode, subheadnode)
		
	simul = kw.get('simul')
	extanttree = kw.get('extanttree', simul.extanttree)
	if simul: llossnodes = kw.get('lossnodes', simul.extincts)
	else: llossnodes = kw.get('lossnodes', extanttree.get_postordertraversal_children())
	if trimLosses: removednodelabels = []
	for node in llossnodes:
		if node.label() in removednodelabels: continue
		# verify node is a loss event node
		nodev = getattr(node, 'event', None)
		if nodev and nodev[1]=='loss':
			headnode, subheadnode = getheadnodes(node, removednodelabels=removednodelabels)
			# tag the the lineage-head node with the SL event (same event id than the lineage-end loss event)
			headnode.event = (nodev[0], 'speciationloss')
			if trimLosses:
				headnode, subheadnode = trimlosses(headnode, subheadnode)
				if headnode.is_root():
					# all lineages were ending by losses; tree is now void
					return None
	return extanttree

# could be in simulators
def traceback_DTLevent_chain(event, simul, **kw):
	"""return the chain of DTL event ids that connects an event's node with the next extant lineage, that is the chain of event that occured over a gene tree branch
	
	takes as argument a models.DTLevent instance and a simulators.DTLtreeSimulator instance by which it was generated, and from which the following attributes will be fetched: 
	a gene tree pruned of its extinct lineages ('extanttree'), a map of gene tree node labels to event ids ('eventsmap') and the list of events affecting extant lineages ('extantevents').
	"""
	endrecnode = simul.extanttree[event.recipient().label()]
	eventchain = []
	recnode = endrecnode
	frecnode = recnode.go_father()
	while frecnode and frecnode.nb_children()==1:
		upevtid = simul.eventsmap.get(recnode.label())
		if upevtid and (upevtid in simul.extantevents):
			eventchain.append(upevtid)
		recnode = frecnode
		frecnode = recnode.go_father()
	return eventchain
				
			
	
