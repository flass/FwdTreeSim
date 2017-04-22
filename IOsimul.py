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

			

class SimulProfile(object):
	"""define generic class of profiles for simulator instances"""
	def __init__(self, **kwargs):
		print 'invoke simulator.SimulProfile.__init__()'
		self.rateschedule = kwargs.get('rateschedule', {}) 	# expects a dictionary with:
														#  keys = FIRST time slice of a range (the rates are checked for update at every time slice)
														#  values = a dictionary with: keys = rate names as in target model; values = the rate value for the range.
														# e.g. {0:{'rdup':.1, 'rtrans':.01, 'rloss':.02}, 100:{'rdup':.01, 'rtrans':.01, 'rloss':.02}}
														# starts at t=0 with relatively high duplication rate, leading to family expansion,
														# then from t=100 settles for evolution with an (expected) constant genome size

class DTLSimulProfile(SimulProfile):
	"""define class of profiles for DTL simulator instances; specifically provides the expected root frequency of the gene/element
	
	use shorthands such as 'core', 'accessory-slow', or 'orfan-fast' to access pre-defined profiles. Profiles can be created and stored in the class attribute
	"""
	dtypes = { \
	           'core' : { 'rootfreq':1, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'accessory-slow' : { 'rootfreq':0.5, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'accessory-fast' : { 'rootfreq':0.5, 'rateschedule': {0:{'rdup':.001, 'rtrans':.001, 'rloss':.002}} }, \
	           'orfan-slow' : { 'rootfreq':0.01, 'rateschedule': {0:{'rdup':.0001, 'rtrans':.0001, 'rloss':.0002}} }, \
	           'orfan-fast' : { 'rootfreq':0.01, 'rateschedule': {0:{'rdup':.001, 'rtrans':.001, 'rloss':.002}} } \
	          }
	
	def __init__(self, **kwargs):
		print 'invoke simulator.DTLSimulProfile.__init__()'
		print 'kwargs', kwargs
		super(DTLSimulProfile, self).__init__(**kwargs) 
		self.type = kwargs.get('type')
		self.rootfreq = kwargs.get('rootfreq')
		if self.type:
			if not self.rateschedule: self.rateschedule = DTLSimulProfile.dtypes[self.type]['rateschedule']
			if not self.rootfreq: self.rootfreq = DTLSimulProfile.dtypes[self.type]['rootfreq']

	
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
		"""reads profile from JSON file or string; see example for format"""
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
			prof = {'rootfreq':dprof['rootfreq'], 'rateschedule':{}}
			assert len(dprof['times'])==len(dprof['rates'])
			for i in range(len(dprof['times'])):
				if dprof['times'][i] in prof['rateschedule']:
					# verify that there is not redundant times listed for scheduled rate change
					raise ValueError, "input profile include several instance of the same time slice (t=%d) marked for changing evolutionary rates"%dprof['times'][i]
				else:
					prof['rateschedule'][dprof['times'][i]] = dprof['rates'][i]	
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
		if simultype:
			if 'DTLtreeSimulator' in simultype:
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

	def DTLeventLog(self, evt):
		self.foutdict['event_record'].write('\t'.join([models.DTLevent.etshorts[evt.eventtype], evt.donrefnode.nodeid()])+'\n')
		
		
	
