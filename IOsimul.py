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
from FwdTreeSim import deleteGenneratorAttr, _byteify
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
		
# for dumppickle() funciton; special case of simulators.BaseTreeSimulator class
simuldumpprompt =  "cannot pickle generator objects, have to delete event id generator 'self.eventidgen' first\n"
simuldumpprompt += "(will discontinue numeration of events for further simulation).\n"
simuldumpprompt += "Delete generator before pickling? (y/n) "
checkdic = {'simulators.BaseTreeSimulator':{'attrname':'eventidgen', 'prompt':simuldumpprompt}}
		
def dumppickle(obj, fileorpath, autoclosefile=True, prompt=False):
	def checkObj(obj):
		for cls in checkdic:
			if isinstance(obj, eval(cls)):
				checkDeleteGenneratorAttr(obj, attrname=checkdic[cls]['attrname'], prompt=(checkdic[cls]['prompt'] if prompt else False))
		
	if type(fileorpath)==str and os.path.exists(os.path.dirname(fileorpath)):
		fpickle = open(fileorpath, 'w')
		fpathstr = " in file '%s'"%fileorpath
	elif isinstance(fileorpath, file):
		if not os.access(os.W_OK):
			raise ValueError, "argument file object is not writeable"
		else:
			fpickle = fileorpath
		fpathstr = ""
	else:
		raise ValueError, "argument should be a file path or file-like object; '%s' is not"%(repr(fileorpath))
	
	# check that objects don't have annoying abjects (e.g. generators) to delete bfore pickling
	if isinstance(obj, list):
		for subobj in obj:
			checkObj(obj)
	elif isinstance(obj, dict):
		for subobj in obj.values():
			checkObj(obj)
	else:
		checkObj(obj)
			
	pickle.dump(obj, file=fpickle, protocol=2)
	print "saved %s in binary format"%repr(obj) + fpathstr
	if autoclosefile: fpickle.close()

			

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
		clsname = d['simulclass']
		lprofiles = d['profiles']
		lngprof = []
		if clsname.startswith('simulator.'): simprofcls = eval(clsname)
		else: cls = eval('simulator.'+clsname)
		for profile in lprofiles:
			ngenes = freqstrprof['ngenefams']
			dprof = freqstrprof['profile']
			prof = {'rootfreq':dprof['rootfreq'], 'rateschedule':{}}
			assert len(dprof['times'])==len(dprof['rates'])
			for i in range(len(dprof['times'])):
				if dprof['times'][i] in prof['rateschedule']:
					# verify that there is not redundant times listed for scheduled rate change
					raise ValueError, "input profile include several instance of the same time slice (t=%d) marked for changing evolutionary rates"%dprof['times'][i]
				else:
					prof['rateschedule'][dprof['times'][i]] = dprof['rates'][i]	
			lngprof.append(float(ngenes), simprofcls(prof))
		return lngprof
			
	def sampleprofile(self, verbose=False):
		p = random.random()
		for f in self.lfreq:
			if p<f:
				if verbose:
					print "use simulation profile: %s"%(repr(self.dprof[f]))
				return self.dprof[f]
			
class SimLogger(object):
	
	def __init__(self, sim, fout):
		self.sim = sim		# attached simulator object from which to get the output stream 
		self.out = fout		# file-like object to write the log stream
		
		
