#!/usr/bin/python
# -*- coding: utf-8 -*-

__all__ = ["models", "simulators", "genome", "IOsimul"]

__author__ = "Florent Lassalle <f.lassalle@imperial.ac.uk>"
__date__ = "22 August 2016"
__credits__ = """Leonor Palmeira and Laurent Guéguen for initiating the tree2.Node module."""

import os
import inspect
nodelabelprefix = dict(livetip='S', deadtip='E', node='N')
			
def deleteAttr(obj, attrname):
	"""delete object attributes which are generators in order to pickle the object"""
	del obj.__dict__[attrname]


def checkDeleteAttr(obj, attrname=None, prompt=None):
	"""check presence of object attributes to delete,  e.g. generator objects in order to pickle the owner object"""
	# cannot pickle generator objects, have to delete it
	if attrname in obj.__dict__:
		if prompt:
			if prompt==True:
				doit = raw_input("Delete generator attribute '%s' (y/n) "%repr(attrname))
			else:
				doit = raw_input(prompt)
			while not (doit in ['y', 'n']):
				print "answer 'y' (for yes) or 'n' (for no)"
				doit = raw_input(prompt)
			if doit=='y':
				deleteGenneratorAttr(obj, attrname=attrname)
			else:
				return 1
		else:
			deleteAttr(obj, attrname=attrname)
			return 0
			
def checkDeleteGenneratorAttr(obj, attrname=None, prompt=None):
	# check that objects don't have annoying objects (e.g. generators) to delete before pickling
	if isinstance(obj, list):
		for subobj in obj:
			checkDeleteGenneratorAttr(subobj)
	elif isinstance(obj, dict):
		for subobj in obj.values():
			checkDeleteGenneratorAttr(subobj)
	elif isinstance(obj, object):
		for atname in obj.__dict__.keys():
			if inspect.isgenerator(getattr(obj, atname)):
				if checkDeleteAttr(obj, attrname=atname, prompt=prompt):
					# deletion of attribute was not performed
					return 1
		else:
			return 0
	else:
		return -1

# utilitary functions for JSON parsing
def _byteify(data, ignore_dicts = False):
	"""parse multiple types of input and return byte (normal) strings instead of unicode strings"""
	# if this is a unicode string, return its string representation
	if isinstance(data, unicode):
		return data.encode('utf-8')
	# if this is a list of values, return list of byteified values
	if isinstance(data, list):
		return [ _byteify(item, ignore_dicts=True) for item in data ]
	# if this is a dictionary, return dictionary of byteified keys and values
	# but only if we haven't already byteified it
	if isinstance(data, dict) and not ignore_dicts:
		return { _byteify(key, ignore_dicts=True): _byteify(value, ignore_dicts=True) for key, value in data.iteritems() }
	# if it's anything else, return it in its original form
	return data
