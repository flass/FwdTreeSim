#!/usr/bin/python

import sys, os, copy
#~ import argparse
import getopt
from FwdTreeSim import models, simulators, IOsimul
import tree2
import random

def usage():
	crindent = "\t\t\t\t\t\t"
	l =  ["Usage:"]
	l += ["python %s [options]"%(sys.argv[0])]
	l += ["Options:"]
	l += ["  General options:"]
	l += ["\t-o  --outputdir path\t\tdirectory for all output files #\tdefaults to current directory."]
	
	l += ["  Simulation parameters:"]
	l += ["  __Species/Genomes population layer__:"]
	l += ["\t-s  --popsize\t\tint\t\tnumber of species to simulate in the underlying Moran process\t# default: 100."]
	l += ["\t-g  --ngen   \t\tint\t\tnumber of generations for which the evolution is simulated\t# default: 1000."]
	l += ["\t-c  --connect.lineages\tfloat\t\tthe multiple species lineage trees from the Moran process' population will", \
	      crindent+"all be connected at their root. The length of branches of the star-like", \
	      crindent+"root is given by the argument, negative value turns it off.", \
	      crindent+"# default: 0 (on)."]
	
	l += ["  __Gene/Locus layer__:"]
	l += ["\t-p  --profiles  path\t\t\tJSON file containing the (multiple) evolutionary profiles for simulated gene families, ", \
	      crindent+"and their respective weights (sampling probability if <=1 or expected number of gene families if >1).", \
	      crindent+"See FwdTreeSim/example/DTLprofiles.json for file format example.", \
	      crindent+"Example pangenome structure (with 20% core, 20% accessory, 60% orfan gene families)", \
	      crindent+"can be generated by replacing path by an empty string ''. "]
	l += ["\t-n  --ngenes    int\t\t\tnumber of gene families to simulate in the pangenome", \
	      crindent+"# default: 10; overriden by providing gene family profiles."]
	l += ["\t-r  --dtlrates  float[,float[,float]]\tglobal rates of Duplication, Transfer and Loss for ALL simulated gene families.", \
	      crindent+"# default: D=0.001 T=D, L=D+T; overriden by providing gene family profiles."]
	l += ["\t-f  --rootfreq  float\t\t\tglobal presence probability for ALL the gene families", \
	      crindent+"at the root of the simulation of the species population.", \
	      crindent+"# default: 0.5 ; overriden by providing gene family profiles."]
	
	l += ["  Output options:"]
	l += ["\t-e  --sample.extant.species\tint\thow many genomes are sampled in the end? trees are pruned accordingly", \
	      crindent+"# default: all sampled."]
	l += ["\t-l  --sample.larger.trees\tint\thow many lineage gene trees from a gene family population should be written out?", \
	      crindent+"Can be handy to just carry a diagnostic the simulated trees.", \
	      crindent+"# default: all lineage trees written as a single connected tree."]
	return '\n'.join(l)

def main():
	
	# option parsing
	try:
		opts, args = getopt.getopt(sys.argv[1:], "o:n:r:f:e:l:c:p:vh", ["outputdir=", "ngenes=", "popsize=", "ngen=", "connect.lineages=", \
																	"dtlrates=", "rootfreq=", "profiles=", \
																	"sample.larger.trees=", "sample.extant.species=", \
																	"help", "verbose"]) #, "connect.all.trees="
	except getopt.GetoptError as err:
		# print help information and exit:
		print str(err)  # will print something like "option -a not recognized"
		print usage()
		sys.exit(2)
		
	dopt = dict(opts)
	if ('-h' in dopt) or ('--help' in dopt):
		print usage()
		sys.exit(0)
	if ('-v' in dopt) or ('--verbose' in dopt): silent = False
	else: silent = True
	outdir = dopt.get('-o', dopt.get('--outputdir', os.getcwd()))
	popsize = int(dopt.get('-s', dopt.get('--popsize', 100)))
	ngen = int(dopt.get('-g', dopt.get('--ngen', 1000)))

	nfprofiles = dopt.get('-p', dopt.get('--profiles'))
	# define the evolution rates and original frequencies of gene families
	if nfprofiles!=None:
		if nfprofiles=='':
			exsampleprof = [(2, 'core'), (2, 'accessory-slow'), (6, 'orfan-fast')]
			dtlprof = IOsimul.MetaSimulProfile(profiles=[(n, IOsimul.DTLSimulProfile(type=t)) for n,t in exsampleprof])
		else:
			# expects a JSON formating of simulator profiles conforming to the IOsimul.MetaSimulProfile class parsers
			# e.g. a list of dict objects representing the arguments of a simulators.DTLSimulProfile class instance
			dtlprof = IOsimul.MetaSimulProfile(json=nfprofiles)
	else:
		# default global parameters (no pangenome structure), overridden by any provided profile
		rootfreq = dopt.get('-f', dopt.get('--rootfreq', 0.5))
		sdtlrates = dopt.get('-r', dopt['--dtlrates'])
		if sdtlrates:
			dtlrates = [float(s) for s in dopt.get('-r', dopt['--dtlrates']).split(',')]
		else:                                           # by default, rates are:
			dtlrates = [0.001] 							# D = 1e-3
		if len(dtlrates)<2: dtlrates += dtlrates[0]		# T = D
		if len(dtlrates)<3: dtlrates += [sum(dtlrates)] # L = T+D
		dglobalprof = {0:{'rdup':dtlrates[0], 'rtrans':dtlrates[1], 'rloss':dtlrates[2]}}
		globalprof = IOsimul.DTLSimulProfile(rateschedule=dglobalprof, rootfreq=rootfreq)
		dtlprof = IOsimul.MetaSimulProfile(profiles=[(1, globalprof)])
	
	# derive number of gene families to simulate from profile weights, or from dedicated option -n (overrides profiles), or take default value of 10
	if dtlprof.ngenes>1:
		ngenes = dtlprof.ngenes
	else:
		ngenes = int(dopt.get('-n', dopt.get('--ngenes', 10)))

	nlargegenetrees = int(dopt.get('-l', dopt.get('--sample.larger.trees', -1)))
	lentoroot = float(dopt.get('-c', dopt.get('--connect.lineages', 0)))
	samplextant = int(dopt.get('-e', dopt.get('--sample.extant.species', 0)))
	assert samplextant <= popsize

	#~ parser = argparse.ArgumentParser(description='Simulate phylogenic trees describing evolution of a population of bacterial genomes, with species, replicon/locus and gene layers.')
	#~ parser.add_argument('-o', '--outdir', )

	# creating output directories
	for d in ['logs', 'pickles', 'genetrees', 'reftrees']:
		outd = "%s/%s"%(outdir, d)
		if not os.path.exists(outd):
			os.mkdir(outd)

	# simualte species tree
	moranmodel = models.MoranProcess(popsize=popsize)
	moransim = simulators.MultipleTreeSimulator(model=moranmodel, ngen=ngen)
	if lentoroot>=0:
		# connect all roots of the species lineage trees
		conrt = moransim.connecttrees(lentoroot, returnCopy=True)
		conrt.write_newick("%s/reftrees/connected.reftree_full.nwk"%(outdir))
		# prune dead lineages and connect all roots of the species lineage trees
		extconrt = moransim.get_extanttree(compute=True, lentoroot=lentoroot)
		extconrt.write_newick("%s/reftrees/connected.reftree_extant.nwk"%(outdir))
		extantspe = extconrt.get_leaf_labels()
	else:
		# write lineage trees separately
		lextrt = moransim.get_extanttrees(compute=True)
		extantspe = []
		for k, extrt in enumerate(lextrt):
			extconrt.write_newick("%s/reftrees/reftree.%d_extant.nwk"%(outdir, k))
			extantspe += extconrt.get_leaf_labels()
			
	# select sampled species among the N extant
	if samplextant:
		sampledspe = random.sample(extantspe, samplextant)
		refnodeswithdescent = moransim.get_nodes_with_descendants(sample=sampledspe)
	else:
		refnodeswithdescent = moransim.get_nodes_with_descendants()
	# serial simulation of gene families, have to offer a parrallel version
	for k in range(ngenes):
		print "### simulate gene tree", k
		# simulate gene tree under the same reference tree set (= species/organism population history)
		bddtlsim = simulators.DTLtreeSimulator(refsimul=moransim, refnodeswithdescent=refnodeswithdescent, profile=dtlprof.sampleprofile(verbose=True), noTrigger=True)
		bddtlsim.evolve(bddtlsim.ngen)

		# connect all the gene trees in each gene population
		congt = bddtlsim.finish(connecttrees=True)
		
		# save ref and gene tree simulation object together to save space as they share references to same objects
		IOsimul.dumppickle({'refsim':moransim, 'genesim':bddtlsim}, "%s/pickles/simul.%d.pickle"%(outdir, k))
		bddtlsim.shadetreenodes()

		# write out the largest n gene trees and corresponding species trees
		if nlargegenetrees>=0:
			genetreesizes = [(genetree.nb_leaves(), i) for i, genetree in enumerate(bddtlsim.genetrees)]
			genetreesizes.sort(reverse=True)
			isavetrees = (genetreesizes[l][1] for l in range(nlargegenetrees))
		else:
			isavetrees = xrange(len(bddtlsim.genetrees))
		for l in isavetrees:
			genetree = bddtlsim.genetrees[l]
			gtoutrad = "%s/genetrees/simul.%d.all_gt"%(outdir, k)
			genetree.write_newick(gtoutrad+".nwk", mode=('w' if l==0 else 'a'))
			genetree.write_nexus(gtoutrad+".nex", mode=('w' if l==0 else 'a'))
			#~ genetree.ref.write_newick("%s/reftrees/simul.%d.rt.%d.nwk"%(outdir, k, l))
		
		# write out connected trees
		congtoutrad = "%s/genetrees/simul.%d.connected_gt_full"%(outdir, k)
		congt.write_newick(congtoutrad+".nwk")
		congt.write_nexus(congtoutrad+".nex")
		# prune dead lineages
		extcongtoutrad = "%s/genetrees/simul.%d.connected_gt_extant"%(outdir, k)
		extconrt = bddtlsim.get_extanttree(compute=True, lentoroot=lentoroot)
		if extconrt:
			extconrt.write_newick(extcongtoutrad+".nwk")
			extconrt.write_nexus(extcongtoutrad+".nex")
		else:
			with open(extcongtoutrad, "w") as fextcongtout:
				fextcongtout.write("no extant gene tree lineages.\n")
				
if __name__=='__main__':
	main()
