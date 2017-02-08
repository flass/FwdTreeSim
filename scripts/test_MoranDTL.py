#!/usr/bin/python

import sys, os, copy
from FwdTreeSim import models, simulators, IOsimul
#~ import tree2

outdir = sys.argv[1]
dtlrate = [float(s) for s in sys.argv[2].split(',')]
for i in range(2, 4):
	if len(dtlrate)<i: dtlrate += dtlrate[:1]
nsims = int(sys.argv[3])
nlargegenetrees = int(sys.argv[4])

for d in ['pickles', 'genetrees', 'reftrees']:
	outd = "%s/%s"%(outdir, d)
	if not os.path.exists(outd):
		os.mkdir(outd)

# simualte species tree
moranmodel = models.MoranProcess(popsize=100)
moransim = simulators.MultipleTreeSimulator(moranmodel, ngen=1000)

for k in range(nsims):
	print "### simulate gene tree", k
	# simulate gene tree under the same reference tree set (= species population history)
	bddtlmodel = models.BirthDeathDTLModel(rdup=dtlrate[0], rtrans=dtlrate[1], rloss=dtlrate[2])
	bddtlsim = simulators.DTLtreeSimulator(model=bddtlmodel, refsimul=moransim, noTrigger=True)
	bddtlsim.evolve(bddtlsim.ngen)

	# save ref and gene tree simulation object together to save space as they share references to same objects
	IOsimul.dumppickle({'refsim':moransim, 'genesim':bddtlsim}, "%s/pickles/simul.%d.pickle"%(outdir, k), prompt=True)

	# write out the largest n gene trees and corresponding species trees
	genetreesizes = [(genetree.nb_leaves(), i) for i, genetree in enumerate(bddtlsim.genetrees)]
	largesizes = copy.copy(genetreesizes)
	largesizes.sort(reverse=True)
	for l in range(nlargegenetrees):
		il = largesizes[l]
		maxgenetree = bddtlsim.genetrees[il[1]]
		maxgenetree.write_newick("%s/genetrees/simul.%d.gt.%d.nwk"%(outdir, k, l))
		maxgenetree.ref.write_newick("%s/reftrees/simul.%d.rt.%d.nwk"%(outdir, k, l))
