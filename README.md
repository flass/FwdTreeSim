# FwdTreeSim

Forward phylogenetic simulation of bacterial pangenome, under a model with gene Duplication, horizontal Transfer and Loss (DTL). 

Description
-----------

The present DTL evolution model is derived from the xODT model described in [Szollosi, G. et al. 2013, Syst. Biol. 62(3):386-397](http://sysbio.oxfordjournals.org/content/62/3/386), for which a probabilistic inference tool, ALE, is available [here](https://github.com/ssolo/ALE).
In addition, diversity of evolutionary dynamics of gene families within a species' (or any coherent organism group's) pangenome is modelled, with classes of genes such as core, accessory and ORFan genes. Also the evolution of genome structure in terms of replicon units (chromosomes, plasmids) and gene linkage are modelled, in order to track the co-evolution of linked genes.

In this model, a species population evolves following a Moran process (Hey, J. 1992. Evolution, 46(3):627-640), i.e. with a fixed sized population of lineages undergoing one lineage extinction and one complementary speciation per time unit. This generates the species history, which is the bottom layer of the model. Importantly, the co-existence of lineages that went extinct at some point with lineages that have extant descendents is modelled.

Genomes evolve following the species history, but their component genes are subject to additonal evolutionary processes, namely gene duplication and loss within, and horizontal transfer between species lineages. This includes the possibility of horizontal transfer from/to lineages that ultimately will go extinct, thus accounting for the (overwhelming) contribution of extinct or unsamped lineages to the observed genetic diversity.
The second layer thus consists of a collection of gene families evolving in reference to the species history, and independently from each other.

This can be complexified by considering the linkage of genes within genome to incorporate the non-independence of evolution of linked genes. This is modelled through the addition of an intermediate layer tracking the evolution of replication units (replicons), and the possibility that DTL events span several neighbouring genes. While this is a more realistic rendering of how DNA macromolecule evolve, it also allows the user to test hypotheses on the co-selection of genesunder a neutral model where gene co-evolve.

Implementation
--------------

The nested nature of the model with genes evolving (within replicons) within genomes provides the opportunity to simulate their history in a incremental manner:

First, the species tree is simulated (e.g. under a Moran process, but other models are implemented). In the case of the Moran process, *n* disjoint species lineages (trees) evolve in parallel; this collection of lineages can be later connected at their root to provide a full species tree (with extinct lineages).

Second, each gene tree is simulated by copying the species tree and then proceeding from the initial time slice (where lineages have their root) to the final time slice (where are leaves of surviving lineages) and gradually editing this tree. Duplication, transfer and loss events occur stochastically in each gene tree branch crossing the current time slice.

Loss events are realized by delting the the subtree under the point of the event on the branch. Duplications are realized by copying the subtree under the event point and grafting that copy as a sister lineage of its template, hence creating a new node at the even point. Transfers are realized similarly, by copying the subtree under the point on the *receiving* branch and grafting that copy to the *donor* branch at the same time point, hence creating a new node there, under which descendents of the donor and recipients of the transfer form sister clades.

In prokaryotic genomes, gene families have heterogeneous frequencies of gene presence amongst species lineages. This is modelled by sampling a fraction of the lineage trees prior to the gene-level simulation, representing the lineages in which the gene is present at first; only these lineage trees will then be included in the gene-level simulation and be subject to DTL events.

By setting this frequency at the root and using DTL rates that have a null sum (D + T - L = 0, thus expecting constant genome size), it is straightforward to implement classes of genes that evolve while maintaining gene class-specific expected frequencies, e.g.:
  * core genes with freq fluctuating narrowly around 1 (f_root=1, D=.0001, T=.0001, L=.0002),
  * accessory genes with freq fluctuating within 0-1 interval moderately (f_root=.5, D=.0001, T=.0001, L=.0002) or more widely (f_root=.5, D=.001, T=.001, L=.002),
  * ORFan genes typically restricted to one or a few genomes, but often poping in and out genomes (f_root=.01, D=0, T=.002, L=.002),
  * or (breaking the balance of gains and losses) an invasive transposon (f_root=.01, D=.001, T=.001, L=.001).

Requirements 
------------

This library is implemented in pure Python (compatible with Python versions 2.7.*) and depends on the [tree2 library](https://github.com/flass/tree2).

Installation
------------

First open a command-line terminal (here using a `bash` shell) and clone the git repository:
```bash
# replace "/path/to/repo" by the path where you decide to clone the respective repositories
cd /path/to/repo
git clone https://github.com/flass/tree2
git clone https://github.com/flass/FwdTreeSim
```
Then simply add the path to the libraries to your PYTHONPATH environment variable (ideally save this in your `~/.bashrc` or `~/.bash_profile`):
```bash
export PYTHONPATH=$PYTHONPATH:/path/to/repo
```

Usage
-----

The set-up and execution of the simulation is fully scriptable (Python environment), with the possibility of providing various levels of detail and complexity in the simulated pangenome history.

```python
#!/usr/bin/python
from FwdTreeSim import models, simulators, IOsimul

## simulate the species tree
# set the model
moranmodel = models.MoranProcess(popsize=100)
# proceed with simulation over 1,000 generations
moransim = simulators.MultipleTreeSimulator(moranmodel, ngen=1000)

## use global parameters for all gene families
ngenes = 10		# number of independent gene families to simulate
				# global rates applied to all gene families
rdup = 0.0001	# rate of gene duplication
rtrans = 0.0001	# rate of horizontal gene transfer
rloss = 0.0001	# rate of gene loss
rootfreq = 0.5	# frequency of the gene

for k in range(ngenes):
	## simulate one gene family tree
	# set the model
	bddtlmodel = models.BirthDeathDTLModel(rdup=rdup, rtrans=rtrans, rloss=rloss, rootfreq=rootfreq)
	# set the simulator engine, providing the reference species history from which many attributes are inherited
	bddtlsim = simulators.DTLtreeSimulator(model=bddtlmodel, refsimul=moransim, noTrigger=True)
	# proceed with simulation for 100 generations
	bddtlsim.evolve(ngen=100)
	# or by default with the number of generations inherited from the species simulation
	bddtlsim.evolve()	# this call would be automatic if we'd not used the 'noTrigger=True' option above

```

Profiles of gene evolution with rates and expected root frequency can be provided either as a dict object, or through a key name corresponding to a typical gene class, or through an input JSON file ([example](https://github.com/flass/FwdTreeSim/blob/master/examples/DTLprofiles.json)).
Profiles also allow to set schedules of time heterogeneity in the process, e.g. with early family expansion by duplication, or later burst of transfers.
```python
# set a single global set of DTL rate parameters to be applied from t=0 onwards:
dprof =  { "rootfreq":0.5, "rateschedule": {0:{"rdup":0.0001, "rtrans":0.0001, "rloss":0.0002}} }
# or set different DTL rate parameter sets across time:
dprof = { "rootfreq":0.5, "rateschedule": {0:{"rdup":0.0005, "rtrans":0.0001, "rloss":0.0002}, 100:{"rdup":0.0001, "rtrans":0.0001, "rloss":0.0002}, 900:{"rdup":0.0001, "rtrans":0.001, "rloss":0.0002}} }
# or equivalently:
dprof = { "rootfreq":0.5, "times":[0, 100, 900], "rates": [{"rdup":0.0005, "rtrans":0.0001, "rloss":0.0002}, {"rdup":0.0001, "rtrans":0.0001, "rloss":0.0002}, {"rdup":0.0001, "rtrans":0.001, "rloss":0.0002}] }
# one can also specify the model to be implemented
dprof['modeltype'] = 'BirthDeathDTLModel'
# but no need as this is the default for the 'IOsimul.DTLSimulProfile' class

# run the simulation
prof = IOsimul.DTLSimulProfile(**prof)
bddtlsim = simulators.DTLtreeSimulator(refsimul=moransim, profile=prof)	# note that we don't have to specify the model as it is encoded in the profile

# multiple profiles can be specified (with weights) so one is picked at random prior to simulation
# provided as (freq, profile) tuples:
lprof = [(100, { "rootfreq":1.0, "rateschedule": {0:{"rdup":.0001, "rtrans":.0001, "rloss":.0002}} }), \
         (200, { "rootfreq":0.5, "rateschedule": {0:{"rdup":.0001, "rtrans":.0001, "rloss":.0002}} })]
dtlprof = IOsimul.MetaSimulProfile(profiles=[(n, IOsimul.DTLSimulProfile(**dprof)) for n, dprof in lprof])
# or using keyword handles:
pangeneprof = [(200, "core"), (200, "accessory-slow"), (600, "orfan-fast")]
dtlprof = IOsimul.MetaSimulProfile(profiles=[(n, IOsimul.DTLSimulProfile(type=t)) for n,t in pangeneprof])
# or using an external file input (can be useful for batch applications) in JSON format:
dtlprof = IOsimul.MetaSimulProfile(json='examples/DTLprofiles.json')

# simulate with one of those profiles (here the expected frequencies are normalized to get probabilities)
bddtlsim = simulators.DTLtreeSimulator(refsimul=moransim, profile=dtlprof.sampleprofile(verbose=True))

# (IN DEV) one can use similar profiles to generate a pangenome (collection of many gene families) in one go
# (here the expected frequencies are summed to specify the number of gene to simulate)
bddtlsim = multigene_simulators.DTLgenomeSimulator(refsimul=moransim, profile=dtlprof.sampleprofile(verbose=True))
# altermnatively one can explicitely specify the number of gene to simulate
bddtlsim = multigene_simulators.DTLgenomeSimulator(refsimul=moransim, ngenes=2000, profile=dtlprof.sampleprofile(verbose=True))

```

Alternatively, a command-line interface with the simulator is provided byt the script [bacterialGenomeDTL.py](https://github.com/flass/FwdTreeSim/blob/master/scripts/bacterialGenomeDTL.py). 
For details of options, please type:
```sh
python bacterialGenomeDTL.py --help
```
\[snapshot of usage/help message [here](https://github.com/flass/FwdTreeSim/blob/master/scripts/bacterialGenomeDTL.py_help_message.txt)\]
