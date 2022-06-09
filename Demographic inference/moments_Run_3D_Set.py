#!/usr/bin/python

#import os
import numpy
import scipy
import matplotlib
import matplotlib.pyplot as pyplot
import cython
import sys
#sys.path.append("/Users/luqman/Documents/Hirzi/ETHPHD/Moments")
import moments
import pylab
import time
from datetime import datetime

## Set working directory
#os.chdir("/Users/luqman/Desktop/Demographic_Inference_Moments/Moments analysis")
import Optimize_Functions
#import Models_3D
#import Models_3D_DurCivSimTopo
#import Models_3D_CivDurSimTopo
import Models_3D_BalkansApenninesCentralAlpsTopo
#import Models_3D_ApenninesBalkansCentralAlpsTopo

 # To parallelise across CPUs, we allow parallelisation across models
model_num = sys.argv[1]
#model_num = 1
model_num = int(model_num)

'''
Usage: python moments_Run_3D_Set.py 1

This is a modified version of the 'moments_Run_Optimizations.py' script in which
we run optimizations for 3D comparisons for a set of models. These models are stored in the
Models_3D.py script, and will be called directly here. The user can delete or
comment out models to analyze a subset of the models available. 

This script must be in the same working directory as Optimize_Functions.py, which
contains all the functions necessary, as well as the  Models_3D.py script, which
has all the model definitions.


General workflow:
 The optimization routine runs a user-defined number of rounds, each with a user-defined
 or predefined number of replicates. It's also good to run the optimization routine 
 multiple times (chains). The starting parameters are initially random, but after
 each round is complete the parameters of the best scoring replicate from that round are
 used to generate perturbed starting parameters for the replicates of the subsequent round.
 The arguments controlling steps of the optimization algorithm (maxiter) and perturbation
 of starting parameters (fold) can be supplied by the user for more control across rounds.
 The user can also supply their own set of initial parameters, or set custom bounds on the
 parameters (upper_bound and lower_bound) to meet specific model needs. This flexibility
 should allow these scripts to be generally useful for model-fitting with any data set.

 
Outputs:
 For each model run, there will be a log file showing the optimization steps per replicate
 and a summary file that has all the important information. Here is an example of the output
 from a summary file, which will be in tab-delimited format:
 
 Model	Replicate	log-likelihood	AIC	chi-squared	theta	optimized_params(nu1, nu2, m, T)
 sym_mig	Round_1_Replicate_1	-1684.99	3377.98	14628.4	383.04	0.2356,0.5311,0.8302,0.182
 sym_mig	Round_1_Replicate_2	-2255.47	4518.94	68948.93	478.71	0.3972,0.2322,2.6093,0.611
 sym_mig	Round_1_Replicate_3	-2837.96	5683.92	231032.51	718.25	0.1078,0.3932,4.2544,2.9936
 sym_mig	Round_1_Replicate_4	-4262.29	8532.58	8907386.55	288.05	0.3689,0.8892,3.0951,2.8496
 sym_mig	Round_1_Replicate_5	-4474.86	8957.72	13029301.84	188.94	2.9248,1.9986,0.2484,0.3688


Notes/Caveats:
 The likelihood and AIC returned represent the true likelihood only if the SNPs are
 unlinked across loci. If SNPs are linked across loci then the likelihood is actually a composite
 likelihood and using something like AIC is no longer appropriate for model comparisons.
 See your notes for more information on this subject. 

Citations:
 If you use these scripts or the sets of diversification models for your work, please
 cite the following publications:
    Portik, D.M., Leache, A.D., Rivera, D., Blackburn, D.C., Rodel, M.-O.,
    Barej, M.F., Hirschfeld, M., Burger, M., and M.K. Fujita. 2017.
    Evaluating mechanisms of diversification in a Guineo-Congolian forest
    frog using demographic model selection. Molecular Ecology 26: 5245-5263.
    doi: 10.1111/mec.14266
    
    Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, H.C., Nagel, P., Onstein, R.E., 
	Portik, D.M., Streicher, J.W., and S.P. Loader. Vanishing refuge: testing the 
	forest refuge hypothesis in coastal East Africa using genome-wide sequence data 
	for five co-distributed amphibians. In Review, Molecular Ecology.

-------------------------
Written for Python 2.7
Python modules required:
-Numpy
-Scipy
-moments
-------------------------
(Original script by Daniel Portik)
daniel.portik@gmail.com
https://github.com/dportik
Updated May 2018

(Modified to work in Moments and added parallelising capability by Hirzi Luqman)
Feb 2019

'''

#===========================================================================
# Import data to create joint-site frequency spectrum
#===========================================================================

### Create AFS from file (output of realSFS with header)
#data = moments.Spectrum.from_file("3DSFS_DurCivSim_polarised.sfs") 
#data = moments.Spectrum.from_file("3DSFS_KonDolSin_polarised.sfs")
#data = moments.Spectrum.from_file("Kapetanovo_Jezero.Veaux.Bayasse.wHeader.sfs")
data = moments.Spectrum.from_file("Kapetanovo_Jezero.Civita.Polsa.wHeader.sfs")

### Fold the AFS
#fs = data.fold()

### Population IDs
#pop_ids = ["Durmitor", "Civita", "Simplonpass"]
#pop_ids = ["Konitsa", "Mt_Dolcedorme", "Sinlio"]
#pop_ids = ["Kapetanovo_Jezero", "Veaux", "Bayasse"]
pop_ids = ["Kapetanovo_Jezero", "Civita", "Polsa"]

###create a prefix based on the population names to label the output files
#ex. Pop1_Pop2_Pop3
prefix = "_".join(pop_ids)

###print some useful information about the afs or jsfs
print "\n\n============================================================================\nData for site frequency spectrum\n============================================================================\n"
#print "projection", proj
print "sample sizes", data.sample_sizes
sfs_sum = numpy.around(data.S(), 2)
print "Sum of SFS = ", sfs_sum, '\n', '\n'

#%%
#================================================================================
# Calling external 3D models from the Models_3D.py script
#================================================================================
'''
 We will use a function from the Optimize_Functions.py script for our optimization routines:
 
 Optimize_Routine(fs, outfile, model_name, func, rounds, param_number, fs_folded=False, reps=None, maxiters=None, folds=None, in_params=None, in_upper=None, in_lower=None, param_labels=" ")
 
   Mandatory Arguments =
    fs:  spectrum object name
    outfile:  prefix for output naming
    model_name: a label to help label the output files; ex. "no_mig"
    func: access the model function from within 'moments_Run_Optimizations.py' or from a separate python model script, ex. after importing Models_3D, calling Models_3D.no_mig
    rounds: number of optimization rounds to perform
    param_number: number of parameters in the model selected (can count in params line for the model)
    fs_folded: A Boolean value (True or False) indicating whether the empirical fs is folded (True) or not (False).

   Optional Arguments =
     chains: number of times to repeat the optimization routine
     reps: a list of integers controlling the number of replicates in each of the optimization rounds
     maxiters: a list of integers controlling the maxiter argument in each of the optimization rounds
     folds: a list of integers controlling the fold argument when perturbing input parameter values
     in_params: a list of parameter values 
     in_upper: a list of upper bound values
     in_lower: a list of lower bound values
     param_labels: list of labels for parameters that will be written to the output file to keep track of their order

'''

#**************
# Define the lists for optional arguments
# You can change these to alter the settings of the optimization routine

### TEST routine
"""
chains = 2
rounds = 2
reps = [3,3]
maxiters = [3,5]
folds = [3,2]
"""
### FINAL routine
chains = 5
rounds = 4
reps = [10,20,30,40]
maxiters = [3,5,10,15]
folds = [3,2,2,1]
"""
### FINAL routine (comprehensive; for more complex models)
chains = 6
rounds = 6
reps = [10,20,30,40,50,50]
maxiters = [3,5,10,15,20,20]
folds = [3,3,2,2,1,1]
"""
#**************
### Indicate whether your frequency spectrum object is folded (True) or unfolded (False)
fs_folded = False

"""
This script has been modified to allow parallelisation across jobs (CPUs). Simply run this script allocating each model to a different CPU (via the required argument)
"""
###### Dianthus sylvestris divergence set ######
# Note: make sure the model in the lists below are all in the same order!

### Model names
modelname_list = ["split_nomig", "split_nomig_growth", "split_nomig_growthall", "split_symmig_all", "split_symmig_all_growth", "split_symmig_all_growthall", "refugia_adj_1", "refugia_adj_1_growth", "refugia_adj_1_growthall", "refugia_adj_2", "refugia_adj_3", "refugia_adj_4", "refugia_adj_4_nomig", "refugia_adj_4_growth", "ancmig_adj_1", "ancmig_adj_2", "ancmig_adj_3", "sim_split_no_mig", "sim_split_no_mig_size", "sim_split_sym_mig_all", "sim_split_sym_mig_all_size", "refugia_adj_5", "refugia_adj_5_simsplit_isol", "refugia_adj_5_simsplit", "refugia_adj_5_simsplit_3epochs", "refugia_adj_5_full", "refugia_adj_5_simsplit_4epochs_iter1", "refugia_adj_5_simsplit_4epochs_iter2", "refugia_adj_5_simsplit_4epochs_iter3", "refugia_adj_5_simsplit_4epochs_iter4", "refugia_adj_5_simsplit_4epochs_iter5", "refugia_adj_5_full_2_iter1", "refugia_adj_5_full_2_iter2", "refugia_adj_5_full_2_iter3", "refugia_adj_5_full_2_iter4", "refugia_adj_5_full_2_iter5", "split_full_3epochs_iter1", "split_full_3epochs_iter2", "split_full_3epochs_iter3", "split_full_3epochs_iter4", "split_full_3epochs_iter5", "split_full_4epochs_iter1", "split_full_4epochs_iter2", "split_full_4epochs_iter3", "split_full_4epochs_iter4", "split_full_4epochs_iter5", "split_simsplit_3epochs_iter1", "split_simsplit_3epochs_iter2", "split_simsplit_3epochs_iter3", "split_simsplit_3epochs_iter4", "split_simsplit_3epochs_iter5"]

### Model functions
#ModelFunctions = [Models_3D.split_nomig, Models_3D.split_nomig_growth, Models_3D.split_nomig_growthall, Models_3D.split_symmig_all, Models_3D.split_symmig_all_growth, Models_3D.split_symmig_all_growthall, Models_3D.refugia_adj_1, Models_3D.refugia_adj_1_growth, Models_3D.refugia_adj_1_growthall, Models_3D.refugia_adj_2, Models_3D.refugia_adj_3, Models_3D.refugia_adj_4, Models_3D.refugia_adj_4_nomig, Models_3D.refugia_adj_4_growth, Models_3D.ancmig_adj_1, Models_3D.ancmig_adj_2, Models_3D.ancmig_adj_3, Models_3D.sim_split_no_mig, Models_3D.sim_split_no_mig_size, Models_3D.sim_split_sym_mig_all, Models_3D.sim_split_sym_mig_all_size, Models_3D.refugia_adj_5, Models_3D.refugia_adj_5_simsplit_isol, Models_3D.refugia_adj_5_simsplit, Models_3D.refugia_adj_5_simsplit_3epochs, Models_3D.refugia_adj_5_full, Models_3D.refugia_adj_5_simsplit_4epochs, Models_3D.refugia_adj_5_full_2]
#ModelFunctions = [Models_3D_DurCivSimTopo.split_nomig, Models_3D_DurCivSimTopo.split_nomig_growth, Models_3D_DurCivSimTopo.split_nomig_growthall, Models_3D_DurCivSimTopo.split_symmig_all, Models_3D_DurCivSimTopo.split_symmig_all_growth, Models_3D_DurCivSimTopo.split_symmig_all_growthall, Models_3D_DurCivSimTopo.refugia_adj_1, Models_3D_DurCivSimTopo.refugia_adj_1_growth, Models_3D_DurCivSimTopo.refugia_adj_1_growthall, Models_3D_DurCivSimTopo.refugia_adj_2, Models_3D_DurCivSimTopo.refugia_adj_3, Models_3D_DurCivSimTopo.refugia_adj_4, Models_3D_DurCivSimTopo.refugia_adj_4_nomig, Models_3D_DurCivSimTopo.refugia_adj_4_growth, Models_3D_DurCivSimTopo.ancmig_adj_1, Models_3D_DurCivSimTopo.ancmig_adj_2, Models_3D_DurCivSimTopo.ancmig_adj_3, Models_3D_DurCivSimTopo.sim_split_no_mig, Models_3D_DurCivSimTopo.sim_split_no_mig_size, Models_3D_DurCivSimTopo.sim_split_sym_mig_all, Models_3D_DurCivSimTopo.sim_split_sym_mig_all_size, Models_3D_DurCivSimTopo.refugia_adj_5, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_isol, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_3epochs, Models_3D_DurCivSimTopo.refugia_adj_5_full, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_4epochs_iter1, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_4epochs_iter2, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_4epochs_iter3, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_4epochs_iter4, Models_3D_DurCivSimTopo.refugia_adj_5_simsplit_4epochs_iter5, Models_3D_DurCivSimTopo.refugia_adj_5_full_2_iter1, Models_3D_DurCivSimTopo.refugia_adj_5_full_2_iter2, Models_3D_DurCivSimTopo.refugia_adj_5_full_2_iter3, Models_3D_DurCivSimTopo.refugia_adj_5_full_2_iter4, Models_3D_DurCivSimTopo.refugia_adj_5_full_2_iter5, Models_3D_DurCivSimTopo.split_full_3epochs_iter1, Models_3D_DurCivSimTopo.split_full_3epochs_iter2, Models_3D_DurCivSimTopo.split_full_3epochs_iter3, Models_3D_DurCivSimTopo.split_full_3epochs_iter4, Models_3D_DurCivSimTopo.split_full_3epochs_iter5, Models_3D_DurCivSimTopo.split_full_4epochs_iter1, Models_3D_DurCivSimTopo.split_full_4epochs_iter2, Models_3D_DurCivSimTopo.split_full_4epochs_iter3, Models_3D_DurCivSimTopo.split_full_4epochs_iter4, Models_3D_DurCivSimTopo.split_full_4epochs_iter5, Models_3D_DurCivSimTopo.split_simsplit_3epochs_iter1, Models_3D_DurCivSimTopo.split_simsplit_3epochs_iter2, Models_3D_DurCivSimTopo.split_simsplit_3epochs_iter3, Models_3D_DurCivSimTopo.split_simsplit_3epochs_iter4, Models_3D_DurCivSimTopo.split_simsplit_3epochs_iter5]
#ModelFunctions = [Models_3D_CivDurSimTopo.split_nomig, Models_3D_CivDurSimTopo.split_nomig_growth, Models_3D_CivDurSimTopo.split_nomig_growthall, Models_3D_CivDurSimTopo.split_symmig_all, Models_3D_CivDurSimTopo.split_symmig_all_growth, Models_3D_CivDurSimTopo.split_symmig_all_growthall, Models_3D_CivDurSimTopo.refugia_adj_1, Models_3D_CivDurSimTopo.refugia_adj_1_growth, Models_3D_CivDurSimTopo.refugia_adj_1_growthall, Models_3D_CivDurSimTopo.refugia_adj_2, Models_3D_CivDurSimTopo.refugia_adj_3, Models_3D_CivDurSimTopo.refugia_adj_4, Models_3D_CivDurSimTopo.refugia_adj_4_nomig, Models_3D_CivDurSimTopo.refugia_adj_4_growth, Models_3D_CivDurSimTopo.ancmig_adj_1, Models_3D_CivDurSimTopo.ancmig_adj_2, Models_3D_CivDurSimTopo.ancmig_adj_3, Models_3D_CivDurSimTopo.sim_split_no_mig, Models_3D_CivDurSimTopo.sim_split_no_mig_size, Models_3D_CivDurSimTopo.sim_split_sym_mig_all, Models_3D_CivDurSimTopo.sim_split_sym_mig_all_size, Models_3D_CivDurSimTopo.refugia_adj_5, Models_3D_CivDurSimTopo.refugia_adj_5_simsplit_isol, Models_3D_CivDurSimTopo.refugia_adj_5_simsplit, Models_3D_CivDurSimTopo.refugia_adj_5_simsplit_3epochs, Models_3D_CivDurSimTopo.refugia_adj_5_full, Models_3D_CivDurSimTopo.refugia_adj_5_simsplit_4epochs, Models_3D_CivDurSimTopo.refugia_adj_5_full_2]
ModelFunctions = [Models_3D_BalkansApenninesCentralAlpsTopo.split_nomig, Models_3D_BalkansApenninesCentralAlpsTopo.split_nomig_growth, Models_3D_BalkansApenninesCentralAlpsTopo.split_nomig_growthall, Models_3D_BalkansApenninesCentralAlpsTopo.split_symmig_all, Models_3D_BalkansApenninesCentralAlpsTopo.split_symmig_all_growth, Models_3D_BalkansApenninesCentralAlpsTopo.split_symmig_all_growthall, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_1, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_1_growth, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_1_growthall, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_2, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_3, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_4, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_4_nomig, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_4_growth, Models_3D_BalkansApenninesCentralAlpsTopo.ancmig_adj_1, Models_3D_BalkansApenninesCentralAlpsTopo.ancmig_adj_2, Models_3D_BalkansApenninesCentralAlpsTopo.ancmig_adj_3, Models_3D_BalkansApenninesCentralAlpsTopo.sim_split_no_mig, Models_3D_BalkansApenninesCentralAlpsTopo.sim_split_no_mig_size, Models_3D_BalkansApenninesCentralAlpsTopo.sim_split_sym_mig_all, Models_3D_BalkansApenninesCentralAlpsTopo.sim_split_sym_mig_all_size, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_isol, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_3epochs, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_4epochs_iter1, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_4epochs_iter2, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_4epochs_iter3, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_4epochs_iter4, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_simsplit_4epochs_iter5, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full_2_iter1, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full_2_iter2, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full_2_iter3, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full_2_iter4, Models_3D_BalkansApenninesCentralAlpsTopo.refugia_adj_5_full_2_iter5, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_3epochs_iter1, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_3epochs_iter2, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_3epochs_iter3, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_3epochs_iter4, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_3epochs_iter5, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_4epochs_iter1, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_4epochs_iter2, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_4epochs_iter3, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_4epochs_iter4, Models_3D_BalkansApenninesCentralAlpsTopo.split_full_4epochs_iter5, Models_3D_BalkansApenninesCentralAlpsTopo.split_simsplit_3epochs_iter1, Models_3D_BalkansApenninesCentralAlpsTopo.split_simsplit_3epochs_iter2, Models_3D_BalkansApenninesCentralAlpsTopo.split_simsplit_3epochs_iter3, Models_3D_BalkansApenninesCentralAlpsTopo.split_simsplit_3epochs_iter4, Models_3D_BalkansApenninesCentralAlpsTopo.split_simsplit_3epochs_iter5]
#ModelFunctions = [Models_3D_ApenninesBalkansCentralAlpsTopo.split_nomig, Models_3D_ApenninesBalkansCentralAlpsTopo.split_nomig_growth, Models_3D_ApenninesBalkansCentralAlpsTopo.split_nomig_growthall, Models_3D_ApenninesBalkansCentralAlpsTopo.split_symmig_all, Models_3D_ApenninesBalkansCentralAlpsTopo.split_symmig_all_growth, Models_3D_ApenninesBalkansCentralAlpsTopo.split_symmig_all_growthall, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_1, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_1_growth, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_1_growthall, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_2, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_3, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_4, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_4_nomig, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_4_growth, Models_3D_ApenninesBalkansCentralAlpsTopo.ancmig_adj_1, Models_3D_ApenninesBalkansCentralAlpsTopo.ancmig_adj_2, Models_3D_ApenninesBalkansCentralAlpsTopo.ancmig_adj_3, Models_3D_ApenninesBalkansCentralAlpsTopo.sim_split_no_mig, Models_3D_ApenninesBalkansCentralAlpsTopo.sim_split_no_mig_size, Models_3D_ApenninesBalkansCentralAlpsTopo.sim_split_sym_mig_all, Models_3D_ApenninesBalkansCentralAlpsTopo.sim_split_sym_mig_all_size, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_simsplit_isol, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_simsplit, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_simsplit_3epochs, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_full, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_simsplit_4epochs, Models_3D_ApenninesBalkansCentralAlpsTopo.refugia_adj_5_full_2]

### Model parameters (names per model)
params_split_nomig = "nu1, nuA, nu2, nu3, T1, T2"
upper_split_nomig = [100, 1000, 100, 100, 1, 1]

params_split_nomig_growth = "nu1, nuA0, nuA, nu2, nu3, T1, T2"
upper_split_nomig_growth = [100, 1000, 1000, 100, 100, 1, 1]

params_split_nomig_growthall = "nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, T1, T2"
upper_split_nomig_growthall = [100, 100, 1000, 1000, 100, 100, 100, 100, 1, 1]

params_split_symmig_all = "nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2"
upper_split_symmig_all = [100, 1000, 100, 100, 100, 100, 100, 100, 1, 1]

params_split_symmig_all_growth = "nu1, nuA0, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2"
upper_split_symmig_all_growth = [100, 1000, 1000, 100, 100, 100, 100, 100, 100, 1, 1]

params_split_symmig_all_growthall = "nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2"
upper_split_symmig_all_growthall = [100, 100, 1000, 1000, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1]

params_refugia_adj_1 = "nu1, nuA, nu2, nu3, m3_1, m3_2, m3_3, T1, T2, T3"
upper_refugia_adj_1 = [100, 1000, 100, 100, 100, 100, 100, 1, 1, 1]

params_refugia_adj_1_growth = "nu1, nuA0, nuA, nu2, nu3, m3_1, m3_2, m3_3, T1, T2, T3"
upper_refugia_adj_1_growth = [100, 1000, 1000, 100, 100, 100, 100, 100, 1, 1, 1]

params_refugia_adj_1_growthall = "nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, m3_1, m3_2, m3_3, T1, T2, T3"
upper_refugia_adj_1_growthall = [100, 100, 1000, 1000, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1]

params_refugia_adj_2 = "nu1, nuA, nu2, nu3, m2_1, m2_2, m2_3, T1, T2"
upper_refugia_adj_2 = [100, 1000, 100, 100, 100, 100, 100, 1, 1]

params_refugia_adj_3 = "nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1a, T1b, T2"
upper_refugia_adj_3 = [100, 1000, 100, 100, 100, 100, 100, 100, 1, 1, 1]

params_refugia_adj_4 = "nu1, nuA, nu2, nu3, m1_1, m3_2, m3_3, T1a, T1b, T2, T3"
upper_refugia_adj_4 = [100, 1000, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_refugia_adj_4_nomig = "nu1, nuA, nu2, nu3, m1_1, T1a, T1b, T2"
upper_refugia_adj_4_nomig = [100, 1000, 100, 100, 100, 1, 1, 1]

params_refugia_adj_4_growth = "nu1, nuA0, nuA, nu2, nu3, m1_1, m3_2, m3_3, T1a, T1b, T2, T3"
upper_refugia_adj_4_growth = [100, 1000, 1000, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_ancmig_adj_1 = "nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2, T3"
upper_ancmig_adj_1 = [100, 1000, 100, 100, 100, 100, 100, 100, 1, 1, 1]

params_ancmig_adj_2 = "nu1, nuA, nu2, nu3, m1_1, T1, T2"
upper_ancmig_adj_2 = [100, 1000, 100, 100, 100, 1, 1]

params_ancmig_adj_3 = "nu1, nuA, nu2, nu3, m1_1, T1a, T1b, T2"
upper_ancmig_adj_3 = [100, 1000, 100, 100, 100, 1, 1, 1]

params_sim_split_no_mig = "nuA, nu1, nu2, nu3, T1"
upper_sim_split_no_mig = [1000, 100, 100, 100, 100]

params_sim_split_no_mig_size = "nuA, nu1a, nu1b, nu2a, nu2b, nu3a, nu3b, T1, T2"
upper_sim_split_no_mig_size = [1000, 100, 100, 100, 100, 100, 100, 100, 100]

params_sim_split_sym_mig_all = "nuA, nu1, nu2, nu3, m_1, m_2, m_3, T1"
upper_sim_split_sym_mig_all = [1000, 100, 100, 100, 100, 100, 100, 100]

params_sim_split_sym_mig_all_size = "nuA, nu1a, nu1b, nu2a, nu2b, nu3a, nu3b, m_1, m_2, m_3, T1, T2"
upper_sim_split_sym_mig_all_size = [1000, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]

params_refugia_adj_5 = "nu1_1, nu1_2, nuA, nu2, nu3, m1_12, m1_21, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1a, T1b, T2, T3"
upper_refugia_adj_5 = [100, 100, 1000, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_refugia_adj_5_simsplit_isol = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T2, T3"
upper_refugia_adj_5_simsplit_isol = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1]

params_refugia_adj_5_simsplit = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T2, T3"
upper_refugia_adj_5_simsplit = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1]

params_refugia_adj_5_simsplit_3epochs = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3"
upper_refugia_adj_5_simsplit_3epochs = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1]

params_refugia_adj_5_full = "nu1_1a, nu1_1b, nu1_2, nu1_3, nuA_a, nuA_b, nu2_2, nu2_3, nu3_2, nu3_3, m1_12, m1_21, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1a, T1b, T2, T3"
upper_refugia_adj_5_full = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_refugia_adj_5_simsplit_4epochs = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4"
upper_refugia_adj_5_simsplit_4epochs = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_refugia_adj_5_full_2 = "nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4"
upper_refugia_adj_5_full_2 = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1, 1, 1]

params_split_full_3epochs = "nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3"
upper_split_full_3epochs = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1]

params_split_full_4epochs = "nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3"
upper_split_full_4epochs = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1, 1]

params_split_simsplit_3epochs = "nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3"
upper_split_simsplit_3epochs = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 1, 1, 1]


### Model parameters (list)
plabels_list = [params_split_nomig, params_split_nomig_growth, params_split_nomig_growthall, params_split_symmig_all, params_split_symmig_all_growth, params_split_symmig_all_growthall, params_refugia_adj_1, params_refugia_adj_1_growth, params_refugia_adj_1_growthall, params_refugia_adj_2, params_refugia_adj_3, params_refugia_adj_4, params_refugia_adj_4_nomig, params_refugia_adj_4_growth, params_ancmig_adj_1, params_ancmig_adj_2, params_ancmig_adj_3, params_sim_split_no_mig, params_sim_split_no_mig_size, params_sim_split_sym_mig_all, params_sim_split_sym_mig_all_size, params_refugia_adj_5, params_refugia_adj_5_simsplit_isol, params_refugia_adj_5_simsplit, params_refugia_adj_5_simsplit_3epochs, params_refugia_adj_5_full, params_refugia_adj_5_simsplit_4epochs, params_refugia_adj_5_simsplit_4epochs, params_refugia_adj_5_simsplit_4epochs, params_refugia_adj_5_simsplit_4epochs, params_refugia_adj_5_simsplit_4epochs, params_refugia_adj_5_full_2, params_refugia_adj_5_full_2, params_refugia_adj_5_full_2, params_refugia_adj_5_full_2, params_refugia_adj_5_full_2, params_split_full_3epochs, params_split_full_3epochs, params_split_full_3epochs, params_split_full_3epochs, params_split_full_3epochs, params_split_full_4epochs, params_split_full_4epochs, params_split_full_4epochs, params_split_full_4epochs, params_split_full_4epochs, params_split_simsplit_3epochs, params_split_simsplit_3epochs, params_split_simsplit_3epochs, params_split_simsplit_3epochs, params_split_simsplit_3epochs]
upper_list = [upper_split_nomig, upper_split_nomig_growth, upper_split_nomig_growthall, upper_split_symmig_all, upper_split_symmig_all_growth, upper_split_symmig_all_growthall, upper_refugia_adj_1, upper_refugia_adj_1_growth, upper_refugia_adj_1_growthall, upper_refugia_adj_2, upper_refugia_adj_3, upper_refugia_adj_4, upper_refugia_adj_4_nomig, upper_refugia_adj_4_growth, upper_ancmig_adj_1, upper_ancmig_adj_2, upper_ancmig_adj_3, upper_sim_split_no_mig, upper_sim_split_no_mig_size, upper_sim_split_sym_mig_all, upper_sim_split_sym_mig_all_size, upper_refugia_adj_5, upper_refugia_adj_5_simsplit_isol, upper_refugia_adj_5_simsplit, upper_refugia_adj_5_simsplit_3epochs, upper_refugia_adj_5_full, upper_refugia_adj_5_simsplit_4epochs, upper_refugia_adj_5_simsplit_4epochs, upper_refugia_adj_5_simsplit_4epochs, upper_refugia_adj_5_simsplit_4epochs, upper_refugia_adj_5_simsplit_4epochs, upper_refugia_adj_5_full_2, upper_refugia_adj_5_full_2, upper_refugia_adj_5_full_2, upper_refugia_adj_5_full_2, upper_refugia_adj_5_full_2, upper_split_full_3epochs, upper_split_full_3epochs, upper_split_full_3epochs, upper_split_full_3epochs, upper_split_full_3epochs, upper_split_full_4epochs, upper_split_full_4epochs, upper_split_full_4epochs, upper_split_full_4epochs, upper_split_full_4epochs, upper_split_simsplit_3epochs, upper_split_simsplit_3epochs, upper_split_simsplit_3epochs, upper_split_simsplit_3epochs, upper_split_simsplit_3epochs]

### Model parameter length
num_params_list = [(p.count(',') + 1) for p in plabels_list]

### Run Moments optimization script
for model in range(model_num-1, model_num):
    #print(modelname_list[model])
    for chain in range(1,chains+1):   # To run the optimization routine multiple times
        prefix_c = prefix + "_{}".format(chain)
        Optimize_Functions.Optimize_Routine(data, prefix_c, modelname_list[model], ModelFunctions[model], rounds, num_params_list[model], fs_folded=fs_folded, reps=reps, maxiters=maxiters, folds=folds, in_upper=upper_list[model], param_labels=plabels_list[model])
   
