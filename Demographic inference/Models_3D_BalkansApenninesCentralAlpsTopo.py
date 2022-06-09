#!/usr/bin/python

### This script contains 2-population model definitions that reflect potential Central Alps recolonisation scenarios.

# Hirzi Luqman
# Feb 2019

# Import modules for running Moments
#import os
import numpy
import scipy
import matplotlib
import matplotlib.pyplot as pyplot
import cython
import sys
#sys.path.append("/Users/luqman/Documents/Hirzi/ETHPHD/Moments")
import moments
#from dadi import Numerics, PhiManip, Integration
#from dadi.Spectrum_mod import Spectrum
#import pylab
import time

## Set working directory
#os.chdir("/Users/luqman/Desktop/Demographic_Inference_Moments/Moments analysis")

#%%
'''
Models for testing three population scenarios.
'''


#def test_demo((nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs), (n1,n2,n3)):
#def test_demo(nuAf, nuB, nuEu0, nuEu, nuAs0, nuAs, mAfB, mAfEu, mAfAs, mEuAs, TAf, TB, TEuAs, ns):
def test_demo(params, ns):
    # Define parameters
    nu0, nu1, nu2, nu12, T1, T2 = params
    # Start  with an equilibrium population of non-zero size.
    #sts = moments.LinearSystem_1D.steady_state_1D(n1+n2+n3)
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    # Make a spectrum object
    fs = moments.Spectrum(sts)
    # Size change in this ancestral population
    #fs.integrate([nuAf], TAf, 0.05)
    # Population split 1
    #fs = moments.Manips.split_1D_to_2D(fs, n1, n2+n3)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    # Define symmetric migration
    #mig1 = numpy.array([[0, mAfB],[mAfB, 0]]) 
    # Define assymmetric migration
    #mig1 = numpy.array([[0, mAfB],[mBAf, 0]])
    # Integrate over time dependent continuous process (e.g. genetic drift at different population sizes or migrations)    
    #fs.integrate([nuAf, nuB], TB, 0.05, m=mig1)
    fs.integrate([nu0, nu12], T1)
    # Population split 2
    #fs = moments.Manips.split_2D_to_3D_2(fs, n2, n3)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    # Define exponential growth in population 1
    #nuEu_func = lambda t: nuEu0*(nuEu/nuEu0)**(t/TEuAs)
    # Define xponential growth in population 2
    #nuAs_func = lambda t: nuAs0*(nuAs/nuAs0)**(t/TEuAs)
    # Define Poverall population size change function
    #nu2 = lambda t: [nuAf, nuEu_func(t), nuAs_func(t)]
    # Define 3 population migration matrix
    #mig2 = numpy.array([[0, mAfEu, mAfAs], [mAfEu, 0, mEuAs], [mAfAs, mEuAs, 0]]) 
    # Intergrate over population size change function and migration matrix
    #fs.integrate(nu2, TEuAs, 0.05, m=mig2)
    fs.integrate([nu0, nu1, nu2], T2)
    return fs


def split_nomig(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #6 parameters	
    nu1, nuA, nu2, nu3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function for T1
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)    
    return fs


def split_nomig_growth(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu1: Size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #7 parameters	
    nu1, nuA0, nuA, nu2, nu3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    fs.integrate(nu_T1_func, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)    
    return fs


def split_nomig_growthall(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    nu10: Initial size of population 1 after split.
    nu1: Final size of population 1 after split.   
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu20: Initial size of population 2 after split.
    nu2: Finall size of population 2 after split.   
    nu30: Initial size of population 3 after split.
    nu3: Final size of population 3 after split.   
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #10 parameters	
    nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    fs.integrate(nu_T1_func, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T2
    nu1_func = lambda t: nu10 * (nu1/nu10)**(t/T2) 
    nu2_func = lambda t: nu20 * (nu2/nu20)**(t/T2)    
    nu3_func = lambda t: nu30 * (nu3/nu30)**(t/T2)   
    nu_T2_func = lambda t: [nu1_func(t), nu2_func(t), nu3_func(t)]
    fs.integrate(nu_T2_func, T2)    
    return fs


def split_symmig_all(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #10 parameters 
    nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1, nuA]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig2)    
    return fs


def split_symmig_all_growth(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu1: Size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #11 parameters 
    nu1, nuA0, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    ## Population function for T1
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1_func, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig2)    
    return fs


def split_symmig_all_growthall(params, ns):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    nu10: Initial size of population 1 after split.
    nu1: Final size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu20: Initial size of population 2 after split.
    nu2: Finall size of population 2 after split.   
    nu30: Initial size of population 3 after split.
    nu3: Final size of population 3 after split.  
    m1_1: Migration rate between population 1 and population (2,3)
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #14 parameters 
    nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    ## Population function for T1
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1_func, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu1_func = lambda t: nu10 * (nu1/nu10)**(t/T2) 
    nu2_func = lambda t: nu20 * (nu2/nu20)**(t/T2)    
    nu3_func = lambda t: nu30 * (nu3/nu30)**(t/T2)   
    nu_T2_func = lambda t: [nu1_func(t), nu2_func(t), nu3_func(t)]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2_func, T2, m=mig2)    
    return fs


def refugia_adj_1(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, gene flow does not occur. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
    'longest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m3_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m3_2: Migration rate between populations 2 and 3 
    m3_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #10 parameters 
    nu1, nuA, nu2, nu3, m3_1, m3_2, m3_3, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1, nu2, nu3]
    mig3 = numpy.array([[0, m3_1, m3_3],[m3_1, 0, m3_2], [m3_3, m3_2, 0]])  
    fs.integrate(nu_T3, T3, m=mig3)   
    return fs


def refugia_adj_1_growth(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, gene flow does not occur. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
    'longest isolation'
    nu1: Size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m3_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m3_2: Migration rate between populations 2 and 3 
    m3_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #11 parameters 
    nu1, nuA0, nuA, nu2, nu3, m3_1, m3_2, m3_3, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    ## Population function for T1
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    fs.integrate(nu_T1_func, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1, nu2, nu3]
    mig3 = numpy.array([[0, m3_1, m3_3],[m3_1, 0, m3_2], [m3_3, m3_2, 0]])  
    fs.integrate(nu_T3, T3, m=mig3)   
    return fs


def refugia_adj_1_growthall(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, gene flow does not occur. Period of symmetric secondary contact occurs between 
    adjacent populations (ie 1<->2, 2<->3, but not 1<->3) after all splits are complete. 
    'longest isolation'
    nu10: Initial size of population 1 after split.
    nu1: Final size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu20: Initial size of population 2 after split.
    nu2: Finall size of population 2 after split.   
    nu30: Initial size of population 3 after split.
    nu3: Final size of population 3 after split.
    m3_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m3_2: Migration rate between populations 2 and 3 
    m3_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #14 parameters 
    nu10, nu1, nuA0, nuA, nu20, nu2, nu30, nu3, m3_1, m3_2, m3_3, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1) 
    ## Population function for T1
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    fs.integrate(nu_T1_func, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu1_func = lambda t: nu10 * (nu1/nu10)**(t/T3) 
    nu2_func = lambda t: nu20 * (nu2/nu20)**(t/T3)    
    nu3_func = lambda t: nu30 * (nu3/nu30)**(t/T3)   
    nu_T3_func = lambda t: [nu1_func(t), nu2_func(t), nu3_func(t)] 
    mig3 = numpy.array([[0, m3_1, m3_3],[m3_1, 0, m3_2], [m3_3, m3_2, 0]])  
    fs.integrate(nu_T3_func, T3, m=mig3)   
    return fs


def refugia_adj_2(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur. Split between pops
    2 and 3, with gene flow. After appearance of 2 and 3, gene flow also occurs between 1 
    and 2.
    'shorter isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #9 parameters 
    nu1, nuA, nu2, nu3, m2_1, m2_2, m2_3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig2)
    return fs


def refugia_adj_3(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3 occurs with gene flow, and gene flow
    happens between 1 and 2 as well.
    'shortest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #11 parameters 
    nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1a, T1b, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1a)
    ## Population function and migration matrix for T1b
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig2)
    return fs


def refugia_adj_4(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3, gene flow does not occur, but then 
    secondary contact occurs from pop 1 to pop 3, and pop 2 to pop 3.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
	m3_2: Migration rate between populations 2 and 3 
    m3_3: Migration rate between populations 1 and 3
    T1a: The scaled time from the split of pops 1 vs 2 and 3, and the end of initial isolation (in units of 2*Na generations).
    T1b: The scaled time from the end of isolation (secondary contact) to the second split (in units of 2*Na generations).
    T2: The scaled time from the second split following which populations are isolated (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #11 parameters 
    nu1, nuA, nu2, nu3, m1_1, m3_2, m3_3, T1a, T1b, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1a)
    ## Population function and migration matrix for T1b
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1, nu2, nu3]
    mig3 = numpy.array([[0, 0, 0],[0, 0, 0], [m3_3, m3_2, 0]])  
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def refugia_adj_4_nomig(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3, gene flow does not occur.
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    T1a: The scaled time from the split of pops 1 vs 2 and 3, and the end of initial isolation (in units of 2*Na generations).
    T1b: The scaled time from the end of isolation (secondary contact) to the second split (in units of 2*Na generations).
    T2: The scaled time from the second split following which populations are isolated (in units of 2*Na generations).
    """
    #8 parameters 
    nu1, nuA, nu2, nu3, m1_1, T1a, T1b, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1a)
    ## Population function and migration matrix for T1b
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    return fs


def refugia_adj_4_growth(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs. Split between pops 2 and 3, gene flow does not occur.
    nu1: Size of population 1 after split.
    nuA0: Size of population (2,3) directly after split from 1.
    nuA: Size of population (2,3) before it splits into populations 2 and 3.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    T1a: The scaled time from the split of pops 1 vs 2 and 3, and the end of initial isolation (in units of 2*Na generations).
    T1b: The scaled time from the end of isolation (secondary contact) to the second split (in units of 2*Na generations).
    T2: The scaled time from the second split following which populations are isolated (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #12 parameters 
    nu1, nuA0, nuA, nu2, nu3, m1_1, m3_2, m3_3, T1a, T1b, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nuA_func = lambda t: nuA0 * (nuA/nuA0)**(t/T1a) 
    nu_T1_func = lambda t: [nu1, nuA_func(t)]
    fs.integrate(nu_T1_func, T1a)
    ## Population function and migration matrix for T1b
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    nu_T1 = [nu1, nuA]
    fs.integrate(nu_T1, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1, nu2, nu3]
    mig3 = numpy.array([[0, 0, 0],[0, 0, 0], [m3_3, m3_2, 0]])  
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def ancmig_adj_1(params, ns):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3 with gene flow, then all gene flow ceases.
    'shortest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    m2_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m2_2: Migration rate between populations 2 and 3 
    m2_3: Migration rate between populations 1 and 3
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #11 parameters 
    nu1, nuA, nu2, nu3, m1_1, m2_1, m2_2, m2_3, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1, nuA]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    mig2 = numpy.array([[0, m2_1, m2_3],[m2_1, 0, m2_2], [m2_3, m2_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1, nu2, nu3]
    fs.integrate(nu_T3, T3)   
    return fs

def ancmig_adj_2(params, ns):
    """
    Model with split between pop 1 and (2,3), with gene flow. Split 
    between pops 2 and 3, and all gene flow ceases.
    'shorter isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #7 parameters 
    nu1, nuA, nu2, nu3, m1_1, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1, nuA]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    return fs

def ancmig_adj_3(params, ns):
    """
    Model with split between pop 1 and (2,3), with gene flow, which then stops. Split 
    between pops 2 and 3, gene flow does not occur at all.
    'longest isolation'
    nu1: Size of population 1 after split.
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_1: Migration rate between population 1 and population (2,3)
    T1a: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T1b: The scaled time for no gene flow between pops 1 and (2, 3) (in units of 2*Na generations).
    T2: The scaled time between the cessation of gene flow and the split of pops 2 and 3 (in units of 2*Na generations).
    """
    #8 parameters 
    nu1, nuA, nu2, nu3, m1_1, T1a, T1b, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1 = [nu1, nuA]
    mig1 = numpy.array([[0, m1_1],[m1_1, 0]])
    fs.integrate(nu_T1, T1a, m=mig1)
    fs.integrate(nu_T1, T1b)    
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1, nu2, nu3]
    fs.integrate(nu_T2, T2)
    return fs


def sim_split_no_mig(params, ns):
    """
    Model with simulataneous splitting of 3 populations. 
    Migration does not occur between any population pair.
    nuA: Size of ancestral population prior to split
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the population split and the present (in units of 2*Na generations).
    """
    #5 parameters	
    nuA, nu1, nu2, nu3, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T1
    nu_T1 = [nu1, nu2, nu3]
    fs.integrate(nu_T1, T1)    
    return fs


def sim_split_no_mig_size(params, ns):
    """
    Model with simulataneous splitting of 3 populations, followed by single discrete size change event.
    Migration does not occur between any population pair.
    nuA: Size of ancestral population prior to split
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    T1: The scaled time between the population split and the discrete population size change event (in units of 2*Na generations).
    T2: The scaled time between the discrete population size change event and the present (in units of 2*Na generations).
    """
    #9 parameters	
    nuA, nu1a, nu1b, nu2a, nu2b, nu3a, nu3b, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    fs.integrate(nu_T1, T1)
    ## Population function for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)   
    return fs


def sim_split_sym_mig_all(params, ns):
    """
    Model with simulataneous splitting of 3 populations.
    Symmetric migration allowed between all pairs.
    nuA: Size of ancestral population prior to split
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m_2: Migration rate between populations 2 and 3 
    m_3: Migration rate between populations 1 and 3
    T1: The scaled time between the population split and the present (in units of 2*Na generations).
    """
    #8 parameters	
    nuA, nu1, nu2, nu3, m_1, m_2, m_3, T1 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T1
    nu_T1 = [nu1, nu2, nu3]
    mig = numpy.array([[0, m_1, m_3],[m_1, 0, m_2], [m_3, m_2, 0]])  
    fs.integrate(nu_T1, T1, m=mig) 
    return fs


def sim_split_sym_mig_all_size(params, ns):
    """
    Model with simulataneous splitting of 3 populations, followed by single discrete size change event.
    Symmetric migration allowed between all pairs.
    nuA: Size of ancestral population prior to split
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m_1: Migration rate between populations 1 and 2 (2*Na*m) 
	m_2: Migration rate between populations 2 and 3 
    m_3: Migration rate between populations 1 and 3
    T1: The scaled time between the population split and the discrete population size change/secondary contact event (in units of 2*Na generations).
    T2: The scaled time between the discrete population size change/secondary contact event and the present (in units of 2*Na generations).
    """
    #12 parameters	
    nuA, nu1a, nu1b, nu2a, nu2b, nu3a, nu3b, m_1, m_2, m_3, T1, T2 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    fs.integrate(nu_T1, T1)
    ## Population function for T2    
    nu_T2 = [nu1b, nu2b, nu3b]
    mig = numpy.array([[0, m_1, m_3],[m_1, 0, m_2], [m_3, m_2, 0]])  
    fs.integrate(nu_T2, T2, m=mig) 
    return fs


def refugia_adj_5(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs (w/ assymmetric migration). Split between pops 2 and 3, gene flow does not occur, but then 
    secondary contact occurs between all pairs (assymmetric migration)
    nu1_1: Size of population 1 after split, at time epoch T1.
    nu1_2: Size of population 1 after split, at time epoch T1.   
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_12: Migration rate from pop 1 to pop 2, during time epoch T1b.
    m1_21: Migration rate from pop 2 to pop 1, during time epoch T1b.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1a: The scaled time from the split of pops 1 vs 2 and 3, and the end of initial isolation (in units of 2*Na generations).
    T1b: The scaled time from the end of isolation (secondary contact) to the second split (in units of 2*Na generations).
    T2: The scaled time from the second split following which populations are isolated (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #17 parameters 
    nu1_1, nu1_2, nuA, nu2, nu3, m1_12, m1_21, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1a, T1b, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1 = [nu1_1, nuA]
    fs.integrate(nu_T1, T1a)
    ## Population function and migration matrix for T1b
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1_2, nu2, nu3]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1_2, nu2, nu3]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def refugia_adj_5_simsplit_isol(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Gene flow does not occur immediately after split, but then 
    secondary contact occurs (assymmetric migration)
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T2: The scaled time from the trichotomous split to end of period of initial isolation (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #14 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1a, nu2a, nu3a]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1b, nu2b, nu3b]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def refugia_adj_5_simsplit(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Two time epochs; gene flow (assymmetric migration) and population size change can occur 
    can occur in both.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T2: The scaled time from the trichotomous split to end of period of initial isolation (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #20 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1a, nu2a, nu3a]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1b, nu2b, nu3b]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def refugia_adj_5_simsplit_3epochs(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Three time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch.
    can occur in both.
    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T2: The scaled time from the trichotomous split to end of period of initial isolation (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


def refugia_adj_5_full(params, ns):
    """
    Model with split between pop 1 and (2,3), gene flow does not occur, but then 
    secondary contact occurs (w/ assymmetric migration & popullation size change). Split between pops 2 and 3, gene flow does not occur, but then 
    secondary contact occurs between all pairs (assymmetric migration & population size change)
    nu1_1: Size of population 1 after split, at time epoch T1.
    nu1_2: Size of population 1 after split, at time epoch T1.   
    nuA: Size of population (2,3) after split from 1.
    nu2: Size of population 2 after split.
    nu3: Size of population 3 after split.
    m1_12: Migration rate from pop 1 to pop 2, during time epoch T1b.
    m1_21: Migration rate from pop 2 to pop 1, during time epoch T1b.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1a: The scaled time from the split of pops 1 vs 2 and 3, and the end of initial isolation (in units of 2*Na generations).
    T1b: The scaled time from the end of isolation (secondary contact) to the second split (in units of 2*Na generations).
    T2: The scaled time from the second split following which populations are isolated (in units of 2*Na generations).
    T3: The scaled time between end of isolation (onset of secondary contact) and the present. (in units of 2*Na generations).
    """
    #22 parameters 
    nu1_1a, nu1_1b, nu1_2, nu1_3, nuA_a, nuA_b, nu2_2, nu2_3, nu3_2, nu3_3, m1_12, m1_21, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1a, T1b, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1a
    nu_T1a = [nu1_1a, nuA_a]
    fs.integrate(nu_T1a, T1a)
    ## Population function and migration matrix for T1b
    nu_T1b = [nu1_1b, nuA_b]    
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1b, T1b, m=mig1)     
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1_2, nu2_2, nu3_2]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3
    nu_T3 = [nu1_3, nu2_3, nu3_3]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3) 
    return fs


### Models below with "iter" suffix indicate copies of same model. This allows easy parallelisation of these models (which have many parameters and thus slow to run)
def refugia_adj_5_simsplit_4epochs_iter1 (params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #28 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_simsplit_4epochs_iter2 (params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #28 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_simsplit_4epochs_iter3 (params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #28 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_simsplit_4epochs_iter4 (params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #28 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_simsplit_4epochs_iter5 (params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #28 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_full_2_iter1 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #33 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_full_2_iter2 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #33 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_full_2_iter3 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #33 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_full_2_iter4 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #33 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def refugia_adj_5_full_2_iter5 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #33 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, nu1d, nu2d, nu3d, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m3_12, m3_13, m3_21, m3_23, m3_31, m3_32, T0, T1, T2, T3, T4 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1 (to reflect sum effect of all previous glacial-interglacial cycles)
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2 (to reflect period of isolation during last glacial)
    nu_T2 = [nu1b, nu2b, nu3b]
    fs.integrate(nu_T2, T2)
    ## Population function and migration matrix for T3 (to reflect inter-glacial expansion)
    nu_T3 = [nu1c, nu2c, nu3c]
    mig3 = numpy.array([[0, m3_12, m3_13],[m3_21, 0, m3_23], [m3_31, m3_32, 0]]) 
    fs.integrate(nu_T3, T3, m=mig3)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T4 = [nu1d, nu2d, nu3d]
    fs.integrate(nu_T4, T4) 
    return fs


def split_full_3epochs_iter1 (params, ns):
    """
    Model with split between pop 1 and (2,3), then between pops 2 and 3, with asymmetric gene flow.
    Final epoch given to capture expected bottleneck effect expected when having a single population represent lineage (no migration).
    nuA: Size of population (2,3) after split from 1.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #19 parameters 
    nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nuA]
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3)   
    return fs


def split_full_3epochs_iter2 (params, ns):
    """
    Model with split between pop 1 and (2,3), then between pops 2 and 3, with asymmetric gene flow.
    Final epoch given to capture expected bottleneck effect expected when having a single population represent lineage (no migration).
    nuA: Size of population (2,3) after split from 1.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #19 parameters 
    nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nuA]
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3)   
    return fs


def split_full_3epochs_iter3 (params, ns):
    """
    Model with split between pop 1 and (2,3), then between pops 2 and 3, with asymmetric gene flow.
    Final epoch given to capture expected bottleneck effect expected when having a single population represent lineage (no migration).
    nuA: Size of population (2,3) after split from 1.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #19 parameters 
    nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nuA]
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3)   
    return fs


def split_full_3epochs_iter4 (params, ns):
    """
    Model with split between pop 1 and (2,3), then between pops 2 and 3, with asymmetric gene flow.
    Final epoch given to capture expected bottleneck effect expected when having a single population represent lineage (no migration).
    nuA: Size of population (2,3) after split from 1.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #19 parameters 
    nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nuA]
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3)   
    return fs


def split_full_3epochs_iter5 (params, ns):
    """
    Model with split between pop 1 and (2,3), then between pops 2 and 3, with asymmetric gene flow.
    Final epoch given to capture expected bottleneck effect expected when having a single population represent lineage (no migration).
    nuA: Size of population (2,3) after split from 1.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time between the split of pops 1 vs 2 and 3 (in units of 2*Na generations).
    T2: The scaled time between the first split and the split of pops 2 and 3 (in units of 2*Na generations).
    T3: The scaled time between the split of pops 2 and 3 and the present (in units of 2*Na generations).	
    """
    #19 parameters 
    nu1a, nuA, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_21, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nuA]
    mig1 = numpy.array([[0, m1_12],[m1_21, 0]])
    fs.integrate(nu_T1, T1, m=mig1)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
     ## Population function and migration matrix for T3
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3)   
    return fs


def split_full_4epochs_iter1 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #29 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3  = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_full_4epochs_iter2 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #29 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3  = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_full_4epochs_iter3 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #29 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3  = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_full_4epochs_iter4 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #29 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3  = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_full_4epochs_iter5 (params, ns):
    """
    Model with split between pop 1 and (2,3), and then between pops 2 and 3. This model is identical to the refugia_adj_5_simsplit_4epochs model, but including a separate time 
    epoch (with assymmetric migration) between the first and second splits.
    expected when having a single population represent lineage.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T0: The scaled time from the first split to the second split (in units of 2*Na generations).    
    T1: The scaled time from the second split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).   
    """
    #29 parameters 
    nu1x, nuA, nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m0_12, m0_21, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T0, T1, T2, T3  = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])
    ## Population function and migration matrix for T0 (initial split; the definition of this time epoch differentiates this model from refugia_adj_5_simsplit_4epochs)
    nu_T0 = [nu1x, nuA]
    mig0 = numpy.array([[0, m0_12],[m0_21, 0]])
    fs.integrate(nu_T0, T0, m=mig0)
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_simsplit_3epochs_iter1(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_simsplit_3epochs_iter2(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_simsplit_3epochs_iter3(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_simsplit_3epochs_iter4(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs


def split_simsplit_3epochs_iter5(params, ns):
    """
    Model with trichomtomous split resulting simultaneously with pops 1, 2 and 3. Four time epochs; gene flow (assymmetric migration) occurs in the first and third
    time epoch, but not in the second (which represents a period of isolation). Population size change can occur at each time epoch. Final epoch given to capture expected bottleneck effect
    expected when having a single population represent lineage. This model is a nested model within the refugia_adj_5_full_2 below.
    nu1: Size of population 1 at indicated time epoch.
    nu2: Size of population 2 at indicated time epoch.
    nu3: Size of population 3 at indicated time epoch.
    m3_12: Migration rate from pop 1 to pop 2, during time epoch T3.
    m3_21: Migration rate from pop 2 to pop 1, during time epoch T3.
    mT_ij: Migration rate from pop j to pop i, during time epoch TT
    T1: The scaled time from the trichotomous split to the start of the last glacial, to represent sum effect of all previous glacial-interglacial cycles  (in units of 2*Na generations).
    T2: The scaled time from the start to the end of the last glacial, to represent period of isolation during last glacial (in units of 2*Na generations).    
    T3: The scaled time from the start to the end of the last interglacial, to represent inter-glacial expansion (in units of 2*Na generations).
    T4: The scaled time near the present, to represent bottleneck to capture single population representation of lineage due to sampling (in units of 2*Na generations).    
    """
    #24 parameters 
    nu1a, nu2a, nu3a, nu1b, nu2b, nu3b, nu1c, nu2c, nu3c, m1_12, m1_13, m1_21, m1_23, m1_31, m1_32, m2_12, m2_13, m2_21, m2_23, m2_31, m2_32, T1, T2, T3 = params
    sts = moments.LinearSystem_1D.steady_state_1D(ns[0] + ns[1] + ns[2])
    fs = moments.Spectrum(sts)
    fs = moments.Manips.split_1D_to_2D(fs, ns[0], ns[1] + ns[2])  
    fs = moments.Manips.split_2D_to_3D_2(fs, ns[1], ns[2])
    ## Population function and migration matrix for T1
    nu_T1 = [nu1a, nu2a, nu3a]
    mig1 = numpy.array([[0, m1_12, m1_13],[m1_21, 0, m1_23], [m1_31, m1_32, 0]]) 
    fs.integrate(nu_T1, T1, m=mig1)
    ## Population function and migration matrix for T2
    nu_T2 = [nu1b, nu2b, nu3b]
    mig2 = numpy.array([[0, m2_12, m2_13],[m2_21, 0, m2_23], [m2_31, m2_32, 0]]) 
    fs.integrate(nu_T2, T2, m=mig2)
    ## Population function and migration matrix for T3 (bottleneck to capture single population representation of lineage)
    nu_T3 = [nu1c, nu2c, nu3c]
    fs.integrate(nu_T3, T3) 
    return fs

#%%