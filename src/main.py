# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


 # Import module
import os
path = '/home/nguyetrt/rbo-rank-genelist/src/'
#path="C:\\Users\\ttran\\OneDrive - Indiana University\\2022-AsusZenbook14-sync\\Documents\\2022 - FALL SEMESTER IU\\2022_I552 Independent studies\\RBO\\rbo-rank-genelist\\src\\"

os.chdir(path)

import numpy as np
import argparse  # check input arguments
import pandas as pd  # handle dataframe
import math   # display progress in simulation
import sys


from random import sample


from generatedata import comparegeneratedlist
from generatedata import simulation
#from generatedata import generatedatadefault
#from generatedata import similarityscorealldata as SMD


from scipy.stats import skellam
from scipy.stats import triang 
from scipy.stats import shapiro

from rbo import rbo_modified as SM
from rbo import wg_geom as wgm, wg_binomial as wbi, wg_poisson as wpo, wg_skellam as wsk, wg_triangular as wtr



from distfit import distfit
import fitter

import matplotlib.pyplot as plt

# Set working directory

if __name__ == '__main__':

    # Add some arguments
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-nbXk', type=int, 
                        default = '10', 
                        help="Number of samples to generate from initial sample (positive INT)")
    
    parser.add_argument('-nbrep', type=int,
                        default = '10', 
                        help="Number of replicates to generate for each sample \
                            (positive INT)")
    
    parser.add_argument('-N', type=int, 
                        default = '1000', 
                        help="Total number of indiv/ gene for sample generation \
                            (positive INT)")
    
    parser.add_argument('-cat', 
                        type=str, 
                        default = 'H', 
                        help="Category. H or MH or M or LM or L")
    
    parser.add_argument('-expid', 
                        type=str, 
                        default = '2', 
                        help="The id of experiment to test on (STR)")
    
    parser.add_argument('-topk', 
                        type=int, 
                        default = '10', 
                        help="Top k first elements within the \
                            category to check(positive INT)")

    args = parser.parse_args()









#------------------------------------------------------------------------
# RETAIN ARGUMENTS

    nbXk = args.nbXk # Number of experiments k
    nbrep = args.nbrep  # Number of replicates from each Xk
    N = args.N
    condition = args.cat
    expid = args.expid
    topk = args.topk

    
#------------------------------------------------------------------------

# FIX VARIABLES


    seed = 0
    list_Xk = []  # list of all transition matrix for each Xk
    
    # Define weighting schemes
    mu1 = 2*topk
    mu2 = topk
    geomw = {'geometric': None, 'p': 0.99}
    binomw = {'binomial': None, 'n': 20,'p': 0.5}
    poisw = {'poisson': None, 'lambda': 10}
    skelw = {'skellam': None, 'mu1': mu1,'mu2': mu2}
    triaw = {'triangular': None, 'topk': topk}
    
    pgeom = geomw['p']
    n = binomw['n']
    pbino = binomw['p']
    lamb = poisw['lambda']
    
    
    
    
    # Check how different each expr is from baseline condition
    group0 = 'Condition_0_state'
    Avgswitch = {}  # dictionary for number of elements that were switched places
    
    

#------------------------------------------------------------------------
    #  Call function to generate transition matrices for each Xk
    #dataframe = generatedatadefault(nbXk, nbrep, N, sd=seed)
    # Alternatively, call an already generated dataframe
    path = os.path.dirname(os.getcwd())
    dataframe = pd.read_csv(path+'/data/testdata.csv', index_col=0)


    # Check how different each expr is from baseline condition
    print('-'*50)
    
    for k in range(1, nbrep+1): 
        group1 = f'Condition_{k}_state'
        Avgswitch[group1] = comparegeneratedlist(dataframe, group0, group1)
        print(f'Condition_{k}:')
        print(f'Average number of modified gene expressions: {np.round(Avgswitch[group1], 2)}')
        print()
    
    
    
# --------------------------------------------------------------------------
     
    # Measure similarity between two lists using different metrics
    # Compare different similarity scores
     
    # For RBO scores
    dissim_pos = [] # List of first disagrreement position
    c = (dataframe[f'Condition_{expid}_state'] == condition)
    l = list(dataframe[c].sort_values(f'Condition_{expid}', ascending=False).index)
     
    print()
    print('-'*70)
    print(f'Looking at Condition_{expid}')
    print('-'*70)
    print()
     
     
    for i in range(1, nbrep+1):
         
         print()
         print('-'*50)
         print(f'For Condition_{expid}_prime_{i}')
         print(f'Top {topk} ranked genes')
         print('.'*30)
         cprime = (dataframe[f'Condition_{expid}_prime_state_{i}'] == condition)    
         lprime = list(dataframe[cprime].sort_values(f'Condition_{expid}_prime_{i}', 
                                                    ascending=False).index)
        
         # Display the topk elements
         for k in range(topk):
             print(f'Rank {k+1} : {l[k]} --- {lprime[k]} --> Equal: {l[k] == lprime[k]} ')
         
         dissim = topk+5
         for k in range(0, topk):
             if (l[k] != lprime[k]):
                 dissim = k+1
                 break
         dissim_pos.append(dissim)

         print('-'*50)
         print()
         
         geom = np.round(SM(wgm, geomw, l, lprime), 3)
         binom = np.round(SM(wbi, binomw, l, lprime), 3)
         pois = np.round(SM(wpo, poisw, l, lprime), 3)
         skel = np.round(SM(wsk, skelw, l, lprime), 3)
         tria = np.round(SM(wtr, triaw, l, lprime), 3)
         
         print(f'score rbo (p={pgeom}) =', geom)    
         print(f'score binomial (n={n}, p={pbino}) =', binom)
         print(f'score poisson (lamb={lamb})=', pois)
         print(f'score skellam (mu1={mu1}, mu2={mu2}) =', skel)
         print(f'score triangular (max_pos={topk})=', tria)

         
# ----------------------------------------------------------------------------
    # Generate dataframes of similarity scores
    
    # List of weighting functions, parameters
    weightfuncs = [wgm, wbi, wpo, wsk, wtr]
    weightparameters = [geomw, binomw, poisw, skelw, triaw]
    weightnames = ["GeomSim.csv", "BinomSim.csv", "PoissonSim.csv", \
                   "SkellamSim.csv", "TriangSim.csv"]
    # WARNING!
    # Run these lines will take time, rather import data directly
    # RUN THIS: Import data from already generated scores
    #print('-'*50)
    #for i in range(len(weightfuncs)):
    #    print('Computing similarity score for data...')
    #    print(f'Weighting scheme: {weightparameters[i]}')
    #    SMD(nbXk, nbrep, dataframe, condition, weightfuncs[i], weightparameters[i],\
    #        weightnames[i])
    #    print()
    df = pd.read_csv(path+'/data/SkellamSim.csv')
    
    
# -----------------------------------------------------------------------------
    # Simulate distribution
# =============================================================================
#     dist = distfit(todf=True)
#     output, output_scores = simulation(withreplacement=True)
#     results = dist.fit_transform(np.array(output_scores))
#     dist.plot()
#     dist.plot_summary()
#     
#     f = fitter.Fitter(output_scores)
#     f.fit()
#     f.summary()
# =============================================================================
    for i in range(10, 51, 10):
        for j in [0.05, 0.2, 0.5, 1.]:
            
            output, output_scores = simulation(withreplacement=True, size_l=i, percent_repl=j)
            dist = distfit(todf=True)
            results = dist.fit_transform(np.array(output_scores))
            dist.plot()
            #dist.plot_summary()
            
            f = fitter.Fitter(output_scores)
            f.fit()
            f.summary()

        
             
        
    
    
        
    