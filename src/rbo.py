
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 08:58:37 2022

@author: ttran
"""

from scipy.stats import poisson
from scipy.special import binom
import matplotlib.pyplot as plt

# Implementation of original rbo

def rbo(
        S: list, 
        T: list,
        p=1.0, 
        k=None, 
        verbose=False):
    """
        Arguments :
        
            k    -- Depth of evaluation
            p    -- probability that user will move forward to rank k
        
            setS -- Ranked List 1
            setT -- Ranked List 2
            d    -- Depth into the comparision
        
        Implementation of Eq.4
        
        Assume:
            -> setS & setT can be of different length
            -> Each element in l1, l2 occurs at most once in each list -->Need a function to evaluate this
        
        Return:
        
            RBO at depth k
            
    """
    
    # evaluate length of each list
    sS = len(S)  
    sT = len(T)
    
    
    
    # if rank k isn't given
    # return the minimum between the length of each list and k
    if k is None:
        k = float("inf")
    k = min(k, sS, sT)
    
       
       
    # initialize list of common elements at each rank
    A = [0]*k  # proportion of agreement up to rank k 
    AO = [0]*k  # average overlap of rank k
    
    
    # if list S or T contains no element
    if (sS == 0) & (sT == 0) : return 1  # both lists are empty
    if (sS == 0) | (sT == 0) : return 0  # either of the list is empty
    
    
    
    # make sure that the value of p is correct
    assert (p > 0.0) & (p < 1.0), "Please input value of p in ]0, 1["
    
    # define weight vectors
    weights = [1.0 * (1 - p) * p**d for d in range(k)]

    # case of 1st element in list
    A[0] = 1.0 if S[0] == T[0] else 0
    AO[0] = weights[0] if S[0] == T[0] else 0
    
    
    
    for d in range(1, k):  # Go through each depth to level d
        
        # define the intersection of subset of S[:d] w/ the subset of T[:d]
        Id = set(S[:d+1]).intersection(set(T[:d+1]))
        
        # define the number of intersected elements
        Xd = len(Id)
        
        # update agreement array
        A[d] = 1.0 * (Xd / (d + 1))
        
        
        # update the average overlap array
        if p == 1.0:
            AO[d] = ((AO[d-1] * d) + A[d]) / (d + 1)
        
        else:  
            AO[d] = AO[d - 1] + weights[d] * A[d]
            
            
        # display the comparision between 2 lists at each position
        if verbose == True:
            print(f'Intersection set = {Id}')
            
        
    # bound the result between 0.0 & 1.0
    #if AO[-1] > 1.0 : 
    #    AO[-1] = 1.0
    #elif AO[-1] < 0.0:
    #    AO[-1] = 0.0
        
    # display the result
    if verbose == True:
        print('-'*50)
        print(f'Proportion of agreement between two lists at each position is {A}')
        print(f'Cumulative average overlap at each position is {AO}')
        print('-'*50)
    return(AO[-1])
    




# Implementation of rbo with a given function

def rbo_modified(
                wg_func,
                params_wg_func : dict,
                lstS: list, 
                lstT: list,
                k = None, 
                verbose = False):
    """
         Arguments :
        
            k          -- Depth of evaluation
            wg_func    -- a weighting function k
            params_wg_func -- {'weight-function-name' : None, 
                                'parameter_1' : , 
                                'parameter_2' : }
        
            setS -- Ranked List 1
            setT -- Ranked List 2
               
        
        Implementation of Eq.4
        
        Assume:
            -> setS & setT can be of different length
            -> Each element in l1, l2 occurs at most once in each list -->Need a function to evaluate this
        
        Return:
        
            RBO at depth k
            
            
    """
    
    
    S, T = lstS.copy(), lstT.copy()
    sS = len(S)  
    sT = len(T)
    
    
    
    # if rank k isn't given
    if k is None:
        k = float("inf")
    k = min(k, sS, sT)   
       
    
    # initialize list of common elements at each rank
    A, AO = [0]*k, [0]*k
    
    
    # if list S, T contain any element
    if (sS == 0) & (sT == 0) : return 1
    if (sS == 0) | (sT == 0) : return 0  # Either of the list is empty
    
    
    # create weight vectors
    if params_wg_func == {}: 
        weights = [1.0 for _ in range(k)]  
    
    elif 'poisson' in params_wg_func:
        weights = [wg_func(params_wg_func, d) for d in range(k)]
        # generate Poisson distribution with sample size 10000
        #t = range(k)
        #plt.plot(t, weights, ':')
        #plt.show()
        
    elif 'geometric' in params_wg_func:
        weights = [wg_func(params_wg_func, d) for d in range(k)]
        
    elif 'hypergeometric' in params_wg_func:
        weights = [wg_hypergeom(params_wg_func, d) for d in range(k)]
        #print(weights)
        #t = range(k)
        #plt.plot(t, weights, ':')
        #plt.show()
        
    elif 'binomial' in params_wg_func:
        weights = [wg_binomial(params_wg_func, d) for d in range(k)]
        #print(weights)
        #t = range(k)
        #plt.plot(t, weights, ':')
        #plt.show()
    
    # case of first element in list
    A[0] = 1.0 if S[0] == T[0] else 0
    AO[0] = weights[0] if S[0] == T[0] else 0

    

    for d in range(1, k):  # Go through each depth to level d
        
        # define the intersection of subset of S[:d] w/ the subset of T[:d]
        Id = set(S[:d+1]).intersection(set(T[:d+1]))
        
        # define the number of intersected elements
        Xd = len(Id)
        
        # update agreement array
        A[d] = 1.0 * (Xd / (d + 1))
        
        
        # update the average overlap array
        AO[d] = AO[d - 1] + weights[d] * A[d]
            
            
        # display the comparision between 2 lists
        if verbose == True:
            print(f'Intersection set ={Id}')
        
    # display the result
    if verbose == True:
        print('-'*50)
        print(f'Proportion of agreement between two lists at each position is {A}')
        print(f'Cumulative average overlap at each position is {AO}')
        print('-'*50)

    return(AO[-1])
    
    

# Poisson distribution
def wg_poisson(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : [lambda]
        d = current rank of evaluation
        
    """
    
    lamb = params_wg_func['lambda']
    res = poisson.pmf(d, lamb)
    
    return(res)


# Geometrical progression
def wg_geom(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : [p]
        d = current rank of evaluation
        
    """
    
    p = params_wg_func['p']
    
    # make sure that the value of p is correct
    assert (p > 0.0) & (p < 1.0), "Please input value of p in ]0, 1["
    
    
    res = 1.0 * (1 - p) * p**d
    return(res)


# Negative hypergeometric distribution
# https://en.wikipedia.org/wiki/hypergeometric_distribution
def wg_hypergeom(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : [N: total number of elements,
                          K : total number of overlap]
        d = current rank of evaluation
        
    """
    
    N = params_wg_func['N']
    K = params_wg_func['K']
    n = params_wg_func['n']
    res = binom(K, d)*binom(N-K, n-d)/binom(N, n)
    
    return(res)



# Binomial distribution

def wg_binomial(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : [n: nb of elements to be considered
                          p : probability of success]
        d = current rank of evaluation
        
    """
   
    n = params_wg_func['n']
    p = params_wg_func['p']
    res = binom(n, d) * (p**d) * (1-p)**(n-d)
    
    return(res)

    
    
