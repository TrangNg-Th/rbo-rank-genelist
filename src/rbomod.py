
# -*- coding: utf-8 -*-
"""

@author: ttran
"""

from scipy.stats import poisson
from scipy.stats import skellam
from scipy.stats import triang
import numpy as np



# Implementation of rbo with a given function
def rbo_modified(
                lstS : list, 
                lstT : list,
                wg_func_name : str = 'poisson',
                params_wg_func : dict = None,
                k : int = 5, 
                strict_k : bool = False):
    
    """
    Parameters
    ----------
    wg_func_name : string  
        Name of the weighting function chosen
        Chosse between 
        ['geometric', ''poisson', 'skellam', 'triangular']
    

    params_wg_func : dict. (OPTIONAL)
                     Dictionary of the chosen weighting function and the chosen parameters
                     {'name': name_of_weighting_function, 
                     Chosse between 
                     ['geometric', ''poisson', 'skellam', 'triangular']
                      'param1' : parameter 1 associated with the weighting function,.... }
         
         
    lstS : list
        Ranked list S
    
    lstT : list
        Ranked list T
    
    k : TYPE, optional
        The number of elements that the user want to evaluate.
        The default is 5.
    
    strict_k : bool, optional
        If True, then only cares about the top k, not the tail.
        The default is False.

    Returns
    -------
    The RSM score for the given weighting function

    """

    
    S, T = lstS.copy(), lstT.copy()
    sS = len(S)  
    sT = len(T)
    
    # Set the weighting function
    name = wg_func_name
    assert name in ['geometric', 'poisson', 'skellam', 'triangular'], \
        "Please choose wegithing scheme in ['geometric', 'poisson', 'skellam', \
            'triangular']"
    if name == 'poisson':
        wg_func = wg_poisson
    elif name == 'skellam':
        wg_func = wg_skellam
    elif name == 'geometric':
        wg_func = wg_geom
    elif name == 'triangular':
        wg_func = wg_triangular        
        
    
    
    # if the user only cares about the top k
    if strict_k == True: 
        n = min(k, sS, sT) 
        A, AO = [0]*n, [0]*n  # initialize list of common elements

        
    elif strict_k == False:
        n = min(sS, sT) 
        A, AO = [0]*n, [0]*n       
    
    
    
    # if list S, T contain any element
    if (sS == 0) & (sT == 0) : return 1
    if (sS == 0) | (sT == 0) : return 0  # Either of the list is empty
    
    
    # Set up weighting parameters 
    default = {'geometric': {'p': 0.5}, 
               'poisson'  : {'lambda': k},
               'skellam'  : {'mu1': 2*k, 'mu2': k},
               'triangular': {'topk': k}}
    
    
    
    if params_wg_func is None :
        params_wg_func = default[name]
    
    
      
    
    # create weight vectors
    weights = [wg_func(params_wg_func, d) for d in range(n)]
    wgnormed = [i/ sum(weights) for i in weights]
    
    
    
    # case of first element in list
    A[0] = 1.0 if S[0] == T[0] else 0
    AO[0] = weights[0] if S[0] == T[0] else 0

    
    for d in range(1, n):  # Go through each depth to level d
        
        # define the intersection of subset of S[:d] w/ the subset of T[:d]
        Id = set(S[:d+1]).intersection(set(T[:d+1]))
        
        # define the number of intersected elements
        Xd = len(Id)
        
        # update agreement array
        A[d] = 1.0 * (Xd / (d + 1))
        
        
        # update the average overlap array
        AO[d] = AO[d - 1] + wgnormed[d] * A[d]
    
    
    
    return(np.round(AO[-1], 3))
    
    

# Poisson distribution
def wg_poisson(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : {'lambda' : the mean of the Poisson distribution}
        d = current rank of evaluation
        
    """
    
    lamb = params_wg_func['lambda']
    res = poisson.pmf(d, lamb)
    
    return(res)


# Geometrical progression
def wg_geom(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : {'p' : in [0, 1] parameter in Binomial distribution}
        d = current rank of evaluation
        
    """
    
    p = params_wg_func['p']
    
    # make sure that the value of p is correct
    assert (p > 0.0) & (p < 1.0), "Please input value of p in ]0, 1["
    
    
    res = 1.0 * (1 - p) * p**d
    return(res)



# Skellam distribution

def wg_skellam(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : {mu1: mean of the first Poisson, > 0
                          mu2 : mean of the second Poisson > 0 }
        d = current rank of evaluation
        
    """
   
    mu1 = params_wg_func['mu1']
    mu2 = params_wg_func['mu2']
    res = skellam.pmf(d, mu1, mu2)
    
    return(res)



# Triangular distribution

def wg_triangular(params_wg_func : dict, d):
    
    """
    Arguments:
        params_wg_func : {'topk' : the top k most important element}
        d = current rank of evaluation
        
    """
    
    topk = params_wg_func['topk']
    c = 1.0
    scale = topk
    loc = 0
    res = triang.pdf(d, c=c, scale=scale, loc=loc)
    
    return(res)


# =============================================================================
# # Implementation of original rbo
# def rbo(
#         S: list, 
#         T: list,
#         p=1.0, 
#         k=None, 
#         verbose=False):
#     """
#         Arguments :
#         
#             k    -- Depth of evaluation
#             p    -- probability that user will move forward to rank k
#         
#             setS -- Ranked List 1
#             setT -- Ranked List 2
#             d    -- Depth into the comparision
#         
#         Implementation of Eq.4
#         
#         Assume:
#             -> setS & setT can be of different length
#             -> Each element in l1, l2 occurs at most once in each list -->Need a function to evaluate this
#         
#         Return:
#         
#             RBO at depth k
#             
#     """
#     
#     # evaluate length of each list
#     sS = len(S)  
#     sT = len(T)
#     
#     
#     
#     # if rank k isn't given
#     # return the minimum between the length of each list and k
#     if k is None:
#         k = float("inf")
#     k = min(k, sS, sT)
#     
#        
#        
#     # initialize list of common elements at each rank
#     A = [0]*k  # proportion of agreement up to rank k 
#     AO = [0]*k  # average overlap of rank k
#     
#     
#     # if list S or T contains no element
#     if (sS == 0) & (sT == 0) : return 1  # both lists are empty
#     if (sS == 0) | (sT == 0) : return 0  # either of the list is empty
#     
#     
#     
#     # make sure that the value of p is correct
#     assert (p > 0.0) & (p < 1.0), "Please input value of p in ]0, 1["
#     
#     # define weight vectors
#     weights = [1.0 * (1 - p) * p**d for d in range(k)]
# 
#     # case of 1st element in list
#     A[0] = 1.0 if S[0] == T[0] else 0
#     AO[0] = weights[0] if S[0] == T[0] else 0
#     
#     
#     
#     for d in range(1, k):  # Go through each depth to level d
#         
#         # define the intersection of subset of S[:d] w/ the subset of T[:d]
#         Id = set(S[:d+1]).intersection(set(T[:d+1]))
#         
#         # define the number of intersected elements
#         Xd = len(Id)
#         
#         # update agreement array
#         A[d] = 1.0 * (Xd / (d + 1))
#         
#         
#         # update the average overlap array
#         if p == 1.0:
#             AO[d] = ((AO[d-1] * d) + A[d]) / (d + 1)
#         
#         else:  
#             AO[d] = AO[d - 1] + weights[d] * A[d]
#             
#             
#         # display the comparision between 2 lists at each position
#         if verbose == True:
#             print(f'Intersection set = {Id}')
#             
#         
#     # bound the result between 0.0 & 1.0
#     #if AO[-1] > 1.0 : 
#     #    AO[-1] = 1.0
#     #elif AO[-1] < 0.0:
#     #    AO[-1] = 0.0
#         
#     # display the result
#     if verbose == True:
#         print('-'*50)
#         print(f'Proportion of agreement between two lists at each position is {A}')
#         print(f'Cumulative average overlap at each position is {AO}')
#         print('-'*50)
#     return(AO[-1])
#     
# 
# =============================================================================



    
    
