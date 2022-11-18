# -*- coding: utf-8 -*-
import numpy as np
import os


# =============================================================================
# from scipy.optimize import newton
# from random import random
# 
# from simRSM import simulation as simRSM
# =============================================================================


import warnings
warnings.filterwarnings("ignore")

def maximaperbloc(X: list, m = 10):
    """
    Parameters
    ----------
    X : list
        List of RBO scores
        
    m : int
        The 

    Returns
    -------
    Z : list
        List of maxima of m blocs of RBO scores, where Z_i = max(X_m*i,..., X_m*(i+1))
        

    """
    
    n = len(X)
    
    if n >= 1000:
        m = 100
    
    Z = []
    
    inc = n // m
    for i in range(m-1):
        z = max(X[i*inc : (i+1) * inc])
        Z.append(z)
    Z.append(max(X[(i+1) * inc:]))
    return(Z, m)
    



# Paper of Theorie des valeurs extremes univariées
# p. 28


def plus(x):
    return(np.max([x, 0], axis=0))

def maxlogvraissemblance(x):
    Z, m = maximaperbloc(X)
    
    gamma = x[0]
    a_m = x[1]
    b_m = x[2] 
    
    #print('gamma=',gamma)
    #print('a_m=', a_m)
    #print('b_m=', b_m)
    
    
    
    if gamma == 0:
      
        term_1 = m * np.log(a_m)
        term_2 = sum([ np.exp( - ( z - b_m ) / a_m) for z in Z ])
        term_3 = sum([( z - b_m ) / a_m for z in Z ])    
        
        res = - term_1 - term_2 - term_3
        return(res)
    
    else:
        
        g = gamma
        
        term_1 = m * np.log(a_m)
        term_2 = (1 + 1 / g) * sum([ np.log(1 + g * plus(( z - b_m ) / a_m)) for z in Z ])
        term_3 = sum([ plus(1 + g * ( z - b_m ) / a_m)**(-1 / g) for z in Z ])    
        #print("term_1=", term_1, "term_2=", term_2, "term_3=", term_3)
        res = - term_1 - term_2 - term_3
        
        return(res)
    
    
# Equation 1.32
def quantileExtreme(a_m, b_m, g_m, p_m, m):
    q_m = b_m + (a_m / g_m) * (1 / (m * p_m)**g_m - 1)
    q_m = min(1.0, q_m)
    return(q_m)
    
# Solve for solution of logmaximumvraisemblance
# =============================================================================
# =============================================================================

# Import RSM scores generated
# Export data generated
source = os.path.dirname(os.getcwd())
path = source + '/data/'
filename = path + 'RSM_sim.txt'
with open(filename) as file:
    X = [float(line.rstrip()) for line in file]
#_, X = simRSM()


# =============================================================================
# g_0 = 0
# a_0 = 0.1
# b_0 = 0.9
# g_m, a_m, b_m = newton(maxlogvraissemblance, [g_0, a_0, b_0], maxiter=300)
# c = 1
# while sum(np.isfinite([g_m, a_m, b_m])) < 3:
#     print('Initizalization: ', c,'-th')
#     g_m, a_m, b_m = newton(maxlogvraissemblance, [g_0, a_0, b_0], maxiter=300)
#     g_0 = random()
#     a_0 = random()
#     b_0 = random()
#     
#     c += 1
#     #print(g_0, a_0, b_0)
# 
#     
#     
# print('Result of approximation')
# print('gamma=', g_m)
# print('a_m=', a_m)
# print('b_m=', b_m)
# 
# if g_m == -1:
#     print('Warning, this estimation is not consistent')
# p_m = 0.01
# m = 10
# q_m = quantileExtreme(a_m, b_m, g_m, p_m, m)
# print('Estimated extreme quantile q=', q_m, 'for p=', p_m)
#     
# 
# =============================================================================
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    