# -*- coding: utf-8 -*-
import numpy as np


from scipy.optimize import root
import matplotlib.pyplot as plt

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
    
    assert n > m, "The sample size is not big enough, please enter a bigger sample"
    assert n % m == 0, "Please enter a multiple of 10"
    
    Z = []
    
    inc = n// m
    for i in range(m-1):
        z = max(X[i*inc : (i+1) * inc])
        Z.append(z)
    return(Z)
    



# Paper of Theorie des valeurs extremes univariées
# p. 28


def plus(x):
    return(np.max([x, 0], 0))


def maxlogvraissemblance(gamma, a_m, b_m):
    
    Z = maximaperbloc(X)
    
    if gamma < -1:
        print('There is no consistent estimator')
        return(None)
    
    
    elif gamma == 0:
      
        term_1 = m * np.log(a_m)
        term_2 = sum([ np.exp( -( z - b_m ) / a_m) for z in Z ])
        term_3 = sum([ - (z - b_m ) / a_m for z in Z ])    
        
        res = - term_1 - term_2 - term_3
        return(res)
    
    else:
        g = gamma
        
        term_1 = m * np.log(a_m)
        term_2 = (1 + 1 / g) * sum([ np.log(1 + g * plus(( z - b_m ) / a_m)) for z in Z ])
        term_3 = sum([ 1 + g * plus((z - b_m ) / a_m)**(1 / g) for z in Z ])    
        
        res = - term_1 - term_2 - term_3
        return(res)
    
    
    
    
# Solve for solution of logmaximumvraisemblance
def sovlemaxlogvraisemblance(fun, x0, method = 'hybr'):
    sol = root(fun, x0, method=method)
    return(sol)


# Equation 1.32
def quantileExtreme(a_m, b_m, g_m, p_m, m):
    
    
    
    q_m = b_m + (a_m / g_m) * (1 / (m * p_m)**g_m - 1)
    return(q_m)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    