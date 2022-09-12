# -*- coding: utf-8 -*-
"""
Created on Sat Sep 10 08:58:37 2022

@author: ttran
"""




# Implementation of original rbo

def rbo(
        S: list, 
        T: list,
        p=1.0, 
        k=None, 
        verbose=True)
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
    
    
    # define weight vectors
    if p == 1.0: 
        weights = [1.0 for _ in range(k)]  
    else:
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
    if AO[-1] > 1.0 : 
        AO[-1] = 1.0
    elif AO[-1] < 0.0:
        AO[-1] = 0.0
        
    # display the result
    if verbose == True:
        print('-'*50)
        print(f'Proportion of agreement between two lists at each position is {A}')
        print(f'Average overlap at each position is {AO}')
        print('-'*50)
        
        
    return(AO[-1])
    



