# Implementation of rbo with a given function

def rbo_modified(
                wg_func,
                lstS: list, 
                lstT: list,
                k=None, 
                verbose=False):
    """
         Arguments :
        
            k          -- Depth of evaluation
            wg_func    -- a weighting function k
        
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
    #if p == 1.0: 
    #    weights = [1.0 for _ in range(k)]  
    #else:
    weights = [wg_func(d) for d in range(k)]
    
    
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
            
        
    # bound the result
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
    
    
    
