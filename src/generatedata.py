# MODULE TO PLOT OUT MARKOV CHAIN
# Simulate lists of elements to test
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle
from random import random, choices, uniform, seed, sample, randint
import pandas as pd
from rbo import rbo_modified as SM
from rbo import wg_geom as wgm, wg_binomial as wbi, wg_poisson as wpo, wg_skellam as wsk, wg_triangular as wtr
import os
import sys
import itertools



def transitionMatrix(
    
    p0: list = None,
    states = ['H', 'MH', 'M', 'LM', 'L'],
    matrix = None,
    d : dict = None,
    error_rate : float = None, 
    nb_replicates: int = 1,
    sd = 1):
    """
    

    Parameters
    ----------
    p0 : list, optional
        LIST OF GENE EXPRESSIONS RANGING FROM 1 TO 100. The default is None.
    
    states : TYPE, optional
        LIST OF STATES : HIGH, MEDIUM HIGH, MEDIUM, LOW MEDIUM, LOW.
        The default is ['H', 'MH', 'M', 'LM', 'L'].
    
    matrix : TYPE, optional
        MARKOV TRANSITION MATRIX. The default is None.
    
    d : dict, optional
        DICTIONARY OF EACH STATES WITH THEIR RANGE. The default is None.
    
    error_rate : float, optional
        IF NONE, USE THE DEFAULT TRANSITION MATRIX (HARD CODE IN THIS CODE). The default is None.
    
    nb_replicates : int, optional
        NUMBER OF REPLICATES THE USER WANTS TO PRODUCE. The default is 1.
    
    sd : TYPE, optional
        ARGUMENT FOR RANDOMIZATION. The default is 1.

    Returns
    -------
    pk : list
        LIST OF GENE EXPRESSIONS PRODUCED FROM pO * matrix
    pkprimelist : list
        LIST OF ALL REPLICATED GENE EXPRESSIONS

    """
    
    seed(sd)
    
    if p0.all() == None : return None
    
    
    if d == None:
        d = {'H': [80, 100], 'MH': [70, 79], 'M': [30, 69], 'LM': [20, 29], 'L': [0.0, 19]}
        
         
    if (np.all(matrix) == None):
        matrix = np.array([
            [0.1, 0.2, 0.1, 0.1, 0.5],
            [0.1, 0.3, 0.1, 0.05, 0.45],
            [0.25, 0.02, 0.05, 0.65, 0.03],
            [0.4, 0.2, 0.15, 0.15, 0.1],
            [0.35, 0.2, 0.15, 0.15, 0.15]])
    
    
    #print(matrix)
    tmat = pd.DataFrame(matrix, index=states, columns=states)  # transition matrix from condition 0 to Xk
    #print(tmat)
    # give pk for condition k
    pk, state = [] , []
    for i in range(len(p0)):
        if p0[i] >= 80 : 
            state_i = choices(states, weights=tmat[states[0]].values)[0]
           
        elif (p0[i] >= 70) & (p0[i] < 80):
            state_i = choices(states, weights=tmat[states[1]].values)[0]
            
        elif (p0[i] >= 30) & (p0[i] < 70):
            state_i = choices(states, weights=tmat[states[2]].values)[0]
           
        elif (p0[i] >= 20) & (p0[i] < 30):
            state_i = choices(states, weights=tmat[states[3]].values)[0]
            
        elif (p0[i] < 20) :
            state_i = choices(states, weights=tmat[states[4]].values)[0]
        pk.append(uniform(d[state_i][0], d[state_i][1]))
        state.append(state_i) 
    

    # --------------------------
    # give list of pkprime from pk   
    pkprimelist, stateprimelist = [], []
    
    
    
    # if noise was requested  
    ntmatlist = []
    for rep in range(nb_replicates):
        pkprime, stateprime = [], []
        
        if error_rate is None: 
            nmatrix = ToleranceTransMat()
        
        else:
            nmatrix = ToleranceTransMat(error_rate = error_rate)
        ntmat = pd.DataFrame(nmatrix, index=states, columns=states)  # noise transition matrix
        
        # Save the transition matrices
        ntmatlist.append(ntmat)
        for i in range(len(pk)):
            if pk[i] >= 80 : 
                state_i = choices(states, weights=ntmat[states[0]].values)[0]
                
            elif (pk[i] >= 70) & (pk[i] < 80):
                state_i = choices(states, weights=ntmat[states[1]].values)[0]
            
            elif (pk[i] >= 30) & (pk[i] < 70):
                state_i = choices(states, weights=ntmat[states[2]].values)[0]
            
            
            elif (pk[i] >= 20) & (pk[i] < 30):
                state_i = choices(states, weights=ntmat[states[3]].values)[0]
           
            elif (pk[i] < 20) :
                state_i = choices(states, weights=ntmat[states[4]].values)[0]
            
            # if there's a change in the state of the current gene
            if state_i != state[i]:
                pkprime.append(uniform(d[state_i][0], d[state_i][1]))
            else:
                pkprime.append(pk[i])
            stateprime.append(state_i)
        
        pkprimelist.append(pkprime)  # list of all the noisy replicates
        stateprimelist.append(stateprime)
        
    return(pk, state, tmat, pkprimelist, stateprimelist, ntmatlist)

 


# Print out histograms of gene expressions   
def histogram(p, title):
    
    # default params
    bins=10
    density=True
    edgecolor='white'
    
    plt.figure(figsize=(10, 3))
    N, bins, patches = plt.hist(p, 
                                bins=bins, 
                                density=density, 
                                edgecolor=edgecolor)
    
    
    cmap = plt.get_cmap('jet')
    low = cmap(0.1)
    lmed = cmap(0.25)
    medium = cmap(0.5)
    medh = cmap(0.75)
    high = cmap(0.95)


    for i in range(0, 2):
        patches[i].set_facecolor(low)
    for i in range(2 , 3):
        patches[i].set_facecolor(lmed)
    for i in range(3, 7):
        patches[i].set_facecolor(medium)
    for i in range(7, 8):
        patches[i].set_facecolor(medh)
    for i in range(8, 10):
        patches[i].set_facecolor(high)

    # create legend
    handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in [low, lmed, medium, medh, high]]
    labels = ["low", "lmedium", "medium", "mediumh", "high"]
    plt.legend(handles, labels)
    
    plt.title = (f'Histogram of gene expression label obtained from {title}')
    plt.xlabel("Gene expression level")
    plt.ylabel("Frequency")
    plt.xticks(fontsize=10)  
    plt.yticks(fontsize=10)

    plt.gca().spines["top"].set_visible(False)  
    plt.gca().spines["right"].set_visible(False)
    
    plt.show()
    







def printout(condition: str, p, n):
    """ Function to print out condition and the associated vector """
    print()
    print(condition)
    print('-'*30)
    print(f'Propotion of state H: {sum(p >= 80) / n}')
    print(f'Propotion of state MH: {sum(((p >= 71) & (p <= 79))) / n}')
    print(f'Propotion of state M: {sum(((p >= 30) & (p <= 70))) / n}')
    print(f'Propotion of state LM: {sum(((p >= 21) & (p <= 29))) / n}')
    print(f'Propotion of state L: {sum((p < 20)) / n}')
    print()


    
    

########### Generate data
def generatelist(
    n: int = 1000,
    Xk: str = 'Condition_1',
    mk = None,  # transition matrix to go from condition 0 to condition k
    sd = 1, 
    nb_replicates: int = 1,
    conditionX0: list = None):
    """
    Design : 
    ---------
    We are interested in comparing two lists of genes in 'comparable contexts' (or comparable conditions).
    This function generates lists of modified conditioned based from one initial condition.
    
        Condition 0 : 
            The condition 0 = condition when we found these genes in the nature
            We simulate a set of n expression values p1,..., pn of n genes g1, g2, ... gn
        
                Denote p = (p1, ..., pn)
        
                80 <= pj <= 100 : higly expressed --> state : H
                70 <= pj <= 80 :                  --> state : unk1
                30 <= pj <= 70 : mediumly expressed --> state : M
                20 <= pj <= 30                      --> state unk2
                0 <= pj <= 20 : lowly expressed    --> state : L
                
        
        Denote Xk - an xperimental condition.
        Condition Xk creates a different condition where one / multiple genes can change their expression value
            Each gene has a probability to switch from one state to another state by Markov model
            After putting our subjects through Xk, we retrieve a vector of gene expression values.
            Let's pk = the gene expression vector
            
            

    Parameters
    ----------
    n : int, optional
        NUMBER OF GENES IN THE SAMPLE. The default is 1000.
        
    Xk : str, optional
        NAME OF THE EXPERIMENT. The default is 'Condition_1'.
        
    mk : TYPE, optional
        TRANSITION MATRIX TO GO FROM CONDTION 0 TO CONDITION K. The default is None.
        
    sd : TYPE, optional
        ARGUMENT FOR REPRODUCIBILITY OF RANDOMIZATION. The default is 1.
        
    nb_replicates : int, optional
        NUMBER OF REPLICATES FOR THE CURRENT EXPERIMENT. The default is 1.
        
    conditionX0 : list, optional
        IF NOT NONE, IT SHOULD BE A LIST OF GENE EXPRESSIONS OF SIZE n
        IF NONE, USE THE DEFAULT LIST. The default is None.

    Returns
    -------
    p0 : LIST OF GENE EXPRESSIONS OF conditionX0
    dk : DATAFRAME OF GENE EXPRESSIONS OF conditionXk, conditionX0
    """
    
    seed(sd)
    dk = pd.DataFrame(index=[f'gene_{i}' for i in range(n)])
    
    # PART 1 :    
    # generate baseline condition (if not given)
    if conditionX0 == None:
                
        #print("Using default condition 0")
        conditionX0 = "Condition_0"
        
        pH = np.array(choices(range(80, 101), k = int(0.1 * n)))
        pM = np.array(choices(range(30, 70), k = int(0.8 * n)))
        pL = np.array(choices(range(0, 20), k = int(0.1 * n)))
        p0 = np.concatenate((pH, pM, pL), axis=None)   
        
    else : 
        assert conditionX0 == None, 'Currently, this data generator uses a default condition_0'
   
    
    state0 = [] # Add the column of states for this baseline condition

    for k in range(len(p0)):
        if p0[k] >= 80:
            state0.append('H')
            
        elif (p0[k] >= 70) & (p0[k] < 80):
            state0.append('MH')
            
        elif (p0[k] >= 30) & (p0[k] < 70):
            state0.append('M')
        
        elif (p0[k] >= 20) & (p0[k] < 30):
            state0.append('LM')
            
        elif (p0[k] < 20):
            state0.append('L')
            
    
    # Display info
    #print('_'*50)
    #printout(conditionX0, p0, n)
    #histogram(p0, title = conditionX0)
    dk[conditionX0] = p0
    dk[conditionX0 + '_state'] = state0   
    
    
  
    
    # PART 2:     
    # for the experiment Xk, return the result vector of the experiment
    # simulate results from each condition
    
    
    # transition matrix base from condition 0
    pk, statek, tmatk, pkprimelist, statekplist, ntmatklist = transitionMatrix(\
                p0=p0, matrix=mk, sd=seed, nb_replicates=nb_replicates) 
    dk[Xk] = pk  # memorize gene expressions
    dk[Xk+'_state'] = statek  # memorize states
    
    # Displaying results
    #printout(Xk, dk[Xk], n)
    #histogram(pk, title=Xk)
    
    # Save the transition matrices in a dictionary
    d = {}
    d[Xk] = tmatk
    
    # For the Xkprime
    for i in range(len(statekplist)):
        dk[Xk+'_prime_'+str(i+1)] = pkprimelist[i]
        dk[Xk+'_prime_state_'+str(i+1)] = statekplist[i]
        #printout(Xk+'_prime_'+str(i+1), dk[Xk+'_prime_'+str(i+1)], n)
        d[Xk+'_prime_'+str(i+1)] = ntmatklist[i]
    return(dk, d)




## Given an error rate,  return transition matrix
def ToleranceTransMat(
    error_rate : float = 0.05,
    states : list = ['H', 'MH', 'M', 'LM', 'L']):
    """
    Design
    ----------
    In this simulation Xkprime = transition_matrix*Xk
    where transition_matrix = Markov transition matrix
    
    
    Parameters
    ----------
    error_rate : float, optional
        IN Xkprime THERE IS 1 - ERROR_RATE PERCENTAGE OF GENES THAT REMAIN WITHIN THE STATE.
        The default is 0.05.
    
    states : list, optional
        LIST OF STATES : HIGH, MEDIUM HIGH, MEDIUM, LOW MEDIUM, LOW.
        The default is ['H', 'MH', 'M', 'LM', 'L'].

    Returns
    -------
    matrix : 
        MARKOV TRANSITION MATRIX

    """
    
   

    
    matrix = np.empty([len(states), len(states)])
    for j in range(len(states)):
        for i in range(len(states)):
            if i != j:  # if the current position is not on the diagonal of the matrix
                matrix[i,j] = error_rate / (len(states) - 1)  # divide evenly the error_rate to this position
            else:
                matrix[i, j] = 1 - error_rate  # if we're on the diagonal, we'll have the biggest possibility
    return(matrix)
            
    
    
# Compare generated data between groups
def comparegeneratedlist(dataframe,
                         group0:str = 'Condition_0_state',
                         group:str = 'Condition_1_state'):
    """
    Description :
    ----------
    GIVEN A DATAFRAME OF DIFFERENT GENE EXPRESSIONS,
    RETURN THE AVERAGE NUMBER OF GENES THAT SWITCHED THE CATEGORY

    Parameters
    ----------
    dataframe : TYPE
        DATAFRAME OF DIFFERENT GENE EXPRESSIONS.
        
    group0 : str, optional
        STATES OF A GIVEN CONDITION. The default is 'Condition_0_state'.
    
    group : str, optional
        STATES OF A GIVEN CONDITION. The default is 'Condition_1_state'.

    Returns
    -------
    avg : 
        COUNT NUMBER OF ELEMENTS THAT CHANGED THE CATEGORY

    """
    
    similarity = {}
    for i in dataframe[group].value_counts().index:
        if i in dataframe[group0].value_counts().index:
            res = (dataframe[group0].value_counts()[i] - dataframe[group].value_counts()[i])**2
        else:
            res = dataframe[group].value_counts()[i]
        similarity[i] = res
    
    avg = np.sqrt(sum(similarity.values()))/5
    return(avg)




# Function to generate a dataframe

# Call function to generate transition matrices for each Xk

def generatedatadefault(nbXk, nbrep, N, sd):
    """
    

    Parameters
    ----------
    nbXk : int
        NUMBER OF EXPERIMENTS X
    nbrep : int
        NUMBER OF REPLICATES FOR EACH EXPERIMENTS
    N : int
        SAMPLE SIZE, NUMBER OF INDIVIDUALS
    sd : int
        RANDOMIZATION PARAMETER. SEED

    Returns
    -------
    dataframe:
        DATAFRAME OF ALL GENERATED DATA

    """
    
    list_Xk = []  # list of transition matrices for k experiments X
    
    for n in range(nbXk):
        seed(n)
        res = []
        
        for i in range(5):
            r = [random() for _ in range(5)]
            s = sum(r)
            res.append([k/s for k in r])
            restmp = np.array(res)
        list_Xk.append(restmp)
    
    
    
    print('Generating data....')
    print(f'Sample size : {N}')
    print(f'Number of experiments: {nbXk}')
    print()
    dataframe, transmatrix = generatelist(N, Xk='Condition_1', 
                                          nb_replicates=nbrep, mk=list_Xk[0]) 

    sd = 2
    for mat in list_Xk[1:]:
        name = 'Condition_' + str(sd)
        dk, d = generatelist(N, Xk=name, sd=sd, nb_replicates=nbrep, mk=mat)
        dataframe = pd.concat([dataframe, dk[dk.columns[2:]]], axis=1)
        transmatrix.update(d)
        sd += 1

    # Export data generated
    source = os.path.dirname(os.getcwd())
    path = source + '/data/'
    name = "testdata.csv"
    dataframe.to_csv(path+name)
    
    print(f'Done. Check data in {path}\n Filename: {name}')
    return(dataframe)



def similarityscorealldata(nbXk, nbrep, dataframe, condition, weightfunc, weightparams, namedf = None):
    """
    

    Parameters
    ----------
    nbXk : int
        NUMBER OF EXPERIMENTS X
    nbrep : int
        NUMBER OF REPLICATES FOR EACH EXPERIMENTS
    dataframe:
        DATAFRAME OF ALL GENERATED DATA
    condition : TYPE
        CONDITION WHICH FILTERS OUT THE GENES WHICH BELONG TO A SPECIFIC CATEGORY
    weightfunc : TYPE
        CHOSEN WEIGHTING SCHEMES.
    weightparams : TYPE
        DICTIONARY OF PARAMETERS ASSOCIATED WITH WEIGHTING SCHEMES
    namedf : (optional)
        NAME OF DATAFRAME TO BE GENERATED

    Returns
    -------
    df : DATAFRAME OF PAIRWISE SIMILARITY SCORES

    """
   
    # Create dataframe of pairwiser similarity score for each experiment & replicates
   
    collist = []
    for k in range(1, nbXk+1):
        for j in range(0, nbrep+1):
            if j == 0:
                collist.append(f'Condition_{k}')
            else:
                collist.append(f'Condition_{k}_prime_{j}')
    
    # Create a matrix to compare rbo scores
    df = pd.DataFrame(index=collist, columns=collist)
    
    
    lindex = list(df.index)
    for i in lindex:
        name_i = i
        
        if 'prime' in i:
            state_i = '_'.join(i.split('_')[:-1]) + '_state_' + i.split('_')[-1]
            
        else:
            state_i = i + '_state'
        
        #print('name_i', name_i)
            
        for j in lindex:
            name_j = j
            if 'prime' in j:
                state_j = '_'.join(j.split('_')[:-1]) + '_state_' + j.split('_')[-1]
            else:
                state_j = j + '_state'
            
            c1 = (dataframe[state_i] == condition)   
            l1 = list(dataframe[c1].sort_values(name_i, ascending=False).index)
           
            c2 = (dataframe[state_j] == condition)  
            l2 = list(dataframe[c2].sort_values(name_j, ascending=False).index)

            
            s = SM(weightfunc, weightparams, l1, l2)
            df.loc[name_i, name_j] = np.round(s, 4)
            
    # Export data generated
    if namedf == None:
        name = "SM_data.csv"
    else:
        name = namedf
        
    source = os.path.dirname(os.getcwd())
    path = source + '\\data\\'
    df.to_csv(path+name)
    df = pd.read_csv(path+name)
    
    print(f'Done. Check data in {path}\n Filename: {name}')
    return(df)
           
    
    
    
# =============================================================================
# def PermutationWORepl(l):
#     """
#     
# 
#     Parameters
#     ----------
#     l : list
#         LIST CONTAINS UNIQUE ELEMENTS
# 
#     Returns
#     -------
#     L_p : list
#         LIST OF ALL PERMUTATIONS OF THIS LIST
# 
#     """
#     
#     n = len(l)
#     T_p = list(itertools.permutations(l, n))
#     L_p = [list(i) for i in T_p]
#     return(L_p)
#     
# =============================================================================

# =============================================================================
# 
# 
# def Smallpermutationwithrepl(size_l: int = 5, 
#                              output_path: str = '/home/nguyetrt/rbo-rank-genelist/data/', 
#                              nb_repl: int = 0,
#                              l: list = None,
#                              weightfunc: list = None, 
#                              weightparams: list = None, 
#                              nb_randomization: int = 100
#                              ):
#     """
#     
# 
#     Parameters
#     ----------
#     size_l : INT
#         SIZE OF THE INPUT LIST. GIVE l IS GIVEN, THIS ARG WILL BE DISCARDED
#         
#     output_path : str
#         PATH WHERE TO SAVE THE GENERATED FILES.
#         
#     nb_repl: int, optional
#         NUMBER OF ELEMENTS TO BE REPLACED. the default is 0
#         
#     l : list, optional
#         LIST OF ELEMENTS TO BE PERMUTED. The default is None.
#         
#    
#     weightfunc : list, optional
#         FUNCTION WHICH PRODUCE APPROPRIATE WEIGHTING SCHEMES. The default is None.
#         
#     weightparams : list, optional
#         DICTIONARY OF PARAMETERS ASSOCIATED WITH WEIGHTING SCHEMES. The default is None.
#         
#     nb_randomization : int, optional
#         NUMBER OF LIST TO BE RETURNED
#     
# 
#     Returns
#     -------
#     None.
# 
#     """
#    
#     # if l is None, generate a list of [1,...,size_l]
#     # l = [0,..., size_l]
#     # l_ref = [0,..., size_l - nb_repl]
#     if l == None:
#         l_ref = list(np.arange(0, size_l-nb_repl, 1))  # generate a random list of size (size_l - nb_repl)
#         l = list(np.arange(0, size_l, 1))
#     
#     else:
#         assert l is None, "This function will generate a list"
#     
#     
#     
#     max_l = max(l)  # set for the biggest number
# 
#     # verify if there is any replacement element
#     print()
#     print(f'Build similarity scores distribution for lists of length {size_l} and {nb_repl} replacements')
#     
#     # in case there's replacement, l
#     while nb_repl > 0:
#         r = randint(max_l+1, max_l+1000)
#         nb_repl -= 1 
#         l_ref += [r]
# 
# 
#     l_permutations = [l_ref]
# 
#     # shuffle the list l, to return a list of random l
#     while nb_randomization - 1> 0:
#         l_shuffled = sample(l_ref, k=len(l_ref))
#         l_permutations.append(l_shuffled)
#         nb_randomization -= 1
#     
#     
#     
#     if (weightfunc is None) | (weightparams is None):
#         print('No weighting scheme chosen, will use the triangular weighting scheme')
#         weightparams = {'triangular': None, 'topk': size_l}
#         weightfunc = wtr
#     
#     
#     # List of all similarity scores for all permutations       
#     j, count = 0, 1
#     div = np.floor(len(l_permutations)/100)
#     
#     l_scores = []
#     
#     for i in range(len(l_permutations)): 
#         l_shuffle = l_permutations[i]
#         e = np.round(SM(weightfunc, weightparams, l, l_shuffle), 3)
#         l_scores.append(e)
#         
#         if (i >= j) :
#             sys.stdout.write('\r')
#             sys.stdout.write("[%-1s] %d%%" % ('='*count, count))
#             sys.stdout.flush()
#             j += div
#             count +=1
# 
#     l_permutations.append(l)
#     l_scores.append(np.round(SM(weightfunc, weightparams, l, l), 3))
#     l_scores_str = str(l_scores)
#     l_scores_str = l_scores_str[1:-1]
#     f = open(output_path + f'perm_list_len_{size_l}.txt',  'w')
#     f.write(l_scores_str)
#     f.close()
#     
#     # Creating histogram
#     #fig, axs = plt.subplots(1, 1, figsize =(10, 7))
#     #axs.hist(l_scores)
#     #plt.suptitle(f'Distribution of similarity scores for lists of length {size_l}')
#     # Show plot
#     #plt.show()
#     print()
#     print(f'Done generating scores, check file /data/perm_list_len_{size_l}.txt for result')
#     return(l_permutations, l_scores)
# 
#     
#         
# =============================================================================

  



def simulation(size_l: int = 10, 
               withreplacement: bool = False,
               output_path: str = '/home/nguyetrt/rbo-rank-genelist/data/', 
               weightfunc: list = None, 
               weightparams: list = None, 
               nb_generatedlist: int = 50,
               size_vocabulary: int = 100,
               Nb_r: int = 1):
    """
    

    Parameters
    ----------
    size_l : int, optional
        The size of the list to be generated. The default is 20.
    withreplacement: bool, optional
        Parameter that controls if the user wants to include replacement of elements in the the simulation
    output_path : str, optional
        The directory to save generated file. The default is '/home/nguyetrt/rbo-rank-genelist/data/'.
    sd : int, optional
        An integer to fix the randomization. The default is 1.
    weightfunc : list, optional
        Weighting function to be applied. The default is None.
    weightparams : list, optional
        Dictionary of parameters associated with the chosen weighting scheme. The default is None.
    nb_generatedlist : int, optional
        Number of desired elements to be generated by this function. The default is 100.
    size_vocabulary: int, optional
        The size of the characters for the list l to be generated from. The default is 100.
    Nb_r: int, optional
        Number of elements that will be replaced if the "withreplacement" parameter is chosen.
        The default is 1

    Returns
    -------
    None.

    """
    
    n = nb_generatedlist
    V = np.arange(0, size_vocabulary, 1)
    
    # part 1: simulate a list of rankings
    
    if withreplacement == False:  # if there's no replacement to be considered in the model
        l = list(np.arange(0, size_l, 1))
        
        # result 
        output = [l]
        
        while n > 0:
            l_shuffle = sample(l, k=len(l))  # generate a shuffled version of l
            output.append(l_shuffle)
            n -= 1 
        
    
    else:  # if we do consider the replacement, then 
        output = []
        while n > 0:
            #seed(n)
            # 1. randomly generate a number Nb_r : this is the number of replacement in the ranking
            Nb_r = randint(0, size_l)
            print('Number of replaced items=', Nb_r)

            # l = [0, 1, ...., size_l - Nb_r]
            l = list(np.arange(0, size_l - Nb_r))   
            
            
            # 2. Complete the rest by generating random number not in l
            for i in range(size_l - Nb_r, size_l):
                #seed(n + i)
                r = randint(V[i], V[-1]) # generate a number to fill in l
                l.append(r)
                
            # 3 . Additional shuffling of the list
            l = sample(l, k=len(l))  
            
            output.append(l)
            n -= 1        

    # Part 2: Similarity ranking measure evaluation
    
    if (weightfunc is None) | (weightparams is None):
        print('No weighting scheme chosen, will use the triangular weighting scheme')
        weightparams = {'triangular': None, 'topk': size_l}
        weightfunc = wtr
    
    
    
    output_scores = []
    
    for i in range(len(output)): 
        for j in range(i+1, len(output)):
            l_i = output[i]
            l_j = output[j]
            e = np.round(SM(weightfunc, weightparams, l_i, l_j), 3)
            output_scores.append(e)
        
        
    
    #output_scores.append(np.round(SM(weightfunc, weightparams, l, l), 3))
    l_scores_str = str(output_scores)
    l_scores_str = l_scores_str[1:-1]
    
    #f = open(output_path + f'perm_list_len_{size_l}.txt',  'w')
    #f.write(l_scores_str)
    #f.close()
    
    # Creating histogram
    fig, axs = plt.subplots(1, 1, figsize =(10, 7))
    axs.hist(output_scores)
    plt.suptitle(f'Distribution of similarity scores for lists of length {size_l} with repl{withreplacement}')
    # Show plot
    plt.show()
    #print()
    #print(f'Done generating scores, check file /data/perm_list_len_{size_l}.txt for result')
    return(output, output_scores)

    

    

        
        
        