#  -*- coding: utf-8 -*-

"""
Created on Tue Nov 15 21:31:30 2022

@author: ttran


This module will generate gene lists
"""

# Import packages
from random import seed
from random import choices
from random import uniform





import numpy as np
import pandas as pd

# Function which takes an error rate and return the transition matrix
def ToleranceTransMat(error_rate : float = 0.05,
                      states : list = ['H', 'MH', 'M', 'LM', 'L']):
    
    
    """
    Design
    ----------
    In this simulation, each ranking obtained from a replicate is supposed to
    be similar to the ranking from the experiment 
    
    The mathematic relationship between the gene expression list of replicate 
    Xkprime and the experiment Xk is :
        
    Xkprime = transition_matrix * Xk
    
    
    where transition_matrix is Markov transition probabilities of having a new 
    gene expression value
    
    
    Parameters
    ----------
    error_rate : float, optional
        
        The difference between Xk and Xkprime is =  error_rate.
        More precisely, 1 - error_rate % of genes remain where they are.
        The default value is 0.05.
    
    states : list, optional
    
        Possible categories of gene expression level :
            high, medium high, medium, low medium, low
        The default is ['H', 'MH', 'M', 'LM', 'L'].

    Returns
    -------
    matrix : 
        Markov transition matrix

    """
    

    
    matrix = np.empty([len(states), len(states)])
    
    for j in range(len(states)):
        
        for i in range(len(states)):
            
            if i != j:  # if we're not on the diagonal
            
                # matrix[i, j] = error_rate / number of cells - 1 
                matrix[i,j] = error_rate / (len(states) - 1)  
            else:
                
                # if we're on the diagonal
                matrix[i, j] = 1 - error_rate  
    
    return(matrix)
            



def transitionMatrix(
    
    p0 : list = None,
    states = ['H', 'MH', 'M', 'LM', 'L'],
    matrix = None,
    d : dict = None,
    error_rate : float = None, 
    nb_replicates : int = 1,
    sd = 1):
    
    """
    

    Parameters
    ----------
    p0 : list, optional
        List of gene expressions ranging from 1 to 100. the default is none.
    
    states : TYPE, optional
        list of states : high, medium high, medium, low medium, low.
        The default is ['H', 'MH', 'M', 'LM', 'L'].
    
    matrix : TYPE, optional
        Markov transition matrix. The default is none.
    
    d : dict, optional
        Dictionary of each states with their range. the default is none.
    
    error_rate : float, optional
        If none, use the default transition matrix (hard code in this code).
        The default is none.
    
    nb_replicates : int, optional
        Number of replicates the user wants to produce. 
        The default is 1.
    
    sd : TYPE, optional
        Argument for randomization. The default is 1.

    Returns
    -------
    pk : list
        List of gene expressions produced from pO * matrix
    pkprimelist : list
        List of all replicated gene expressions

    """
    
    seed(sd)
    
    if p0.all() == None : return None
    
    
    if d == None:
        d = {'H': [80, 100], 'MH': [70, 79], 'M': [30, 69], 'LM': [20, 29], 
             'L': [0.0, 19]}
        
         
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

 




########### Generate data
def generatelist(
    n: int = 1000,
    Xk: str = 'Condition_1',
    mk = None,  # transition matrix to go from condition 0 to condition k
    sd = 1, 
    nb_replicates: int = 1,
    conditionX0: list = None,
    error_rate: float = None
    ):
    
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
        
    error_rate: float
        PROPORTION OF ELEMENTS THAT SWITCH PLACE

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
        conditionX0 = "C0"
        
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
            
    
   
    dk[conditionX0] = p0
    dk[conditionX0 + 's'] = state0   
    
    
  
    
    # PART 2:     
    # for the experiment Xk, return the result vector of the experiment
    # simulate results from each condition
    
    
    # transition matrix base from condition 0
    pk, statek, tmatk, pkprimelist, statekplist, ntmatklist = transitionMatrix(\
                p0=p0, matrix=mk, sd=seed, nb_replicates=nb_replicates,
                error_rate=error_rate) 
    dk[Xk] = pk  # memorize gene expressions
    dk[Xk + 's'] = statek  # memorize states
    
    
    
    # Save the transition matrices in a dictionary
    d = {}
    d[Xk] = tmatk
    
    # For the Xkprime
    for i in range(len(statekplist)):
        dk[Xk + '_p' + str(i+1)] = pkprimelist[i]
        dk[Xk + '_p' + str(i+1) + 's'] = statekplist[i]
        d[Xk + '_p' + str(i+1)] = ntmatklist[i]
        
        
    dk = dk.round(3)
    return(dk, d)



# =============================================================================
# def generatedatadefault(nbXk, nbrep, N, sd):
#     """
#     
# 
#     Parameters
#     ----------
#     nbXk : int
#         NUMBER OF EXPERIMENTS X
#     nbrep : int
#         NUMBER OF REPLICATES FOR EACH EXPERIMENTS
#     N : int
#         SAMPLE SIZE, NUMBER OF INDIVIDUALS
#     sd : int
#         RANDOMIZATION PARAMETER. SEED
# 
#     Returns
#     -------
#     dataframe:
#         DATAFRAME OF ALL GENERATED DATA
# 
#     """
#     
#     list_Xk = []  # list of transition matrices for k experiments X
#     
#     for n in range(nbXk):
#         seed(n)
#         res = []
#         
#         for i in range(5):
#             r = [random() for _ in range(5)]
#             s = sum(r)
#             res.append([k/s for k in r])
#             restmp = np.array(res)
#         list_Xk.append(restmp)
#     
#     
#     
#     print('Generating data....')
#     print(f'Sample size : {N}')
#     print(f'Number of experiments: {nbXk}')
#     print()
#     dataframe, transmatrix = generatelist(N, Xk='Condition_1', 
#                                           nb_replicates=nbrep, mk=list_Xk[0]) 
# 
#     sd = 2
#     for mat in list_Xk[1:]:
#         name = 'Condition_' + str(sd)
#         dk, d = generatelist(N, Xk=name, sd=sd, nb_replicates=nbrep, mk=mat)
#         dataframe = pd.concat([dataframe, dk[dk.columns[2:]]], axis=1)
#         transmatrix.update(d)
#         sd += 1
# 
#     # Export data generated
#     # source = os.path.dirname(os.getcwd())
#     # path = source + '/data/'
#     # name = "testdata.csv"
#     # dataframe.to_csv(path+name)
#     # print(f'Done. Check data in {path}\n Filename: {name}')
#     
# 
#     return(dataframe)
# 
# =============================================================================

# =============================================================================
# 
# # Print out histograms of gene expressions   
# def histogram(p, title):
#     
#     # default params
#     bins=10
#     density=True
#     edgecolor='white'
#     
#     plt.figure(figsize=(10, 3))
#     N, bins, patches = plt.hist(p, 
#                                 bins=bins, 
#                                 density=density, 
#                                 edgecolor=edgecolor)
#     
#     
#     cmap = plt.get_cmap('jet')
#     low = cmap(0.1)
#     lmed = cmap(0.25)
#     medium = cmap(0.5)
#     medh = cmap(0.75)
#     high = cmap(0.95)
# 
# 
#     for i in range(0, 2):
#         patches[i].set_facecolor(low)
#     for i in range(2 , 3):
#         patches[i].set_facecolor(lmed)
#     for i in range(3, 7):
#         patches[i].set_facecolor(medium)
#     for i in range(7, 8):
#         patches[i].set_facecolor(medh)
#     for i in range(8, 10):
#         patches[i].set_facecolor(high)
# 
#     # create legend
#     handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in [low, lmed, medium, medh, high]]
#     labels = ["low", "lmedium", "medium", "mediumh", "high"]
#     plt.legend(handles, labels)
#     
#     plt.title = (f'Histogram of gene expression label obtained from {title}')
#     plt.xlabel("Gene expression level")
#     plt.ylabel("Frequency")
#     plt.xticks(fontsize=10)  
#     plt.yticks(fontsize=10)
# 
#     plt.gca().spines["top"].set_visible(False)  
#     plt.gca().spines["right"].set_visible(False)
#     
#     plt.show()
# =============================================================================

# =============================================================================
# # A function to display basic information about gene rankings
# 
# def printout(condition: str, p, n):
#     """ Function to print out condition and the associated vector """
#     print()
#     print(condition)
#     print('-'*30)
#     print(f'Propotion of state H: {sum(p >= 80) / n}')
#     print(f'Propotion of state MH: {sum(((p >= 71) & (p <= 79))) / n}')
#     print(f'Propotion of state M: {sum(((p >= 30) & (p <= 70))) / n}')
#     print(f'Propotion of state LM: {sum(((p >= 21) & (p <= 29))) / n}')
#     print(f'Propotion of state L: {sum((p < 20)) / n}')
#     print()
# =============================================================================


# =============================================================================
# # Compare generated data between groups
# def comparegeneratedlist(dataframe,
#                          group0:str = 'Condition_0_state',
#                          group:str = 'Condition_1_state'):
#     """
#     Description :
#     ----------
#     GIVEN A DATAFRAME OF DIFFERENT GENE EXPRESSIONS,
#     RETURN THE AVERAGE NUMBER OF GENES THAT SWITCHED THE CATEGORY
# 
#     Parameters
#     ----------
#     dataframe : TYPE
#         DATAFRAME OF DIFFERENT GENE EXPRESSIONS.
#         
#     group0 : str, optional
#         STATES OF A GIVEN CONDITION. The default is 'Condition_0_state'.
#     
#     group : str, optional
#         STATES OF A GIVEN CONDITION. The default is 'Condition_1_state'.
# 
#     Returns
#     -------
#     avg : 
#         COUNT NUMBER OF ELEMENTS THAT CHANGED THE CATEGORY
# 
#     """
#     
#     similarity = {}
#     
#     for i in dataframe[group].value_counts().index:
#         
#         if i in dataframe[group0].value_counts().index:
#             
#             res = (dataframe[group0].value_counts()[i] - dataframe[group].value_counts()[i])**2
#         
#         else:
#             res = dataframe[group].value_counts()[i]
#         
#         similarity[i] = res
#     
#     avg = np.sqrt(sum(similarity.values()))/5
#     
#     return(avg)
# 
# 
# =============================================================================
