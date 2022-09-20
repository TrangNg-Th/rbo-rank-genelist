# MODULE TO PLOT OUT MARKOV CHAIN
# Simulate lists of elements to test
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle
import random 
import pandas as pd



from markovchain import MarkovChain


def transitionMatrix(
    p0: list = None,
    states = ['H', 'MH', 'M', 'LM', 'L'],
    matrix = None,
    d : dict = None,
    error_rate : float = None, 
    nb_replicates: int = 1,
    seed=1):
    """
    
    Function to generate pk given a transition matrix and p0
    Arg:
        p0: list of expression values of baseline condition
        states : list of states
        matrix : transition matrix
        d : dictionary of each states with their associated scale
        error_rate : if error rate if False, use default "transition matrix" to go from Xk to Xk_prime
        nb_replicate : number of replicates the user wants to produce
        
    return:
        pk: list of expression values of baseline condition after the transition matrix
        pkprimelist: list of list of expression values (pkprime) of Xk after the noise transition matrix 
    """
    
    random.seed(seed)
    
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
    tmat = pd.DataFrame(matrix, index=states, columns=states)  # transition matrix from condition 0 to Xk
    
    # give pk for condition k
    pk, state = [] , []
    for i in range(len(p0)):
        if p0[i] >= 80 : 
            state_i = random.choices(states, weights=tmat[states[0]].values)[0]
           
        elif (p0[i] >= 70) & (p0[i] < 80):
            state_i = random.choices(states, weights=tmat[states[1]].values)[0]
            
        elif (p0[i] >= 30) & (p0[i] < 70):
            state_i = random.choices(states, weights=tmat[states[2]].values)[0]
           
        elif (p0[i] >= 20) & (p0[i] < 30):
            state_i = random.choices(states, weights=tmat[states[3]].values)[0]
            
        elif (p0[i] < 20) :
            state_i = random.choices(states, weights=tmat[states[4]].values)[0]
        pk.append(random.uniform(d[state_i][0], d[state_i][1]))
        state.append(state_i) 
    
    
    print()
    print('Displaying the transition matrix')
    print(tmat)
            
    # --------------------------
    # give list of pkprime from pk   
    pkprimelist, stateprimelist = [], []
    
    
    
    # if noise was requested    
    for rep in range(nb_replicates):
        pkprime, stateprime = [], []
        
        if error_rate is None: 
            nmatrix = ToleranceTransMat()
        
        else:
            nmatrix = ToleranceTransMat(error_rate=error_rate)
        ntmat = pd.DataFrame(nmatrix, index=states, columns=states)  # noise transition matrix
     
        for i in range(len(pk)):
            if pk[i] >= 80 : 
                state_i = random.choices(states, weights=ntmat[states[0]].values)[0]
                
            elif (pk[i] >= 70) & (pk[i] < 80):
                state_i = random.choices(states, weights=ntmat[states[1]].values)[0]
            
            elif (pk[i] >= 30) & (pk[i] < 70):
                state_i = random.choices(states, weights=ntmat[states[2]].values)[0]
            
            
            elif (pk[i] >= 20) & (pk[i] < 30):
                state_i = random.choices(states, weights=ntmat[states[3]].values)[0]
           
            elif (pk[i] < 20) :
                state_i = random.choices(states, weights=ntmat[states[4]].values)[0]
            
            # if there's a change in the state of the current gene
            if state_i != state[i]:
                pkprime.append(random.uniform(d[state_i][0], d[state_i][1]))
            else:
                pkprime.append(pk[i])
            stateprime.append(state_i)
        
        pkprimelist.append(pkprime)  # list of all the noisy replicates
        stateprimelist.append(stateprime)
        
        print()
        print('Displaying the noise transition matrix')
        print(ntmat)

    
    # Fix the Markov graph !!!
    #mc = MarkovChain(matrix, states)
    #mc.draw()
    return(pk, state, tmat, pkprimelist, stateprimelist, ntmat)

    
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
    seed = 1, 
    nb_replicates: int = 1,
    conditionX0: list = None):
    
    """
    Design : 
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
            
        
    Argument : 
        n     -- number of gene in the sample
        seed  -- set seed for reproducibility of random choice
        Xk: experiment name
        conditionX0: if not None, it should be a vector (list) of genes expression values of size n.
        mk : transition matrix to go from conditionX0 to conditionXk
        nb_replicates : number of replicates of Xk
        
    Return :
        conditionX0 :   if not None, it should be a vector (list) of genes expression values of size n.
                        if None, it will return the hardcoded baseline condition
    
        p0 : vector of gene expression values of conditionX0
        dk : dataframe of vector of gene expression values of conditionXk and conditionX0
        
    
    """
    
    random.seed(seed)
    dk = pd.DataFrame(index=[f'gene_{i}' for i in range(n)])
    
    # PART 1 :    
    # generate baseline condition (if not given)
    if conditionX0 == None:
        
        print("Using default condition 0")
        conditionX0 = "Condition_0"
        
        pH = np.array(random.choices(range(80, 101), k = int(0.1 * n)))
        pM = np.array(random.choices(range(30, 70), k = int(0.8 * n)))
        pL = np.array(random.choices(range(0, 21), k = int(0.1 * n)))
        p0 = np.concatenate((pH, pM, pL), axis=None)  
        
        
    else : 
        
        p0 = np.array([])
        for k in exppr.keys():
            p0 = np.concatenate((p0, np.array(exppr[k])), axis=None)
            
            
    # Display info
    print('_'*50)
    printout(conditionX0, p0, n)
    #histogram(p0, title = conditionX0)
    dk[conditionX0] = p0
    
    
    
  
    
    # PART 2:     
    # for the experiment Xk, return the result vector of the experiment
    # simulate results from each condition
    pk, statek, _, pkprimelist, statekp, _ = transitionMatrix(p0=p0, matrix=mk, seed=seed, 
                                                              nb_replicates=nb_replicates) # transition matrix base from condition 0
    dk[Xk] = pk
    dk[Xk+'_state'] = statek
    
    # Displaying results
    printout(Xk, dk[Xk], n)
    #histogram(pk, title=Xk)
    
    # For the Xkprime
    for i in range(len(statekp)):
        dk[Xk+'_prime'+str(i+1)] = pkprimelist[i]
        dk[Xk+'_primestate'+str(i+1)] = statekp[i]
        printout(Xk+'prime'+str(i+1), dk[Xk+'_prime'+str(i+1)], n)
    return(dk)




## Given an error rate,  return transition matrix
def ToleranceTransMat(
    error_rate : float = 0.05,
    states : list = ['H', 'MH', 'M', 'LM', 'L']):
    """
    
    Given an error rate, return the transition matrix of size 5 x 5
    
    """
    print()
    print(f'No error rate was given, default value = {error_rate}')
    print()
    matrix = np.empty([len(states), len(states)])
    for j in range(len(states)):
        for i in range(len(states)):
            if i != j:  # if the current position is not on the diagonal of the matrix
                matrix[i,j] = error_rate / (len(states) - 1)  # divide evenly the error_rate to this position
            else:
                matrix[i, j] = 1 - error_rate  # if we're on the diagonal, we'll have the biggest possibility
                
    return(matrix)
            
    
    
    
