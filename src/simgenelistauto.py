# -*- coding: utf-8 -*-
"""
Created on Tue Nov 15 22:12:29 2022

@author: ttran

Module to generate automatically all genelists
"""


import numpy as np
import pandas as pd
import argparse  
import os

from random import random
from random import seed
from simgenelist import generatelist





# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument('-nbXk', type=int, 
                    default = '10', 
                    help="Number of samples to generate from initial sample (positive INT)")

parser.add_argument('-nbrep', type=int,
                    default = '10', 
                    help="Number of replicates to generate for each sample \
                        (positive INT)")

parser.add_argument('-N', type=int, 
                    default = '100', 
                    help="Total number of indiv/ gene for sample generation \
                        (positive INT)")

args = parser.parse_args()

#------------------------------------------------------------------------
# RETAIN ARGUMENTS

nbXk = args.nbXk # Number of experiments k
nbrep = args.nbrep  # Number of replicates from each Xk
N = args.N
    
#------------------------------------------------------------------------



def generatedatadefault(nbXk, nbrep, N, sd, replrate):
    
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
        
    replrate: float
        PROPORTION OF ELEMENTS THAT SWITCH PLACE

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
    
    
    
# =============================================================================
#     print('Generating data....')
#     print(f'Sample size : {N}')
#     print(f'Number of experiments: {nbXk}')
#     print()
# =============================================================================
    dataframe, transmatrix = generatelist(N, Xk='C1', 
                                          nb_replicates=nbrep, mk=list_Xk[0]) 

    sd = 2
    for mat in list_Xk[1:]:
        name = 'C' + str(sd)
        dk, d = generatelist(N, Xk=name, sd=sd, nb_replicates=nbrep, mk=mat,
                             error_rate=replrate)
        dataframe = pd.concat([dataframe, dk[dk.columns[2:]]], axis=1)
        transmatrix.update(d)
        sd += 1

    # Export data generated
    source = os.path.dirname(os.getcwd())
    path = source + '/data/'
    name = "generateddata.csv"
    dataframe.to_csv(path+name)
    print(f'Done data generation.\nCheck data in {path} \n\nFilename: {name}')
    
    # Export the transition matrix
    with open(path+'Tmatrix.txt', 'w') as f:
        count = 1
        for line in list_Xk:
            f.write(f'Transition matrix of replicates of experiment {count}\n')
            f.write("\n")
            f.write(f"{line}\n")    
            f.write("="*70+"\n")
            count += 1
    
    print('Transition matrix: Tmatrix.txt')
    
            
    return(dataframe, list_Xk)



# Run the module (if called)
#dataframe, list_Xk = generatedatadefault(nbXk=nbXk, nbrep=nbrep, N=N, sd=1)