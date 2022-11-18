# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""


 # Import module
import os
import glob
import sys


import numpy as np
import argparse  # check input arguments
import pandas as pd  # handle dataframe



# Packages to generate rankings
from simgenelistauto import generatedatadefault

# from generatedata import simulation

# Call module to calculate RSM
from rbomod import rbo_modified as RSM


# Call module simulation
from simRSM import simulation as simRSM

# Call function to fit distribution to the simulation
from distfit import distfit
import fitter


# Check for skewness
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy.optimize import newton

# check for histogram
import matplotlib.pyplot as plt


# Import random
from random import random
from random import seed




#=============================================================================================================
if __name__=='__main__':

    # Set arguments 
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    

    parser.add_argument('-data', type=str, default = 'N', 
                        help="If the option data is Y, then the user will provide the data (under the csv format). If not, then the user can give parameters nbXk, nbrep, N, to simulate data")

    parser.add_argument('-nbXk', type=int, default = '1', 
                        help="Number of samples to generate from initial sample (positive INT)")

    parser.add_argument('-nbrep', type=int, default = '10', 
                        help="Number of replicates to generate for each sample (positive INT)")

    parser.add_argument('-N', type=int, default = '1000', 
                        help="Total number of indiv/ gene for sample generation (positive INT)")

    parser.add_argument('-cat', type=str, default = 'H', 
                        help="Category. H or MH or M or LM or L")

    parser.add_argument('-expid', type=str, default = '1', 
                        help="The id of experiment to test on (STR)")

    parser.add_argument('-topk', type=int, default = '20', 
                        help="Top k first elements within the \
                        category to check(positive INT)")

    parser.add_argument('-strict_k', type=bool, default = 'False', 
                        help="If True, the user only cares about the top k elements")
                        
    parser.add_argument('-wg_func_name', type=str, default = 'poisson', 
                        help="Weighting function to be chosen from, among the list ['geometric', 'poisson', 'skellam', 'triangular']")
                        
    parser.add_argument('-p_geometric', type=str, default='None', 
                        help="Parameters associated with a certain weighting function \
                        default = {'geometric': {'p': 0.5}, \
                                   'poisson'  : {'lambda': k},\
                                   'skellam'  : {'mu1': 2*k, 'mu2': k}, \
                                   'triangular': {'topk': k}} ")
    parser.add_argument('-lamda_poisson',type=str, default = 'None', 
                        help="Parameters associated with a certain weighting function \
                        default = {'geometric': {'p': 0.5}, \
                                   'poisson'  : {'lambda': k},\
                                   'skellam'  : {'mu1': 2*k, 'mu2': k}, \
                                   'triangular': {'topk': k}} ")

    parser.add_argument('-mu1_skellam',type=str, default = 'None', 
                        help="Parameters associated with a certain weighting function \
                        default = {'geometric': {'p': 0.5}, \
                                   'poisson'  : {'lambda': k},\
                                   'skellam'  : {'mu1': 2*k, 'mu2': k}, \
                                   'triangular': {'topk': k}} ")
                            
    parser.add_argument('-mu2_skellam',type=str, default = 'None', 
                        help="Parameters associated with a certain weighting function \
                        default = {'geometric': {'p': 0.5}, \
                                   'poisson'  : {'lambda': k},\
                                   'skellam'  : {'mu1': 2*k, 'mu2': k}, \
                                   'triangular': {'topk': k}} ")

    parser.add_argument('-topk_triangular', type=str,default = 'None', 
                        help="Parameters associated with a certain weighting function \
                        default = {'geometric': {'p': 0.5}, \
                                   'poisson'  : {'lambda': k},\
                                   'skellam'  : {'mu1': 2*k, 'mu2': k}, \
                                   'triangular': {'topk': k}} ")
                            
                            
    parser.add_argument('-simRSM', type=bool, default = 'True', 
                        help="If True, we will generate random lists from which \
                        we calculate RSM scores, then build RSM scores distribution.")
                        
    parser.add_argument('-withreplacement', type=bool, default = 'True', 
                        help="Parameter that sets if we expect replacement in the sample") 
                        
    parser.add_argument('-percent_repl', type=float, default = '0.25', help="Replacement rate")
                        
    args = parser.parse_args()


    # RETAIN ARGUMENTS
# =============================================================================
    dataopt = args.data
    nbXk = args.nbXk # Number of experiments k
    nbrep = args.nbrep  # Number of replicates from each Xk
    N = args.N
    
    condition = args.cat
    expid = args.expid
    topk = args.topk
    strict_k = args.strict_k


    # Weighting function parameter
    wg_func_name = args.wg_func_name
    p = args.p_geometric
    lamb = args.lamda_poisson
    mu1 = args.mu1_skellam
    mu2 = args.mu2_skellam
    topk_triangular = args.topk_triangular

    # RSM simulation parameters for distribution
    simRSMtest = args.simRSM
    repltest = args.withreplacement
    replrate = args.percent_repl


    # Create the dictionary of parameters for weighting function:
    if (p == 'None') & (lamb == 'None') & (mu1 == 'None') & (mu2 == 'None') & \
    (topk_triangular == 'None') :
        print('User default parameters of weighting function')
        params_wg_func = None

    elif (p != 'None'):
        params_wg_func = {'p': p}
    
    elif (lamb != 'None'):
        params_wg_func = {'lambda': lamb}
    
    elif (mu1 != 'None') | (mu2 != 'None'):
        params_wg_func = {'mu1': mu1, 'mu2': mu2}
    
    elif (topk_triangular != None):
        params_wg_func = {'topk': topk_triangular}
    
# =============================================================================
# Clear folder for data
# =============================================================================



    source = os.path.dirname(os.getcwd())
    path = source + '/data/'
    files = glob.glob(path+'*')
    for f in files:
        os.remove(f)
        
    seed(0)
# =============================================================================
# Load data
# =============================================================================

    if dataopt == "Y":
        print('Please put the data you want to load in the current directory')
        path = os.path.dirname(os.getcwd())
        print('We are currently at:', path)
        print()
        filename = input("Enter the full file name (format .csv required)")
        dataframe = pd.read_csv(path+filename, index_col=0)
    else :
        print('No data loaded, will generate data using original parameter')
        dataframe, list_Xk = generatedatadefault(nbXk=nbXk, nbrep=nbrep, N=N, sd=0, replrate=replrate)

# =============================================================================

# Display the dissimilarity of each pairwise element for the selected sample

    c = (dataframe[f'C{expid}s'] == condition)
    l = list(dataframe[c].sort_values(f'C{expid}', ascending=False).index)


# Display
    for i in range(1, nbrep+1):
     
        print()
        print('-'*50)
        print(f'Condition {expid} vs replicate {i}')
        print(f'Top {topk} ranked genes')
        print('.'*30)
        cprime = (dataframe[f'C{expid}_p{i}s'] == condition)    
        lprime = list(dataframe[cprime].sort_values(f'C{expid}_p{i}', 
                                                ascending=False).index)
    
     # Display the topk elements
    for k in range(topk):
        print(f'Rank {k+1} : {l[k]} --- {lprime[k]} --> Equal: {l[k] == lprime[k]} ') 
        # =============================================================================
        # Calculate RSM scores for the gene rankings

    numcols = [i for i in dataframe.columns if 's' not in i]

    # Create dataframe for RSM scores
    df_rsm = pd.DataFrame(index=numcols, columns = numcols)
    df_rsm = df_rsm.replace(np.nan, 1.)

    # Calculate RSM scores for the given category
    for lstS in numcols:
        for lstT in numcols:
            if lstS != lstT :
                #print(lstS, lstT)
                
                S = list(dataframe[lstS].sort_values(ascending=False).index)
                T = list(dataframe[lstT].sort_values(ascending=False).index)
            
                # Calculate RSM score
                rsm = RSM(lstS=S, lstT=T, wg_func_name=wg_func_name,
                          params_wg_func=params_wg_func, k=topk, strict_k=strict_k)
                df_rsm.loc[lstS, lstT] = rsm
                df_rsm.loc[lstT, lstS] = rsm
                #print('S=', lstS[:10], '\nT=', lstT[:10], '\n', rsm)
                # Add a subset of rsm values
            
    # Export data generated
    name = f"RSM{wg_func_name}.csv"
    df_rsm.to_csv(path+name)
    print(f'\nCheck RSM scores in {path} \n\nFilename: {name}')

    # Plot out a histogram
    plt.figure()
    plt.hist(df_rsm[f'C{expid}'])
    plt.title(f'Histogram of RSM scores obtained from data of Experiment {expid}')
    plt.savefig(path+'HistRSMscoresdata.png')


    # Fit a distribution to the ranked gene
    print(f'\nFitting distribution to the gene ranking of Experiment {expid}')

    distrlist = ['norm', 'expon', 'lognorm', 'gamma', 'dweibull', 
             'logistic','pareto', 't', 'burr', 'chi2', 'frechet_r[x]]', 'cauchy']

    rsm_scoresdata = df_rsm[f'C{expid}']

    dist = distfit(todf=True, distr=distrlist)
    results = dist.fit_transform(np.array(rsm_scoresdata))
    dist.plot()
    dist.plot_summary()
    df_paramsdata = results['summary']

    f = fitter.Fitter(rsm_scoresdata, distributions=distrlist)
    f.fit()
    df_nonparamsdata = f.summary()

    # Display fitted distribution
    print('Parametric fit for data\n', df_paramsdata)
    print()
    print('Non-parametric fit for data\n', df_nonparamsdata)
    # Export data generated
    source = os.path.dirname(os.getcwd())
    path = source + '/data/'
    
    # For parametric estrimation
    name = f"FitDist{wg_func_name}_P_data.csv"
    df_paramsdata.to_csv(path+name)
    print(f'\nCheck fitted distribution in {path} \n\nFilename: {name}')
    
    # For non parametric
    name = f"FitDist{wg_func_name}_NP_data.csv"
    df_nonparamsdata.to_csv(path+name)
    print(f'\nCheck fitted distribution in {path} \n\nFilename: {name}')


    # Check the skewness of generated data
    print('The skewness of the RSM scores from data =', skew(rsm_scoresdata))
    print('The kurtosis of the RSM scores from data =', kurtosis(rsm_scoresdata))

    
  
    




# =============================================================================
# Fit distribution to the generated data
# =============================================================================

    if simRSMtest == True:
        print('Generating data and fitting distribution')
        
        size_l = N 
        distrlist = ['norm', 'expon', 'lognorm', 'gamma', 'dweibull', 
                     'logistic','pareto', 't', 'burr', 'chi2', 'frechet_r[x]]', 'cauchy']
    
        dist = distfit(todf=True, distr=distrlist)
        # Set up parameters for the replacement/ permutation
        percent_repl =  replrate
        percent_shuffle = replrate
        rsm_rnksim, rsm_scoressim = simRSM(topk, strict_k, size_l, repltest, wg_func_name, 
                                           params_wg_func, percent_repl=percent_repl, 
                                           percent_shuffle=percent_shuffle)
        
        results = dist.fit_transform(np.array(rsm_scoressim))
        dist.plot()
        dist.plot_summary()
        df_params = results['summary']
        
        f = fitter.Fitter(rsm_scoressim, distributions=distrlist)
        f.fit()
        df_nonparams = f.summary()
        
        # Display fitted distribution
        print('Parametric fit\n', df_params)
        print()
        print('Non-parametric fit\n', df_nonparams)
        # Export data generated
        source = os.path.dirname(os.getcwd())
        path = source + '/data/'
        
        # For parametric estrimation
        name = f"FitDist{wg_func_name}_P.csv"
        df_params.to_csv(path+name)
        print(f'\nCheck fitted distribution in {path} \n\nFilename: {name}')
        
        # For non parametric
        name = f"FitDist{wg_func_name}_NP.csv"
        df_nonparams.to_csv(path+name)
        print(f'\nCheck fitted distribution in {path} \n\nFilename: {name}')
        
        
        # Check the skewness of generated data
        print('The skewness of the RSM scores generated =', skew(rsm_scoressim))
        print('The kurtosis of the RSM scores generated =', kurtosis(rsm_scoressim))
        
        plt.figure()
        plt.hist(rsm_scoressim, bins=20)
        plt.title('Histogram of the generated RSM scores')
        plt.savefig(path+'HistRSMscoressim.png')
        
        
        
    

# =============================================================================
# Finally, run the estimation on the quantile
# Extreme value distribution

# Import extreme value distribution
    from extremevaluedistr import maxlogvraissemblance
    from extremevaluedistr import quantileExtreme

    g_0 = 0
    a_0 = 0.1
    b_0 = 0.9
    
    g_m, a_m, b_m = newton(maxlogvraissemblance, [g_0, a_0, b_0], maxiter=300)
    c = 1
    

    while (sum(np.isfinite([g_m, a_m, b_m])) < 3):
        
        print('Iteration: ', c,'-th')
        g_m, a_m, b_m = newton(maxlogvraissemblance, [g_0, a_0, b_0], maxiter=300)
        g_0 = random()
        a_0 = random()
        b_0 = random()
        
        c += 1
        #print(g_0, a_0, b_0)
        
        
        
    print('Result of approximation')
    print('gamma=', g_m)
    print('a_m=', a_m)
    print('b_m=', b_m)

    if g_m == -1:
        print('Warning, this estimation is not consistent')
    
    p_m = 0.05
    m = 10
    q_m = quantileExtreme(a_m, b_m, g_m, p_m, m)
    print('Estimated extreme quantile q=', q_m, 'for p=', p_m)
    



