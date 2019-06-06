from __future__ import division #This means division always gives a floating result
import plot_module as plm
import numpy as np
import pandas as pd;
import sys
from scipy.stats import linregress, ks_2samp,poisson;
from math import *
from scipy.interpolate import interp1d



def compute_likelihood_model(directory,results_path, population_data,merged_dataframe,model,low_res=False):
 
 # This function computes the likelihood of a model and the most likely values for the model parameters. Data is stored in results path
 # The procedure can compute the likelihood of three different classes of model: the epidemiological model, a linear model where the infected population is proportional
 # to population density and a constant model where the expected size of the infected population is constant.The input data is in the merged data_frame (Note for Camille: this is the data we gave Eriksson - should be possible to simplify)
 # When the low_res parameter is set to true, the system produces low_res graphs. Used for system testing and exploratory testing
 # Set up the ranges of values to be tested for the three parameters in the model
 # In low_res we explore the same ranges as in high_res but with many fewer samples
 # The parameter values are set up in the info files for the population_dats
 # lambda is the death rate
 # Zetta is the base probability that a  territorial unit) contains at least one site
 # Eps is an error - means that there is a positive probability that a site is present even in a territory with below threshold population
 
 
    if low_res:
        lambda_v=np.linspace(population_data.likelihood_parameters[0],population_data.likelihood_parameters[1],num=24)
        zetta_v=np.exp(np.linspace(log(population_data.likelihood_parameters[2]),log(population_data.likelihood_parameters[3]),num=11,endpoint=False)) 
        eps_v=np.linspace(population_data.likelihood_parameters[4],population_data.likelihood_parameters[5],num=11,endpoint=False)#eriksson only
  #      comm_v=np.linspace(0,2000,num=20)
    else:
        lambda_v=np.linspace(population_data.likelihood_parameters[0],population_data.likelihood_parameters[1],num=101)
        zetta_v=np.exp(np.linspace(log(population_data.likelihood_parameters[2]),log(population_data.likelihood_parameters[3]),num=101,endpoint=False))
        eps_v=np.linspace(population_data.likelihood_parameters[4],population_data.likelihood_parameters[5],num=101,endpoint=False)
  #      comm_v=np.linspace(0,1,num=1)
    if model=='constant':
        eps_v=np.array([1])
 # rho_bins are the bins where we count number of samples and controls. Note we use 3001 bins - probably more than necessary. 
 # this gives the highest possible resolution using the Eriksson data
    rho_bins=np.linspace(0,33,num=3001,endpoint=False)
 # The bins in the python histogram function are open intervals - in matlab they are closed. We add the additional point to make sure the code is
 # identical in both languages
 
    rho_bins_4_python=np.append(rho_bins,33)
# When we show the actual frequencies in Figure 1 - we use smaller bins - with a larger number of samples per bin. This makes the graph easier to read
# bin_boundaries2_4_python contains the second type of bis
    bin_width=2 
    bin_boundaries2_4_python=np.linspace(0,33,num=16,endpoint=False) 
    bin_boundaries2_4_python=np.append(bin_boundaries2_4_python,33)
# rho_bins2 contains actual bins used for graphing actual frequencies - ugly and could be improved
    rho_bins2=bin_boundaries2_4_python[0:len(bin_boundaries2_4_python)-1]+bin_width/2 #This is horribly complicated - needs to be cleaned up
# Count the number of samples and controls in each bin - in Timmermann hi res these are giving zero counts
    samples_counts=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==1],bins=rho_bins_4_python)[0] #Column  vector This is not strict translation of mathlab code. In mathlab the last bin contains 3000. In python it contains 2999-3000
    controls_counts=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==0],bins=rho_bins_4_python) [0] #Column  vector
# Repeat for the larger bins
    control_counts2=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==0],bins=bin_boundaries2_4_python)[0]
    sample_counts2=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==1],bins=bin_boundaries2_4_python)[0]
# Compute total number of controls and samples
    n_controls=np.sum(controls_counts)  
    n_samples=np.sum(samples_counts) 
# Meant to avoid underflow in calculations but not fully understood.
    l_shift=n_samples*(log(float(n_samples)/float(n_controls))-1)
# Computes size of parameter ranges according to range actually chosen - could be cleaned up
    n_lambda=len(lambda_v)
    n_eps=len(eps_v)
    n_zetta=len(zetta_v)
 #   n_comm=len(comm_v)
#  Kills lambda loop for constant and linear models
    if model=='constant' or model=='linear':
        n_lambda=1
# Set up a (population data specific) range of possible values for the likelihood of a given set of observations
    acc_likelihoods=np.linspace(population_data.likelihood_parameters[6],population_data.likelihood_parameters[7],num=2001) 
#  Set up an array representing the accumulated likelihood of a given set of sample and control counts 
#  Across all possible values of the parameters
    acc=np.zeros((len(acc_likelihoods),len(rho_bins)))
    lnL=np.zeros((n_lambda,n_eps,n_zetta))
    sqrt_rho_bins=np.sqrt(rho_bins) #These are the values we are computing - rhobins_4_python are intervals for histogram only. In original program were inside loop. Have moved it outside
    max_LL=-float('inf') 
    bin_zeros=np.zeros(len(sqrt_rho_bins))
# Scan all possible values of the parameters
    for i_lambda in range (0,n_lambda):
        print 'Percentage completed=', i_lambda/float(n_lambda)
        my_lambda=float(lambda_v[i_lambda]) 
# Compute the predicted size of the infected population  as a proportion of population (p_infected) for all possible values of rho_bins, given the value of lambda. Guarantee it is always 0 or greater
        if model=='epidemiological':
            p_infected=np.maximum(bin_zeros,1-my_lambda/sqrt_rho_bins)  #COLUMN VECTOR
        else:
            if model=='richard':
#                i_star=np.maximum(bin_zeros,rho_bins-my_lambda*sqrt_rho_bins)
                p_infected=np.maximum(bin_zeros,1-my_lambda/sqrt_rho_bins)  
            if model=='linear' or model=='constant':
                p_infected=rho_bins/max(rho_bins)
# Scan all possible values of lambda, zetta and eps
        for i_zetta in range(0,n_zetta):
            for i_eps in range (0, n_eps):
                 my_zetta=zetta_v[i_zetta]
                 my_eps=eps_v[i_eps]
# =============================================================================
#                  for i_comm in range(0, n_comm):
#                      my_comm=comm_v[i_comm]
#                      comm_mult=86.655/my_comm
# =============================================================================
# Predicts the probability of finding at least one site in a territory for each of the population densities in rho_bins. Make sure value is never too small
                 if model=='epidemiological':
                     p_predicted=compute_epidemiological_model(p_infected,my_zetta,my_eps)
                 else:
                     if model=='linear':
                         p_predicted=compute_linear_model(p_infected,my_zetta,my_eps)
                     else:
                         if model=='constant':
                             p_predicted=compute_constant_model(p_infected,my_zetta,my_eps)
                         else:
                             if model=='richard':
                                 p_predicted=compute_richard_model(p_infected,rho_bins,my_zetta,my_eps)
 #                                print 'p_predicted=', p_predicted
                                
                             else:
                                 print "Model=",model
                                 print "No model available"
                                 sys.exit()
# Computes the log likelihood of obtaining the OBSERVED number of samples at a given value of rho_bins, given the predicted number of samples
                  
                 log_samples=np.dot(samples_counts,np.log(p_predicted)) 
# The same for controls
# Problem occurs when log_controls=NaN
                 
                 # p_predicted is giving values higher than 1 which causes calculation of negative log)
                 log_controls=np.dot(controls_counts,np.log(1-p_predicted)) 
                 if np.isnan(log_controls):
                     print 'controls counts=',controls_counts
                     print 'p_predicted=', p_predicted
                     print 'my_lambda=',my_lambda
                     print 'my_zetta=',my_zetta
                     print 'my_eps=',my_eps
                     print 'p_infected=',p_infected
                     sys.exit()
# Computes the log likelihood of a certain number of samples AND a certain number of controls (for a given value of rho_bins)
                 LL=log_samples+log_controls
# Finds the parameter values with the highest likelihood
                 if LL>max_LL:
                     max_LL=LL
                     max_lambda=my_lambda
                     max_zetta=my_zetta
                     max_eps=my_eps
 #                    max_comm=my_comm
                     max_likelihood=LL
                 if np.isnan(np.min(LL)):
                     print 'LL is nan'
                     print 'nSamples=',n_samples
                     print 'nControls=',n_controls
                     print 'my_lambda=',my_lambda
                     print 'my_zetta=',my_zetta
                     print 'my_eps=',my_eps
                     print 'log_samples',log_samples
                     print 'log_controls', log_controls
                     print 'p_predicted=', p_predicted
                     print 'samples_counts=', samples_counts
                     print 'controls_counts=',controls_counts
                     sys.exit()
                         
# Stores the log likelihood in an array indexed the position of the parameter values in the parameter ranges
# =============================================================================
#                  if np.isnan(LL)==False:
                 lnL[i_lambda,i_eps,i_zetta]=LL 
# =============================================================================
# Computes the actual likelihood of the observations and applies a left shift to make sure it is not too large (This means values shown are relative only)
                 L=np.exp(LL-l_shift)
                 len_acc_likelihoods=np.array(len(acc_likelihoods))
                 len_acc_likelihoods.fill(len(acc_likelihoods))
# Create a one dimensional array of indexes pointing to values in acc_likelihoods (e.g. possible likelihood values) corresponding to different values of pObs,COMPLEX - WOULD BE NICE TO HAVE EASIER APPROACH. 
                 i_acc=np.minimum(len_acc_likelihoods,np.floor(1+p_predicted/acc_likelihoods[1]).astype(int)) #This yields column vector of indexes corresponding to different values of pObs. ector length =401. Maximum value of index =400 (zero based vector).I am keeping it 1-based
                 for i in range(0,len(i_acc)):
                     x_coord=i_acc[i]-1
                     y_coord=i
# Accumulate likelihood values (x coord) for a each possible value of rho_bins (y_coord) across all values of the parameters
                     acc[x_coord,y_coord]=acc[x_coord,y_coord]+L
# Compute threshold from model - used in grap
 #   opt_threshold=max_lambda**2  #Not sure about this
# Plot maximum likelihood graph
# Plot graphs for most likely values of each parameter
 #   print 'lnl going into plots=',lnL
    interpolated_lambdas=plm.plot_parameter_values(lnL,lambda_v, zetta_v, eps_v,model,directory,results_path)
    opt_threshold=interpolated_lambdas[2]**2  #Not sure about this
    print 'opt_threshold=',opt_threshold
    plm.plot_maximum_likelihood(acc,rho_bins,rho_bins2,acc_likelihoods, lambda_v, opt_threshold, sample_counts2, control_counts2, model,directory,results_path)
    return(max_lambda, max_zetta, max_eps, max_likelihood,interpolated_lambdas)
    
def compute_epidemiological_model(p_infected, my_zetta,my_eps):
    p_predicted=np.zeros(len(p_infected)).astype(float) 
    p_predicted=my_zetta*((1-my_eps)*p_infected+my_eps)
    p_predicted_small=np.zeros(len(p_predicted))
    p_predicted_large=np.zeros(len(p_predicted))
    p_predicted_small.fill(1e-20)
    p_predicted_large.fill(1-0.000000001)
    p_predicted=np.maximum(p_predicted,p_predicted_small)
    p_predicted=np.minimum(p_predicted,p_predicted_large)
    p_predicted=p_predicted.astype(float) #Probably not necessary
    return(p_predicted)
    
def compute_richard_model(p_infected,rho_bins,my_zetta,my_eps):
    p_predicted=np.zeros(len(rho_bins)).astype(float) 
    p_predicted=my_zetta*((1-my_eps)*p_infected*rho_bins)+my_eps
    p_predicted_small=np.zeros(len(p_predicted))
    p_predicted_large=np.zeros(len(p_predicted))
    p_predicted_large.fill(1-0.000000001)
    p_predicted_small.fill(1e-20)
    p_predicted=np.maximum(p_predicted,p_predicted_small)
    p_predicted=np.minimum(p_predicted,p_predicted_large)
    p_predicted=p_predicted.astype(float) #Probably not necessary
    return(p_predicted)
    
def compute_linear_model(p_infected, my_zetta,my_eps):
    p_predicted=np.zeros(len(p_infected)).astype(float) 
    p_predicted=my_zetta*((1-my_eps)*p_infected+my_eps)
    p_predicted_small=np.zeros(len(p_predicted))
    p_predicted_small.fill(1e-20)
    p_predicted=np.maximum(p_predicted,p_predicted_small)
    p_predicted=p_predicted.astype(float) #Probably not necessary
    return(p_predicted)
    
def compute_constant_model(p_infected, my_zetta,my_eps):
    p_predicted=np.zeros(len(p_infected)).astype(float) 
    p_predicted=my_zetta*((1-my_eps)*p_infected+my_eps)
    p_predicted_small=np.zeros(len(p_predicted))
    p_predicted_small.fill(1e-20)
    p_predicted=np.maximum(p_predicted,p_predicted_small)
    p_predicted=p_predicted.astype(float) #Probably not necessary
    return(p_predicted)
                 

def process_dataframe(dataframe):

    conditions = [];
    samples_growth_coefficients = []
    valid_ids = []

    cluster_ids = dataframe.cluster_id.unique();


    for cluster_id in cluster_ids:
        cluster_df = dataframe[dataframe.cluster_id==cluster_id]
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        if np.isnan(sample_cluster_df['density'].median()):
            print("Removing cluster " + str(cluster_id));
            dataframe = dataframe[dataframe.cluster_id != cluster_id];
            continue;

        # extract all periods and all population as arrays from the samples dataframe
        sample_times = sample_cluster_df['period'].values
        sample_populations = sample_cluster_df['density'].values

        # compute growth coefficients for samples
        growth_coefficient_samples=compute_growth_coefficient(sample_times, sample_populations)

        valid_ids.append(cluster_id);
        samples_growth_coefficients.append(growth_coefficient_samples)

    for cluster_id in valid_ids:
        conditions.append((dataframe['cluster_id'] == cluster_id));

    dataframe['samples_growth_coefficient'] = np.select(conditions, samples_growth_coefficients);
    return dataframe;

def compute_growth_coefficient(times, populations):
    if len(times)>=2:
        for i in range(0, len(populations)):
            if np.isnan(populations[i]):
                populations[i] = 0
        slope, intercept, r_value, p_value, std_err = linregress(times, populations)
        return slope 
    else:
        return -1.


def generate_bin_values_dataframe(dataframe, globals_dataframe, population_data, minimum_globals):
    bin_size = population_data.bin_size
    max_population = population_data.max_population
    # minimum_bin=max_for_uninhabited
    minimum_bin = 0
    bins_to_omit=int(minimum_bin/bin_size)

    ######################
    # Create bin columns #
    ######################
    # creating bins according to density of the row

    # main dataframe
    dataframe['bin_index'] = (dataframe.density/bin_size)-bins_to_omit
    dataframe['bin_index'] = dataframe.bin_index.astype(int)
    dataframe = dataframe[dataframe.bin_index >= 0]
    dataframe['bin'] = dataframe.bin_index*bin_size+minimum_bin

    # globals dataframe
    globals_dataframe['bin_index'] = (globals_dataframe.density/bin_size)-bins_to_omit
    globals_dataframe['bin_index'] = globals_dataframe.bin_index.astype(int)
    #we add bin_size/2 to get midpoint of each bin
    globals_dataframe['bin'] = globals_dataframe.bin_index*bin_size+minimum_bin
    bin_array = []
    sample_counts = []
    global_counts = []
    likelihood_ratios = []
    p_samples = []
    p_globals = []

    ##############
    # Get Totals #
    ##############
    # total samples by summing contributions
    # total globals by counting rows
    total_samples = dataframe[dataframe.type=='s']['contribution'].sum()
    total_globals = globals_dataframe['density'].count()
    
    #########################
    # Loop through each bin . data untrimmed - would be better to filter here#
    #########################
    current_bin = minimum_bin
    while(current_bin < max_population):
        
        bin_array.append(current_bin)
        # sample count: for all samples in the bin, sum all contributions
        samples_dataframe = dataframe[dataframe.type=='s']
        current_sample_count = samples_dataframe[samples_dataframe.bin == current_bin]['contribution'].sum()
        if np.isnan(current_sample_count):
            current_sample_count = 0;
        sample_counts.append(current_sample_count)
        
        # global count: count all globals dataframe rows in the bin
        current_global_count = globals_dataframe[globals_dataframe.bin == current_bin]['density'].count()
        if np.isnan(current_global_count):
            current_global_count = 0;
        global_counts.append(current_global_count)
        
        # likelihood ratio: sample_count/global_count - probably no lomger necessary
        likelihood_ratio = -1
        if(current_global_count != 0):
            likelihood_ratio = float(current_sample_count)/current_global_count
        likelihood_ratios.append(likelihood_ratio)
        
        # p_sample: sample_count/total_samples
        p_sample = -1
        if total_samples > 0:
            p_sample = float(current_sample_count)/total_samples
        p_samples.append(p_sample)
        
        # p_global: global_count/total_globals
        p_global = -1
        if total_globals > 0:
            p_global = float(current_global_count)/total_globals
        p_globals.append(p_global)

        current_bin += bin_size

    df = pd.DataFrame({'bin_array': bin_array, 'sample_counts': sample_counts, 'global_counts': global_counts, 'likelihood_ratios': likelihood_ratios, 'p_samples': p_samples, 'p_globals': p_globals})

    return df;

def generate_statistics(dataframe, globals_dataframe, bin_values_df, minimum_globals):

    trimmed_bin_values_df = bin_values_df[bin_values_df.global_counts > minimum_globals];
    trimmed_bin_values_df['cum_p_samples'] = trimmed_bin_values_df.p_samples.cumsum();
    trimmed_bin_values_df['cum_p_globals'] = trimmed_bin_values_df.p_globals.cumsum();

    
    if len(trimmed_bin_values_df.index) < len(bin_values_df.index)/2:
        return None, None;

    stat_dictionary = {};

    stat_dictionary['trimmed_bin_values_df'] = trimmed_bin_values_df;

    stat_dictionary['total_samples'] = dataframe[dataframe.type=='s']['density'].count()
    stat_dictionary['total_globals'] = globals_dataframe ['density'].count()

    stat_dictionary['median_samples'] = dataframe[dataframe.type=='s']['density'].median()
    stat_dictionary['median_globals'] = globals_dataframe ['density'].median()


    stat_dictionary['mean_samples'] = dataframe[dataframe.type=='s']['density'].mean()
    stat_dictionary['mean_globals'] = globals_dataframe ['density'].mean()

    stat_dictionary['std_samples'] = dataframe[dataframe.type=='s']['density'].std()
    stat_dictionary['std_globals'] = globals_dataframe ['density'].std();


    trimmed_p_samples = trimmed_bin_values_df['p_samples'].values;
    trimmed_p_globals = trimmed_bin_values_df['p_globals'].values;


    ks_d,ks_p= ks_2samp(trimmed_p_samples,trimmed_p_globals)

    stat_dictionary['ks_d'] = ks_d;
    stat_dictionary['ks_p'] = ks_p;

    return stat_dictionary, trimmed_bin_values_df;




# =============================================================================
# def write_results(aFile,anIdentifier, aPath,dataframe, globals_dataframe,population_data, min_globals, min_p):
#     
#     bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals = generate_bin_values(dataframe, globals_dataframe, population_data, min_globals);
#     wrm.write_bin_table(aFile, bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals, min_globals)
#     
#     ################################
#     # Compute and write statistics #
#     ################################
#     # - binomial test
#     # - wilcoxon
#     wrm.write_label(aFile, "Statistics")
# 
# 
#     #######################
#     # Statistics
#     #######################
# 
#     wrm.write_label(aFile, "Statistics");
# 
#     total_samples=dataframe[dataframe.type=='s']['density'].count ()
#     total_globals=globals_dataframe ['density'].count()
#     aFile.write('Total sites: '+str(total_samples)+'\n')
#     aFile.write('Total globals: '+str(total_globals)+'\n\n')
# 
#     median_samples=dataframe[dataframe.type=='s']['density'].median()
#     median_globals=globals_dataframe ['density'].median()
#     aFile.write('Median density for sites: '+str(median_samples)+'\n')
#     aFile.write('Median density for globals: '+str(median_globals)+'\n\n')
# 
# 
#     mean_samples=dataframe[dataframe.type=='s']['density'].mean()
#     mean_globals=globals_dataframe ['density'].mean()
#     aFile.write('Mean density for sites: '+str(mean_samples)+'\n')
#     aFile.write('Mean density for globals: '+str(mean_globals)+'\n\n')
# 
#     std_samples=dataframe[dataframe.type=='s']['density'].std()
#     std_globals=globals_dataframe ['density'].std()
#     aFile.write('Mean density for sites: '+str(median_samples)+'\n')
#     aFile.write('Mean density for globals: '+str(median_globals)+'\n')
# 
# 
#     #######################
#     # Test distributions are different
#     #######################
# 
#     wrm.write_label(aFile, "K-S2 Test\n")
# 
#     ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
#     aFile.write( 'KS test  for samples vs globals with full controls:'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
#     if ks_p<min_p:
#         aFile.write('The two distribitions are significantly different p<0.001'+'\n')
# 
#     ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
#     f2.write( 'KS test for p_samples vs p_globals :'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
#     if ks_p<self.min_p:
#          f2.write('The two distribitions are significantly different p<0.001'+'\n')    
#         
#     # plot graphs
#     plm.plot_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals,population_data.bin_size, anIdentifier, aPath)
#     plm.plot_cumulative_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, population_data.bin_size,median_samples,median_globals, anIdentifier, aPath)
#     plm.plot_detection_frequencies (trimmed_bin_array, trimmed_likelihood_ratios, population_data.bin_size, population_data.max_population-population_data.bin_size*2, anIdentifier, "detection_frequencies", aPath)
# 
#     # - plots targets and globals on a map
#     plm.plot_targets_on_map(dataframe, globals_dataframe, aPath, anIdentifier)
#     
# 
#     
# 
# =============================================================================
