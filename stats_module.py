from __future__ import division #This means division always gives a floating result
import plot_module as plm
import numpy as np
import pandas as pd;
from scipy.stats import linregress, ks_2samp,poisson;
from scipy.sparse import spdiags
from math import *
import sys
import cPickle as pkl
import json
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon





def compute_likelihood_model(directory,results_path, merged_dataframe, low_res=False):
 # fix parameter values for scan
 # when the low_res parameter is set to true, the system produces low_res graphs. Used for system testing and exploratory testing
    if low_res:
        lambda_v=np.arange(25,50,1) #same for eriksson and timmermann 
        eps_v=np.linspace(0,0.1,num=11,endpoint=False) #eriksson only
 #       eps_v=np.linspace(0,0.5,num=11,endpoint=False) #timmermann
        zetta_v=np.exp(np.linspace(log(1e-5),log(1e-3),num=11,endpoint=False)) #eriksson only
       # zetta_v=np.exp(np.linspace(log(1e-4),log(1e-1),num=11,endpoint=False)) #timmermann  
    else:
        lambda_v=np.arange(25,50,0.1)
        eps_v=np.linspace(0,0.1,num=101,endpoint=False)
        zetta_v=np.exp(np.linspace(log(1e-5),log(1e-3),num=101,endpoint=False))      
    rho_bins=np.arange(2,3001)# creates numbers (2...3000). Note: in the matlab this is a COLUMN VECTOR - can't see why we need 3000 when using timmermann data
    rho_bins_4_python=np.append(rho_bins,3001)  
   #     bin_width=300 #This is obviously wide for timmermann data
    bin_width=200 #This gives nicer looking graph
    bin_boundaries2_4_python=np.arange(2,3001,bin_width) 
    samples_counts=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==1],bins=rho_bins_4_python)[0] #Column  vector This is not strict translation of mathlab code. In mathlab the last bin contains 3000. In python it contains 2999-3000
    controls_counts=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==0],bins=rho_bins_4_python) [0] #Column  vector
    globals_counts=np.histogram(merged_dataframe['density'],bins=rho_bins_4_python)[0]
    control_counts2=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==0],bins=bin_boundaries2_4_python)[0]
    sample_counts2=np.histogram(merged_dataframe['density'][merged_dataframe.is_sample==1],bins=bin_boundaries2_4_python)[0]
    n_controls=np.sum(controls_counts)  
    n_samples=np.sum(samples_counts) 
    n_globals=n_controls+n_samples # not sure these are needed
    rho_bins2=bin_boundaries2_4_python[0:len(bin_boundaries2_4_python)-1]+bin_width/2 #This is EC2 not sure this is correct - depends on usage later on, needs to be adjusted to take account of nature of bins
    l_shift=n_samples*(log(float(n_samples)/float(n_controls))-1)
    n_lambda=len(lambda_v)
    n_eps=len(eps_v)
    n_zetta=len(zetta_v)
    y_acc=np.linspace(0,1e-4,num=401) #COLUMN VECTOR. In Tindbergen program he had different values that generated a lot of zeros. These in turn created problems when we had to divide by last element in yy
    acc=np.zeros((len(y_acc),len(rho_bins)))
    lnL=np.zeros((n_lambda,n_eps,n_zetta))
    sqrt_rho_bins=np.sqrt(rho_bins) #These are the values we are computing - rhobins_4_python are intervals for histogram only. In original program were inside loop. Have moved it outside
    max_LL=-float('inf')
    for i_lambda in range (0,n_lambda):
        print '.'
        my_lambda=float(lambda_v[i_lambda]) 
        bin_zeros=np.zeros(len(sqrt_rho_bins)) #in python I can't compare a vector with a scalar
        pI=np.maximum(bin_zeros,1-my_lambda/sqrt_rho_bins)  #COLUMN VECTOR
        for i_zetta in range(0,n_zetta):
            for i_eps in range (0, n_eps):
                 pObs=np.zeros(len(pI)).astype(float)                
                 pObs=zetta_v[i_zetta]*((1-eps_v[i_eps])*pI+eps_v[i_eps]) #COLUMN VECTOR (scalars * a column vector)
                 pObs_small=np.zeros(len(pObs))
                 pObs_small.fill(1e-20)
                 pObs=np.maximum(pObs,pObs_small)
                 pObs=pObs.astype(float)
                 log_samples=np.dot(samples_counts,np.log(pObs)) 
                 log_controls=np.dot(controls_counts,np.log(1-pObs)) 
                 LL=log_samples+log_controls
                 if LL>max_LL:
                     max_LL=LL
                     max_lambda=lambda_v[i_lambda]
                     max_zetta=zetta_v[i_zetta]
                     max_eps=eps_v[i_eps]
                 L=np.exp(LL-l_shift)
                 lnL[i_lambda,i_eps,i_zetta]=LL # This is different from original datastructure. Will require change of later code. I could also assign using an array op.
        #         rhs=np.floor(1+pObs/y_acc[1]).astype(int) #This is original code - yields a 1-based index. I dont like y_acc[1]- This is actually y_acc step
                 len_y_acc=np.array(len(y_acc))
                 len_y_acc.fill(len(y_acc))
                 i_acc=np.minimum(len_y_acc,np.floor(1+pObs/y_acc[1]).astype(int)) #vThis yields column vector of indexes corresponding to different values of pObs. ector length =401. Maximum value of index =400 (zero based vector).I am keeping it 1-based
                 for i in range(0,len(i_acc)):
                     x_coord=i_acc[i]-1
                     y_coord=i
                     acc[x_coord,y_coord]=acc[x_coord,y_coord]+L
    scale = (2/sqrt(3))/100 #  convert from hexagon pop size to Timmermann units, inds/100 km^2
    lambda_v = lambda_v*sqrt(scale)
    rho_bins = rho_bins*scale
    rho_bins2=rho_bins2*scale
    max_lambda=max_lambda*sqrt(scale)
    opt_threshold=max_lambda**2  #Not sure about this
    plm.plot_maximum_likelihood(acc,rho_bins,rho_bins2,y_acc, lambda_v, opt_threshold, sample_counts2, control_counts2, scale,directory,results_path)
    plm.plot_parameter_values(lnL,lambda_v, zetta_v, eps_v,directory,results_path)
    return(max_lambda, max_zetta, max_eps, opt_threshold)
    
 

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
        
        # likelihood ratio: sample_count/global_count
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

    stat_dictionary['total_samples'] = dataframe[dataframe.type=='s']['density'].sum()
    stat_dictionary['total_globals'] = globals_dataframe ['density'].sum()

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




def write_results(aFile,anIdentifier, aPath,dataframe, globals_dataframe,population_data, min_globals, min_p):
    
    bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals = generate_bin_values(dataframe, globals_dataframe, population_data, min_globals);
    wrm.write_bin_table(aFile, bin_array, sample_counts, global_counts, likelihood_ratios, p_samples, p_globals, min_globals)
    
    ################################
    # Compute and write statistics #
    ################################
    # - binomial test
    # - wilcoxon
    wrm.write_label(aFile, "Statistics")


    #######################
    # Statistics
    #######################

    wrm.write_label(aFile, "Statistics");

    total_samples=dataframe[dataframe.type=='s']['density'].sum()
    total_globals=globals_dataframe ['density'].sum()
    aFile.write('Total sites: '+str(total_samples)+'\n')
    aFile.write('Total globals: '+str(total_globals)+'\n\n')

    median_samples=dataframe[dataframe.type=='s']['density'].median()
    median_globals=globals_dataframe ['density'].median()
    aFile.write('Median density for sites: '+str(median_samples)+'\n')
    aFile.write('Median density for globals: '+str(median_globals)+'\n\n')


    mean_samples=dataframe[dataframe.type=='s']['density'].mean()
    mean_globals=globals_dataframe ['density'].mean()
    aFile.write('Mean density for sites: '+str(mean_samples)+'\n')
    aFile.write('Mean density for globals: '+str(mean_globals)+'\n\n')

    std_samples=dataframe[dataframe.type=='s']['density'].std()
    std_globals=globals_dataframe ['density'].std()
    aFile.write('Mean density for sites: '+str(median_samples)+'\n')
    aFile.write('Mean density for globals: '+str(median_globals)+'\n')


    #######################
    # Test distributions are different
    #######################

    wrm.write_label(aFile, "K-S2 Test\n")

    ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
    aFile.write( 'KS test  for samples vs globals with full controls:'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
    if ks_p<min_p:
        aFile.write('The two distribitions are significantly different p<0.001'+'\n')

    ks_d,ks_p=ks_2samp(trimmed_p_samples,trimmed_p_globals)
    f2.write( 'KS test for p_samples vs p_globals :'+str(float(ks_d))+ '   p='+str(float(ks_p))+'\n')
    if ks_p<self.min_p:
         f2.write('The two distribitions are significantly different p<0.001'+'\n')    
        
    # plot graphs
    plm.plot_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals,population_data.bin_size, anIdentifier, aPath)
    plm.plot_cumulative_p_graphs(trimmed_bin_array, trimmed_p_samples, trimmed_p_globals, population_data.bin_size,median_samples,median_globals, anIdentifier, aPath)
    plm.plot_detection_frequencies (trimmed_bin_array, trimmed_likelihood_ratios, population_data.bin_size, population_data.max_population-population_data.bin_size*2, anIdentifier, "detection_frequencies", aPath)

    # - plots targets and globals on a map
    plm.plot_targets_on_map(dataframe, globals_dataframe, aPath, anIdentifier)
    

    
