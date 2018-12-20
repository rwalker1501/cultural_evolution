import numpy as np
import math
import copy
import plot_module as plm
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import write_module as wrm
import pandas as pd;
from mpmath import mpf, exp
from classes_module import Target, PopulationData
from scipy.optimize import curve_fit
from scipy import stats;
from scipy.stats import linregress, chisquare, binom_test, levene,wilcoxon,ks_2samp
from sklearn.metrics import r2_score
from sklearn.model_selection import RepeatedKFold
from sklearn.linear_model import LogisticRegression
from operator import itemgetter




def process_dataframe(data, max_for_uninhabited):

    # print(data)
    all_sample_medians=[]
    all_control_medians=[]
    samples_growth_coefficients = []
    controls_growth_coefficients = []
    samples_gt_controls = 0
    growth_samples_gt_controls=0
    n_targets_gt_0 = 0

    max_cluster_id = data['cluster_id'].max()
    if np.isnan(max_cluster_id):
        max_cluster_id = 0


    for i in range(0,max_cluster_id+1):

        ######################
        # For each cluster.. #
        ######################
        cluster_df = data[data.cluster_id==i]
        # consider only data where densities > max_for_uninhabited
        cluster_df = cluster_df[cluster_df.density > max_for_uninhabited]


        #########################################
        # Get the median for samples and globals #
        #########################################
        # Get every sample in the cluster..
        sample_cluster_df = cluster_df[cluster_df.type == 's']
        # take the mean of these samples
        sample_median = sample_cluster_df['density'].median()

        #if the sample mean is zero we do not consider the sample or the related global
        if np.isnan(sample_median): 
            sample_median=0
            # we delete the cluster from the original dataframe
            data = data[data.cluster_id != i]
            continue
        n_targets_gt_0 += 1

        # Get every global in the cluster..
        control_cluster_df = cluster_df[cluster_df.type == 'c']
        # take the mean of these globals
        control_median = control_cluster_df['density'].median()
        if np.isnan(control_median):
            control_median=0
        if sample_median>control_median:
            samples_gt_controls += 1
        all_sample_medians.append(sample_median)
        all_control_medians.append(control_median)


        ###############################
        # Compute growth coefficients #
        ###############################

        # extract all periods and all population as arrays from the samples dataframe
        sample_times = sample_cluster_df['period'].values
        sample_populations = sample_cluster_df['density'].values

        # extract all periods and all population as arrays from the controls dataframe
        control_times = control_cluster_df['period'].values
        control_populations = control_cluster_df['density'].values

        # compute growth coefficients for samples and globals
        growth_coefficient_samples=compute_growth_coefficient(sample_times, sample_populations)
        samples_growth_coefficients.append(growth_coefficient_samples)
    
    return  samples_growth_coefficients, data

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

def get_confounder_analysis_values(keys, values):
    
    sample_bin_array = values[keys[0]][0];

    or_MHs = [];
    MH_statistics = [];
    MH_statistic_p_values = [];

    top_or_MH = [];
    bottom_or_MH = [];
    top_MH_statistics = [];
    bottom_MH_statistics = [];
    for bin in sample_bin_array:
        top_or_MH.append(0);
        bottom_or_MH.append(0);
        top_MH_statistics.append(0);
        bottom_MH_statistics.append(0);
    
    for key in keys:
        # data = [bin_array, sample_counts, global_counts, control_counts, odds_ratios, top_MHs, bottom_MHs, top_test_MHs, bottom_test_MHs]
        data = values[key];
        top_MHs = data[5];
        bottom_MHs = data[6];
        top_test_MHs = data[7];
        bottom_test_MHs = data[8];
        for j in range(0, len(top_MHs)):
            top_or_MH[j] += top_MHs[j];
            bottom_or_MH[j] += bottom_MHs[j]
            top_MH_statistics[j] += top_test_MHs[j];
            bottom_MH_statistics[j] += bottom_test_MHs[j];

    for i in range(0, len(sample_bin_array)):
        or_MH = top_or_MH[i]/bottom_or_MH[i];
        or_MHs.append(or_MH);

        MH_statistic = (top_MH_statistics[i]*top_MH_statistics[i])/bottom_MH_statistics[i]
        MH_statistics.append(MH_statistic);

        MH_statistic_p_values.append(1-stats.chi2.cdf(MH_statistic, 1));

    return or_MHs, MH_statistics, MH_statistic_p_values;






def trim_values(bins, likelihood_ratios, n_samples, n_globals, ratios, p_samples, p_globals, lower_cis, upper_cis):
    
    # trimmed_bins = bins
    # trimmed_likelihood_ratios = likelihood_ratios
    # trimmed_n_samples = n_samples
    # trimmed_n_globals = n_globals
    # trimmed_p_likelihood_ratios = p_likelihood_ratios
    # trimmed_p_samples = p_samples
    # trimmed_p_globals = p_globals
    # trimmed_lower_cis = lower_cis;
    # trimmed_upper_cis = upper_cis;

    trimmed_bins = []
    trimmed_likelihood_ratios = []
    trimmed_n_samples = []
    trimmed_n_globals = []
    trimmed_ratios = []
    trimmed_p_samples = []
    trimmed_p_globals = []
    trimmed_lower_cis = [];
    trimmed_upper_cis = [];
    for x in range(0, len(bins)):
        if not math.isnan(ratios[x]):
            trimmed_bins.append(bins[x])
            trimmed_likelihood_ratios.append(likelihood_ratios[x])
            trimmed_n_samples.append(n_samples[x])
            trimmed_n_globals.append(n_globals[x])
            trimmed_ratios.append(ratios[x])
            trimmed_p_samples.append(p_samples[x])
            trimmed_p_globals.append(p_globals[x])
            if len(lower_cis) > 0:
                trimmed_lower_cis.append(lower_cis[x]);
                trimmed_upper_cis.append(upper_cis[x]);

    return trimmed_bins, trimmed_likelihood_ratios, trimmed_n_samples, trimmed_n_globals, trimmed_ratios, trimmed_p_samples, trimmed_p_globals, trimmed_lower_cis, trimmed_upper_cis

def divide_graph(bins, samples, the_globals, p_likelihood_ratios, divider):
    left_samples_array = []
    left_globals_array = []
    left_p_array = []
    right_samples_array = []
    right_globals_array = []
    right_p_array = []
    for x in range(0, len(bins)):
        if bins[x] < divider:
            left_samples_array.append(samples[x])
            left_globals_array.append(the_globals[x])
            left_p_array.append(p_likelihood_ratios[x])
        else:
            right_samples_array.append(samples[x])
            right_globals_array.append(the_globals[x])
            right_p_array.append(p_likelihood_ratios[x])

    return left_samples_array, left_globals_array, left_p_array, right_samples_array, right_globals_array, right_p_array

def compute_divided_graph_variance(left_array, right_array):
    left_variance = np.var(left_array)
    right_variance = np.var(right_array)
    lev_stat, lev_p = levene(left_array, right_array)
    return left_variance, right_variance, lev_stat, lev_p

def compute_divided_graph_ratio(left_samples_array, left_globals_array, right_samples_array, right_globals_array):

    left_samples_sum = float(sum(left_samples_array))
    right_samples_sum = float(sum(right_samples_array))
    
    left_globals_sum = float(sum(left_globals_array))
    right_globals_sum = float(sum(right_globals_array))


    left_p_samples = left_samples_sum/(left_samples_sum+right_samples_sum)
    right_p_samples = right_samples_sum/(left_samples_sum+right_samples_sum)
    left_p_globals = left_globals_sum/(left_globals_sum+right_globals_sum)
    right_p_globals = right_globals_sum/(left_globals_sum+right_globals_sum)

    left_p_likelihood_ratio = -1
    if left_p_globals != 0:
       left_p_likelihood_ratio = left_p_samples/left_p_globals
    right_p_likelihood_ratio = -1
    if right_p_globals != 0:
        right_p_likelihood_ratio = right_p_samples/right_p_globals

    return left_p_likelihood_ratio, right_p_likelihood_ratio

def generate_p_threshold_and_binomial(p_samples, p_globals, bin_array):
    binomial = 100
    threshold = 0
    threshold_trial_count = 0
    threshold_success_count = 0
    threshold_samples = []
    threshold_globals = []

    #########################
    # For every bin value.. #
    #########################
    for i in range(1, len(bin_array)):
        success = 0
        trials = 0

        ########################
        # Count success/trials #
        ########################
        for x in range(0, i):
            trials += 1
            if p_samples[x] <= p_globals[x]:
                success += 1
    
        ###########################################
        # Compute binomial for success and trials #
        ###########################################
        temp = binom_test(success, trials, 0.5)

        ################################################
        # Save values if binomial is smaller than prev #
        ################################################
        if temp < binomial:
            binomial = temp
            threshold = bin_array[i]
            threshold_trial_count = trials
            threshold_success_count = success
            threshold_samples = p_samples[:i]
            threshold_globals = p_globals[:i]

    return binomial, threshold, threshold_success_count, threshold_trial_count, threshold_samples, threshold_globals


def lms(data_y,predicted_y):
    error=0
    for i in range(0,len(data_y)-1):
        error=error+(data_y[i]-predicted_y[i])**2
    return error 
    
def linear_fit(data_x,data_y):
    p_opt,p_cov=curve_fit(linear_model, data_x,data_y)
    return p_opt, p_cov   

# =============================================================================
# def fit_to_threshold_model(data_x,sample_counts, global_counts):
#     data_y=[]
#     for i in range(0,len(sample_counts)):
#         data_y.append(float(sample_counts[i]/global_counts[i]))
#     if len(data_x)>3:
# #p0=[1000.0,0.1,0.0002],
# # guessed threshold needs to be reset for Timmermann data
#  #       p_opt,p_cov=curve_fit(threshold_model, data_x,data_y, p0=[1400.0,0.001,0.002],bounds=([0, 0, 0]), diag=(1./data_x.mean(),1./data_y.mean()), [np.inf, 25, 1])
#         diag1=1./np.asarray(data_x).mean()
#         diag2=1./(np.asarray(data_y).mean())
#         print "diag1=",str(diag1)
#         print "diag2=",str(diag2)
#         p_opt,p_cov=curve_fit(threshold_model, data_x,data_y,p0=[1400.0,0.0000001])
#     else:
#         p_opt=[0,0]
#         p_cov=0
#     return p_opt, p_cov 
# =============================================================================

def generate_tuples_from_dataframe(merged_dataframe):
    unique_densities = np.sort(merged_dataframe.density.unique());
    zero_array = [0 for i in unique_densities];
    sample_base_dict = dict(zip(unique_densities,zero_array));
    global_base_dict = dict(zip(unique_densities,zero_array));

    merged_dataframe['count'] = 0;
    sample_dataframe = merged_dataframe[merged_dataframe.type == 's'];
    global_dataframe = merged_dataframe[merged_dataframe.type == 'g'];
    
    sample_tuples = generate_sorted_dataframe_tuples(sample_base_dict, sample_dataframe);
    global_tuples = generate_sorted_dataframe_tuples(global_base_dict, global_dataframe);

    return sample_tuples, global_tuples;

def logit_values(sample_tuples, global_tuples):

    x = [];
    y = [];
    for i in range(len(sample_tuples)):
        density = sample_tuples[i][0];
        sample_count = sample_tuples[i][1];
        global_count = global_tuples[i][1];

        # sample > 0 check, all densities
        x.append([density]);
        y.append(1 if sample_count > 0 else 0);

    x = np.array(x);
    y = np.array(y);

    logit = LogisticRegression();
    logit.fit(x,y);
    coef = logit.coef_;
    intercept = logit.intercept_;
    loss = logistic_model_2(x*coef + intercept).ravel();
    return x.ravel(), y, coef, intercept, loss;


    
def logistic_model(x):
    return 1 / (1 + np.exp(-x));

def logistic_model_2(x):
    return (np.exp(x))/(1+(np.exp(x)))

def generate_sorted_dataframe_tuples(base_dict, dataframe):
    print("Generate sorted dataframe tuples")
    dataframe['count'] = 0;
    dataframe = dataframe.groupby(['density']).count().reset_index();
    # print(dataframe);
    dataframe_dict = dict(zip(dataframe['density'], dataframe['count']));
    base_dict.update(dataframe_dict)
    sorted_dataframe_tuples = sorted(base_dict.iteritems(), key=itemgetter(0));

    return sorted_dataframe_tuples;


def generate_cumulated_densities(base_dict, dataframe):
    sorted_dataframe_tuples = generate_sorted_dataframe_tuples(base_dict, dataframe);

    cum_values = [0 for x in base_dict.keys()];
    cum_values[0] = sorted_dataframe_tuples[0][1];
    for i in range(1, len(sorted_dataframe_tuples)):
        cum_values[i] = cum_values[i-1] + sorted_dataframe_tuples[i][1];

    return cum_values;



def generate_cumulated_detection_frequency(merged_dataframe):
    unique_densities = np.sort(merged_dataframe.density.unique());
    zero_array = [0 for x in unique_densities];

    samples_base_dict = dict(zip(unique_densities,zero_array));
    globals_base_dict = dict(zip(unique_densities,zero_array));

    samples = merged_dataframe[merged_dataframe.type == 's']
    cum_samples = generate_cumulated_densities(samples_base_dict, samples)

    the_globals = merged_dataframe[merged_dataframe.type == 'g']
    cum_globals = generate_cumulated_densities(globals_base_dict, the_globals);


    print("unique_densities")
    print(unique_densities)
    print("Cumulative samples")
    print(cum_samples);
    print("Cumulative globals")
    print(cum_globals);
    cum_det_freq = [];
    for i in range(len(cum_globals)):
        cum_det_freq.append(float(cum_samples[i])/float(cum_globals[i]));

    return unique_densities, cum_samples, cum_globals, cum_det_freq;
   

def linear_model(x,beta):
# This function implements a linear model with no intercept
    results=[]
    for an_x in x: 
        result=float(an_x*beta)
        results.append(result)
    return (results)

def generate_statistics(dataframe, globals_dataframe, bin_values_df, minimum_globals):

    trimmed_bin_values_df = bin_values_df[bin_values_df.global_counts > minimum_globals];

    
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
    

    
